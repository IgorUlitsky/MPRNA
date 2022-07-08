# Ulitsky lab MPRNA counting suite
Tools for quantifying and visualizing MPRNA datasets
This suite of tools contains methods for analysis of MPRNA and MPRNA-RIP datasets desribed in:
* Lubelsky & Ulitsky Nature 2018
* Zuckerman et al. Mol Cell 2020
* Lubelsky et al. EMBO J 2021
* Ron et al. Nature Communications 2022
* von Kügelgen et al. bioRxiv 2021


The counting is done using Java code, and downstream analysis using R

## Initial counting of reads
You will need to make sure Java is installed, and the FASTA file of your tiles without adapters is available in `LIBRARY.fa`. You will need to make a list of your sample names in `names.txt`. The commannds use LSF, but can be easily adapted to other computing setups.

Count reads using the simple counting, splitting each FASTQ into 100 slices that are counted indepedently:
```cat names.txt | while read p; do for i in {0..99}; do bsub -o $p.$i.log -R rusage[mem=5000] java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.patch.AnalyzeConservedPatches count_reads LIBRARY.fa ${p}_1.fastq.gz -slice 100 $i $p.counts.$i.txt 1000000000; done; done```

Additional parameters:
* `-fq2 FILE` - if sequencing was with paired-end, use this option to provide the name of the 2nd FASTQ file
* `-adapter ADAPTER` and `-adapter2 ADAPTER` - by default adapter1 is TTGATTCGATATCCGCATGCTAGC and adapter2 is CGGCTTGCGGCCGCACTAGT, adjust if needed (adapter2 is expected at the beginning of read2
* `-min_read2_len NUMBER` - the minimal length of expected read2

This will generate `XXX.counts.YYY.txt` files for each slice, and `XXX.counts.YYY.txt.noMatch.fa.gz` files that contain unmapped sequences in a table format. 

##  Optional BLAST step 
Build FASTA files that contain the UMI in the header, and that can be used for BLAST:
```cat names.txt | while read p; do for i in {0..99}; do zcat $p.countsAB.$i.txt.noMatch.fa.gz | cut -f 1,2,7 | awk '{print ">" $1 $2 ":" $4 "\n" $3}' > $p.$i.not_matching.fa; done; done```

### Prepare a BLAST database with the tile sequences:
```makeblastdb -dbtype nucl -in NucLibAB_noMut.fa```

### Run BLAST for each slice:
```cat names.txt | while read p; do for i in {0..99}; do bsub -R rusage[mem=5000] -o $p.$i.log blastn -task blastn -word_size 12 -query $p.$i.not_matching.fa -db NucLibAB_noMut.fa -strand plus -out $p.$i.blast.txt -outfmt 7 -evalue 1e-10; done; done```

(wait for all jobs to finish)

### Gzip the FASTA files and the blast outputs
```for FILE in *not_matching.fa *.blast.txt; do bsub -o log.log gzip $FILE; done```

### Process the BLAST results:
```cat names.txt | while read p; do for i in {0..99}; do bsub -o $p.$i.log java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.FilterTilingBLASTMapping $p.$i.blast.txt.gz 30 $p.blast_counts.$i $p.$i.not_matching.fa.gz; done; done```

(wait for all jobs to finish)

### Combine the regular and the BLAST counts:
```cat names.txt | while read p; do bsub -R rusage[mem=2000] java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.patch.AnalyzeConservedPatches combine_count_files $p.counts,$p.blast_counts 0 99 $p.combined_blast.txt; done```

(wait for all jobs to finish - optional BLAST step over)

## Combine counts

```cat names.txt | while read p; do java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.patch.AnalyzeConservedPatches combine_count_files $p.counts 0 99 $p.combined.txt; done```

This will generate `XXX.combined.txt` files that contain the combined counts for each tile in each library

## Quality control
Uses the `names.txt` file and R version 3.6.0

Rscript qcTwist.R ./ names.txt combined.txt OUT_FILE.pdf 

where `combined.txt` is the suffix of your files built by combining the counts
This will generate `OUT_FILE.pdf` with some plots

## Filtering and computing ratios

The next script filters the data and computes ratios both old-fashioned/simple/pseudocout way and with DESeq2. It needs a configuration file `config.txt` which has 2 columns (sample in `samples/config.txt`), key and value:
* `rootDir` → Directory where the files sit
* `groupsFile` → name of a 4 column file: FILE (from names.txt) GROUP (e.g. “Nuc”) REPLICATE (e.g “1”) USE_IN_DESEQ (T or F). For groups use “WCE” for WCE samples and “Input” for plasmid samples (for later filtering)
* `ratiosFile` → name of a 3 column file : RATIO_NAME (name of the column in the output file) FILE1 FILE2 → this will be used to compute ‘simple’ ratios (e.g., per replicate). The ratios will be FILE1/FILE.
* `groupRatiosFile` → name of a 3 column file : COMPARISON_NAME GROUP1 GROUP2 (this will be used by DESeq2 to compare GROUP1 (case) to GROUP2 (control). Note that the ratios will be GROUP1/GROUP2.
* `suffix` → suffixes of the files, as in the QC script
* `outBase` → base of the filename for which output files will be formed
* `repeatFile` → used for annotating repeats (can be empty)
* `featuresFile` → used for annotating features (can be empty)
* `minWCEReads` → minimum average number of reads in all the samples in groups with “WCE” in their name. Filtering will only take place if there are such groups.
* `minInputReads` → minimum average number of reads in all the samples in groups with “Input” in their name. Filtering will only take place if there are such groups.
* `minRatioReads` → minimum number of reads in FILE1+FILE2 for the ratio to be not NA
* `pseudo` → pseudocount used to compute ratios
* `useUMIs` → whether UMIs (T) or reads (F) should be used
* `useReps` → whether replicate number should be used in DESeq (only use if there is >1!
* `useDESeq2` → whether DESeq2 should be used (default - T)


Run the script:
```Rscript processTwist.R config.txt```

## Generating context plots

The ratios can then be ploted accross the genes represented in the library. This script needs a configuration file that specifies the output of the previous step, the condition that needs to be plotted and the number of individual replicates. This file also has 2 columns (sample in `samples/config.txt`), key and value:
* `rootDir` → Directory where the files sit
* `finalTab` → Table that contains the ratios (output from the previous step)
* `condition` → The ratio that will be plotted
* `numOfReplicates` → The number of replicates plotted
* `outFile` → The name of the output file

Run with `Rscript ContextPlot.R CONFIG_FILE`
A sample configuration file is available at samples/contextPlot.configFile.txt

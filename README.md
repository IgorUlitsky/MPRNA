# Ulitsky lab MPRNA designing counting suite
Tools for designing, quantifying and visualizing MPRNA datasets
This suite of tools contains methods for analysis of MPRNA and MPRNA-RIP datasets desribed in:
* Lubelsky & Ulitsky Nature 2018
* Zuckerman et al. Mol Cell 2020
* Lubelsky et al. EMBO J 2021
* Ron et al. Nature Communications 2022
* von Kügelgen et al. Nature Neuroscience 2023

The counting is done using Java code, and downstream analysis using R. You will need to make sure Java is installed.

## MPRNA library design

```java -Xmx48000m -cp jar/compbioLib.jar:jar/compbio.jar:jar/picard.jar scripts.lincs.patch.AnalyzeConservedPatches design_patches OUT_FILE BED_FILE GENOME_2BIT_FILE DESCRIPTION_FILE TILE_LEN OFFSET_LEN_COVERED OFFSET_LEN_CONSERVED EXCLUDED_SEQS N_CONTROLS CUSTOM_SEQ_FILES CUSTOM_PRIMERS PATCH_FILE```

Parameters:
* ``OUT_FILE`` base name of the output file
* ``BED_FILE`` BED format file (can be gzipped) that contains transcript co-orindates (use `NONE` if not needed)
* ``GENOME_FILE`` Genome sequence file in 2bit format
* ``DESCRIPTION_FILE`` Tab-deliminated library description file (see below)
* ``TILE_LEN`` length of the designed library tiles
* ``OFFSET_LEN_COVERED`` length of offset for regular sequences
* ``OFFSET_LEN_CONSERVED`` length of offset for conserved patches (provided in the ``PATCH_FILE``)
* ``EXCLUDED_SEQS`` comma-separated list of sequences that have to be avoided in tiles (e.g., restriction enzyme cut sites)  (use `NONE` if not needed)
* ``N_CONTROLS`` number of dinucleotide-shuffled controls to design for each sequence
* ``CUSTOM_SEQ_FILES`` FASTA file of custom sequences for which tiles will be designed (separated by colons, use `NONE` if not needed)
* ``CUSTOM_PRIMERS`` pairs of primers to use for the custom sequences, separated by : (separated by colons, use `NONE` if not needed)
* ``PATCH_FILE`` BED file of conserved patches for design

The ``DESCRIPTION_FILE`` is the key file as it contains information about the different subsets of the library.
It contains the following 6 columns:
* *Subset type* - as described below
* *Basic information* - e.g., file name, described below
* *Subset name* - offsets used for tiling this subset
* *Offset length* - offsets used for tiling this subset
* *Prefix* - prefix (forward primer) appended to this sequence
* *Suffix* - prefix (typically reverse complement of the reverse primer) appended to this sequence
  
The subset types categories are available:
* ``BED`` - a name of the BED file is expected in the *Basic information*, the exonic sequences of each BED element are loaded and tiled with the offset lengths specified for this subset
* ``File`` - a name of the FASTA file is expected in the *Basic information*, which can be followed with ":" and a name of a corresponding BED FILE
* ``FileCirc`` - same as ``File``, but the sequences are circular and so back-splicing junction will also be tiled.
* ``FileMutate`` - a name of the FASTA file is expected in the *Basic information*, which can be followed with ":" and a name of a corresponding BED FILE, and in the *Offset length* the expected format is OFFSET_LEN:START:END, and positions between START and END will be systematically mutated
* ``BEDcirc`` - same as ``BED`` but for circular sequences (adding tiling also over the backspliced junction
* ``Fragment`` - a specific genomic position (``CHR:START-END`` format) is expected in *Basic information*, and it is tiled
* ``RefSeq`` - a specific transcript id is expected in *Basic information*, and it is extracted from the ``BED_FILE`` and tiled
* ``RefSeqFind`` - a specific paired ``TRANCRIPT_ID:SEQ`` is expected in *Basic information* and a pair of numbers ``UPSTREAM:TOTAL`` is expected in *Offset length*. ``TRANSCRIPT_ID`` is extracted from the ``BED_FILE``, then ``SEQ`` sequence is located in it. Then ``UPSTREAM`` bases are appended to the ``SEQ`` and additional bases are added to meet the final length of ``TOTAL``. These tiles are called *Context* tiles.
* ``RefSeqFindMutate`` - a specific paired ``TRANCRIPT_ID:SEQ`` is expected in *Basic information* and a pair of numbers ``UPSTREAM:TOTAL:OLD:NEW`` is expected in *Offset length*. ``TRANSCRIPT_ID`` is extracted from the ``BED_FILE``, then ``SEQ`` sequence is located in it. Then ``UPSTREAM`` bases are appended to the ``SEQ`` and additional bases are added to meet the final length of ``TOTAL``. Then, all instances of ``OLD`` sequence are replaced with ``NEW`` sequence.
* ``Repeat`` - a FASTA file is expected in *Basic information* and START-END-TOTAL in *Offset length*, the part of the sequence between START and END will be repeated in the tiles
* ``Delete`` - a FASTA file is expected in *Basic information* and START-END-FROM-TO in *Offset length*, the part of the sequence starting between START and END and in length between MIN_LEN and MAX_LEN will be deleted (and moved to the end of the title) in the tiles
* ``RepeatSpecific`` - a kmer sequence is expected in *Basic information* and total length of the tile in *Offset length*, the kmer will be repeated to total length
* ``MutateAndRepeatSpecific`` - as ``RepeatSpecific`` but systematic mutations of the kmer will also be produced
* ``InsertMer:KMER`` - a FASTA file is expected in *Basic information* and START-END-STEP in *Offset length*. The ``KMER`` will be inserted in each position between START and END, with the requested STEPs between them
* ``ReplaceMer:KMER`` - as ``InsertMer`` but replacing instead of inserting the ``KMER``.
* ``Mutate`` - a FASTA file is expected in ``BasicInformation`` and START-END-STEP in *Offset length*. Positions between START and END in each sequence, with the requested STEPs, will be systematically mutated
* ``MutatePairs`` - as in ``Mutate``, but *Offset length* has format START-END-EXCLUDE_START-EXCLUDE_END and all pairs of positions between START and END but exclusing EXCLUDE_START-EXCLUDE_END are mutated
* ``MutateMer:KMER_LEN`` - as in ``Mutate``, but instead of point mutations, kmers of length ``KMER_LEN`` are systematically mutated (A<->T, G<->C)
* ``MutateShuffle:N_SHUFFLES`` - as in ``Mutate``, but instead of point mutations, the region between ``START`` and ``END`` is shuffled, and ``N_SHUFFLES`` random sequences are introduced in the same place

## Initial counting of reads
You need to make sure that the FASTA file of your tiles without adapters is available in `LIBRARY.fa`. You will need to make a list of your sample names in `names.txt`. The commannds use LSF, but can be easily adapted to other computing setups.

Count reads using simple counting, splitting each FASTQ into 100 slices that are counted indepedently:
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

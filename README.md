# Ulitsky lab MPRNA counting suite
Tools for quantifying and visualizing MPRNA datasets
This suite of tools contains methods for analysis of MPRNA and MPRNA-RIP datasets desribed in:
* Lubelsky & Ulitsky Nature 2018
* Zuckerman et al. Mol Cell 2020
* Lubelsky et al. EMBO J 2021
* Ron et al. 2021

The counting is done using Java code, and downstream analysis using R

## Initial counting of reads
You will need to make sure Java is installed, and the FASTA file of your tiles without adapters is available in `LIBRARY.fa`. You will need to make a list of your sample names in `names.txt`. The commannds use LSF, but can be easily adapted to other computing setups.

Count reads using the simple counting, splitting each FASTQ into 100 slices that are counted indepedently:
```cat names.txt | while read p; do for i in {0..99}; do bsub -o $p.$i.log -R rusage[mem=5000] /home/labs/ulitsky/shared/lib/run.csh scripts.lincs.patch.AnalyzeConservedPatches count_reads LIBRARY.fa ${p}_1.fastq.gz -slice 100 $i $p.counts.$i.txt 1000000000; done; done```

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
```cat names.txt | while read p; do for i in {0..99}; do bsub -o $p.$i.log /home/labs/ulitsky/shared/lib/run.csh scripts.lincs.FilterTilingBLASTMapping $p.$i.blast.txt.gz 30 $p.blast_counts.$i $p.$i.not_matching.fa.gz; done; done```

(wait for all jobs to finish)

### Combine the regular and the BLAST counts:
```cat names.txt | while read p; do bsub -R rusage[mem=2000] /home/labs/ulitsky/shared/lib/run.csh scripts.lincs.patch.AnalyzeConservedPatches combine_count_files $p.counts,$p.blast_counts 0 99 $p.combined_blast.txt; done```

(wait for all jobs to finish - optional BLAST step over)

## Combine counts

```cat names.txt | while read p; do /home/labs/ulitsky/shared/lib/run.csh scripts.lincs.patch.AnalyzeConservedPatches combine_count_files $p.counts 0 99 $p.combined.txt; done```

This will generate `XXX.combined.txt` files that contain the combined counts for each tile in each library

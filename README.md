# LexoNSseq2015
Contains notes and custom scripts used to analyze Illumina reads for our NS-seq paper.


GCcontentInReads: contains script for calculating GC content in reads.

quadparsersuite: contains program for getting G4 location BED file, G4-CPMR, G4-start-site-CPMR, etc


See Methods and especially the Supp Methods to our paper, "Characterizing and controlling intrinsic biases of Lambda exonuclease in nascent strand sequencing reveals phasing between nucleosomes and G-quadruplex motifs around a subset of human replication origins", for more information on analyzing peaks and reads.
Get reads at the NCBI Short Read Archive (SRA): http://www.ncbi.nlm.nih.gov/sra

Search for: SRP045284

This should return the following 7 results:

1. LexoG0 Reads (MCF7)- Rep1
1 ILLUMINA (Illumina HiSeq 2000) run: 162.1M spots, 8.1G bases, 5Gb downloads
Accession: SRX669870
Select item 1239341

2. NS-seq Reads (MCF7)- Rep 3
1 ILLUMINA (Illumina HiSeq 2000) run: 136.2M spots, 6.8G bases, 4.2Gb downloads
Accession: SRX865372
Select item 1239340

3. NS-seq Reads (MCF7)- Rep 2
1 ILLUMINA (Illumina HiSeq 2000) run: 139.3M spots, 7G bases, 4.2Gb downloads
Accession: SRX865371
Select item 1239339

4. LexoG0 Reads (MCF7) - Rep3
1 ILLUMINA (Illumina HiSeq 2000) run: 174.9M spots, 8.7G bases, 5.3Gb downloads
Accession: SRX865370
Select item 1239338

5. LexoG0 Reads (MCF7) - Rep2
1 ILLUMINA (Illumina HiSeq 2000) run: 196.5M spots, 9.8G bases, 6.1Gb downloads
Accession: SRX865369
Select item 939888

6. NS-seq Reads (MCF7)- Rep 1
1 ILLUMINA (Illumina HiSeq 2000) run: 128.2M spots, 6.4G bases, 3.8Gb downloads
Accession: SRX669871
Select item 939886

7. G0gDNA Control Reads (MCF7)
1 ILLUMINA (Illumina HiSeq 2000) run: 193.6M spots, 9.7G bases, 6.1Gb downloads
Accession: SRX669869



Finding other data relevant to the paper:
1. Gap locations in hg19

2. CpG island locations in hg19

3. human rDNA sequence

4. 

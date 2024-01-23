# 1. Running fastqc 
## illumina matepair reads
getting illumina matepair reads
```bash
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR335/007/SRR3359937/SRR3359937_1.fastq.gz 
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR335/007/SRR3359937/SRR3359937_2.fastq.gz 
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR335/007/SRR3359936/SRR3359936_1.fastq.gz 
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR335/007/SRR3359936/SRR3359936_2.fastq.gz 
```

running fastqc
```bash
    fastqc -o '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC'' -f fastq '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC/SRR3359936_1_fastqc.zip'

    fastqc -o '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC'' -f fastq '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC/SRR3359936_2_fastqc.zip'

    fastqc -o '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC' -f fastq '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC/SRR3359937_1_fastqc.zip'

    fastqc -o '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC' -f fastq '../../../data/genome/01_qualityCheck/files_genomeV1/matepairReadsQC/SRR3359937_1_fastqc.zip'
```

## current reference genome
running fastqc
```bash
    fastqc -o '../../../data/genome/01_qualityCheck/files_genomeV1' -f fastq '../../../data/genome/01_qualityCheck/files_genomeV1/Asp11-DH88_14_fastqc.zip'
```

## pacBio reads
running fastqc
```bash
    fastqc -o '../../../data/genome/01_qualityCheck/files_genomeV2' -f fastq '../../../data/genome/01_qualityCheck/files_genomeV2/REF88_DNA_m84074_230730_135048_s4.hifi_reads_fastqc.zip'
```

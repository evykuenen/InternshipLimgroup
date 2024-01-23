get test data 
```bash
 wget http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R1_001.fastq.gz
 wget http://juicerawsmirror.s3.amazonaws.com/opt/juicer/work/HIC003/fastq/HIC003_S2_L001_R2_001.fastq.gz
 
 wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta
 wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.amb
 wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.ann
 wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.bwt
 wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.pac
 wget https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.fasta.sa
```

get right folder structure
```bash 
cp '../../../data/genome/04_Hi-c/reference_genome/chr19HomoSapiens.fasta references/chr19HomoSapiens.fasta'
cp '../../../data/genome/04_Hi-c/Hi-C_test_reads/HIC003_S2_L001_R1_001.fastq.gz /home/student/projects/asp-pan-genome-evy/data/genome/04_Hi-c/juicer/fastq'
gunzip '../../../data/genome/04_Hi-c/juicer/fastq/HIC003_S2_L001_R1_001.fastq.gz'
cp '../../../data/genome/04_Hi-c/Hi-C_test_reads/HIC003_S2_L001_R2_001.fastq.gz /home/student/projects/asp-pan-genome-evy/data/genome/04_Hi-c/juicer/fastq'
gunzip '../../../data/genome/04_Hi-c/juicer/fastq/HIC003_S2_L001_R2_001.fastq.gz'
```

index genome
```bash
bwa index '../../../data/genome/04_Hi-c/reference_genome/chr19HomoSapiens.fasta'
```

make input json test (json needs to have specific structure and heap size and ram_gb must be doubled if program stops)
```bash
python scripts/make_input_json_from_portal.py -a ENCSR645ZPH --enzymes HindIII --outfile '../../../data/genome/04_Hi-c/juicer/test.json'
```

run hic-pipeline
```bash
caper run '../../../data/genome/04_Hi-c/juicer/hic-pipeline/make_restriction_site_locations.wdl' -i '../../../data/genome/04_Hi-c/juicer/hic-pipeline/test_make_restriction.json'

caper run '../../../data/genome/04_Hi-c/juicer/hic-pipeline/hic.wdl' -i '../../../data/genome/04_Hi-c/juicer/test.json'
```






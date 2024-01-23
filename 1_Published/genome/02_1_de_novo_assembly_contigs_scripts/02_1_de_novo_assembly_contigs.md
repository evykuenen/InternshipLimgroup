# 1. Getting new contigs:
Getting contigs flye:
```bash
  time '../../../data/genome/02_deNovoAssembly/contigs/tools/Flye/bin/' --pacbio-hifi '../../../data/PacBio/WGS/DH88_14/HiFi_data/REF88_DNA_m84074_230730_135048_s4.hifi_reads.bam' \
      --genome-size 1.3g \
      --threads 128 \
      --no-alt-contigs \
      --scaffold \
      --out-dir '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput'
```

# 2. Getting current contigs
get agp with contigs on current refseq:
```bash
  wget --recursive -e robots=off --reject "index.html" --no-host-directories --cut-dirs=6 https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/4686/100/GCF_001876935.1_Aspof.V1/GCF_001876935.1_Aspof.V1_assembly_structure/Primary_Assembly/
```

getting all contigs in current refseq (needs to be done per chromosome so paths need to be changed in code):
```bash
  python '../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/getAllContigsOfHarkess.py'
```

combining output of script 
```bash
  cd '../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/'| cat *.fa > '../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/AllCurrentContigs.fa'
```


# 3. mapping current contigs to new contigs
activate mummer
```bash
conda activate mummer4
```

mapping contigs to refseq:
```bash
  nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq' \
  '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' -l 200 -c 100 \
  '../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/AllCurrentContigs.fa'
```

filter delta file
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.filtered.delta'
```

getting mapping info
```bash
  show-coords -qcl ../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.delta > ../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.coord
```

plot
```bash
  mummerplot -p ../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq ../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.filtered.delta
```

# 4. Getting length distribution of new contigs and current contigs
get length distribution of new contigs
```r
---
title: "lengthDistributionFlyeContigs"
output: html_document
date: "2023-09-12"
---
```{r}
library(ggplot2)
library(gridExtra)
library("scales")

# Eerste plot
data1 <- read.table("../../../data/02_deNovoAssembly/contigs/Flye_results/results_flye_V1/40-polishing/filtered_stats.txt", header = TRUE, sep = "\t", col.names = c("seq_name", "length", "coverage"))
plot1 <- ggplot(data1, aes(x = length)) +
    geom_histogram(fill = "blue", color = "black", bins = 100) +
    geom_vline(xintercept = 594446, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Length Distribution Flye contigs", x = "Length of contig", y = "Number of contigs") +
  scale_x_continuous(labels = comma)+
    theme_minimal()

# Tweede plot
data2 <- read.table("../../../data/02_deNovoAssembly/annotation_releases/annotation_releases/4686/100/GCF_001876935.1_Aspof.V1/GCF_001876935.1_Aspof.V1_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP/combined_contig_lengths.txt", header = FALSE, sep = ",", col.names = c("Contig", "Length"))
data2$Length <- as.numeric(gsub("[^0-9]", "", data2$Length))
plot2 <- ggplot(data2, aes(x = Length)) +
    geom_histogram(fill = "blue", color = "black", bins = 100) +
    labs(title = "Length Distribution old reference genome contigs", x = "Length of contig", y = "Number of contigs") +
    geom_vline(xintercept = 10042, color = "red", linetype = "dashed", size = 1) +
  scale_x_continuous(labels = comma)+
    theme_minimal()

grid.arrange(plot1, plot2, ncol = 2)
```

# 5. Getting BUSCO score of current reference genome (files must be copied to docker paths don't work)
get BUSCO score of contigs current reference genome
```bash
  cd '../../../data/genome/02_deNovoAssembly/contigs/BUSCO-docker'
  cp '../../../refseq/AsparagusCHR_V1.1.fa'  AsparagusCHR_V1.1.fa

  docker run -it --rm -v '../../../data/02_deNovoAssembly/contigs/BUSCO-docker:/home/working' -w /home/working chrishah/busco-docker run_BUSCO.py --in 'AsparagusCHR_V1.1.fa' --out 'run_genome.BUSCO' -l 'liliopsida_odb10' -m genome -f
```

# 6. Getting BUSCO score of new contigs (files must be copied to docker paths don't work)
get BUSCO score of new contigs 
```bash
  cd '../../../data/genome/02_deNovoAssembly/contigs/BUSCO-docker'
  cp '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta contigs.fasta'

  docker run -it --rm -v '../../../data/02_deNovoAssembly/contigs/BUSCO-docker:/home/working' -w /home/working chrishah/busco-docker run_BUSCO.py --in 'contigs.fasta' --out 'run_FlyeV1contigs.BUSCO' -l 'liliopsida_odb10' -m genome -f
```


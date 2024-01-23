# 1. Flye contigs ordered oriented filtered mapped to itself
mapping new refseq to new refseq for overlap
```bash
nucmer -p '../../../data/genome/03_evaluation_genome/newGenomeToNewGenome/NewToNew' \
'../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa' -l 50 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa'
```

plotten mummer:
```bash
mummerplot -p '../../../data/genome/03_evaluation_genome/newGenomeToNewGenome/NewToNew' '../../../data/genome/03_evaluation_genome/newGenomeToNewGenome/NewToNew.delta'
```

# 2. evaluate marker order with allmaps
## 2.1 index contigs
get copy of contigs
```bash
cp '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/genome/allmapsGetContigOrder/0_index/currentGenome/scaffolds.fasta'

cp '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/genome/allmapsGetContigOrder/0_index/newGenome/scaffolds.fasta'
```

index copy current genome (not my code)
```bash
./build_index_hisat2.sh scaffolds '../../../data/genome/allmapsGetContigOrderCurrentContig/0_index/scaffolds.fasta'
```

index copy new genome (not my code)
```bash
./build_index_hisat2.sh scaffolds '../../../data/genome/allmapsGetContigOrderNewContigs/0_index/scaffolds.fasta'
```

## 2.2 make files for allmaps from genetic map (not my code)
give script rights
```bash
chmod +x ./ExtractAndMapMarkersCurrentRefseq.sh
chmod +x ./ExtractAndMapMarkersNewGenome.sh
```

make vcf file with python (new genome has other name that matches with name in fasta so chr1 and current genome has asparagusV1_01 change in r17 and change input file in r6):
```bash
python '../../../data/genome/03_evaluation_genome/allmapsGetContigOrderCurrentContig/1_onemap/make_vcf.py'
```

change ./om5 r16 to own refseq (not my code)
current genome 
```bash 
./ExtractAndMapMarkersCurrentRefseq.sh '../../../data/genome/allmapsGetContigOrderCurrentContig/0_index/scaffolds' '../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/genetische_kaart_K397.map.csv' '../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/asparagusV2.vcf scaffolds'

python '../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/getFilesNextStep.py'
```

change ./om5 r16 to own refseq (not my code)
new genome 
```bash 
./ExtractAndMapMarkersNewGenome.sh '../../../data/genome/allmapsGetContigOrderNewContigs/0_index/scaffolds' '../../../data/genome/allmapsGetContigOrderNewContigs/1_onemap/genetische_kaart_K397.map.csv' '../../../data/genome/allmapsGetContigOrderNewContigs/1_onemap/asparagusV2.vcf' scaffolds

python '../../../data/genome/allmapsGetContigOrderCurrentContig/1_onemap/getFilesNextStep.py'
```

## 2.3 run allmaps

run allmaps jcvi current genome
```bash
python -m jcvi.assembly.allmaps path '../../../data/genome/allmapsGetContigOrder/3_allmaps/scaffolds_mapped_onemap-cM.bed' '../../../data/genome/allmapsGetContigOrder/0_index/scaffolds.fasta' 
```

run allmaps jcvi new genome
```bash
python -m jcvi.assembly.allmaps path '../../../data/genome/allmapsGetContigOrderNewContigs/3_allmaps/scaffolds_mapped_onemap-cM.bed' '../../../data/genome/allmapsGetContigOrderNewContigs/0_index/scaffolds.fasta'
```

# 3. evaluate MeDuSa scaffolds
mapping MeDuSa scaffolds to new refseq of Flye contigs
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsNewGenome/medusaOrderedToNewRefseq' \
'../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsNewGenome/orderedOrientedScaffoldsNew.fa'
```

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsNewGenome/medusaOrderedToNewRefseq.delta' >  '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsNewGenome/medusaOrderedToNewRefseq.filtered.delta'
```


plotting mummer:
```bash
mummerplot -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsNewGenome/medusaOrderedToNewRefseq' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsNewGenome/medusaOrderedToNewRefseq.filtered.delta'
```

# 4. evaluate telomere locations
running tidk:
```bash
tidk search  --string TTAGGG  --output asparagusV2 --dir '../../../data/genome/03_evaluation_genome/tidk' -w 100000 '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa'

tidk plot -o '../../../data/genome/03_evaluation_genome/tidk/asparagusV2' --tsv '../../../data/genome/03_evaluation_genome/tidk/asparagusV2_telomeric_repeat_windows.tsv' 
```

making graph of tidk results:
```bash
python '../../../data/genome/03_evaluation_genome/tidk/makeGraphTelomeres.py'
```

# 5. evaluate centromere locations
running tandem repeat finder:
```bash
./trf '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa' 2 5 7 80 10 50 2000 -h 
```

running trf2gff to make gff file from dat file
```bash
trf2gff -i '../../../data/genome/03_evaluation_genome/trFinder/TRF-4.09.1/src/100verplaatsingenMetUnplaced.fa.2.5.7.80.10.50.2000.dat'
```

```bash
tr -d '\r' < '../../../data/genome/03_evaluation_genome/trFinder/TRF-4.09.1/src/100verplaatsingenMetUnplaced.fa.2.5.7.80.10.50.2000.gff3' > '../../../data/genome/03_evaluation_genome/trFinder/TRF-4.09.1/src/100verplaatsingenMetUnplaced_unix.gff3'
```

index new refseq
```bash
samtools faidx '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa' 
```

get bed file with window size of 1.000.000 base pairs to make graph
```bash
awk -v window_size=1000000 'BEGIN {FS="\t"}; {for (i=0; i<$2; i+=window_size) print $1 FS i FS i+window_size}' '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa.fai' > '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100_windows.bed'
```

make file with coverage of tandemrepeats in window of 1.000.000 base pairs
```bash
bedtools coverage -a '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100_windows.bed' -b '../../../data/genome/03_evaluation_genome/trFinder/TRF-4.09.1/src/100verplaatsingenMetUnplaced_unix.gff3' > '../../../data/genome/03_evaluation_genome/trFinder/coverage.txt'
```

plot centromere coverage on chromosomes
```bash
python '../../../data/genome/03_evaluation_genome/trFinder/makeCoverageGraph.py'
```

# 6. check BUSCO score
get BUSCO score of contigs new reference genome
```bash
  docker run -it --rm \
    -v '/home/student/projects/asp-pan-genome-evy/data/genome/02_deNovoAssembly/contigs/BUSCO-docker:/home/working' \
    -w /home/working \
    chrishah/busco-docker run_BUSCO.py \
    --in '100verplaatsingenMetUnplaced.fa' \
    --out 'run_newgenome.BUSCO' \
    -l 'liliopsida_odb10' \
    -m genome -f
```

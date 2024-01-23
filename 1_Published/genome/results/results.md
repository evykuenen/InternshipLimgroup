agp file with genome placed and unplaced, mitochondrial genome and chloroplast contigs
```bash
less '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenWithUnplaced.reindexed.agp'
```

command from agp to fasta in contigs
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/getContigsInOrder.py'
```

command from agp to fasta chromosomes
```bash
agpfile='/home/student/projects/asp-pan-genome-evy/data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenWithUnplaced.reindexed.agp'
componentfasta='/home/student/projects/asp-pan-genome-evy/data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'
targetfasta='./data/results/genome/Asparagus_V2_EVY.fasta'

python -m jcvi.formats.agp build $agpfile $componentfasta $targetfasta
```

fasta file with genome placed and unplaced, mitochondrial genome and chloroplast contigs
```bash
less '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa'
```

mapping contigs flye new refseq to harkess refseq in contigs 
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs' \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/AllCurrentContigs.fa' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.fa'
```

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.delta'  > '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.filtered.delta'
  ```
  
plotting
```bash
mummerplot -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs'  '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.filtered.delta' 
```
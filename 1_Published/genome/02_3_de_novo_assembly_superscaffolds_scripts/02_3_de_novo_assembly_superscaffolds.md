# getting order of flye contigs
activate mummer
```bash
conda activate mummer4
```

mapping of harkess contigs to flye contigs:
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq' \
'../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/AllCurrentContigs.fa'
```

getting coord file
```bash
show-coords -qcl '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.coord'
```

getting order and orientation of new contigs and filter them on more then 1 percent match
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/getOrderedContigs.py'
```

convert AGP to fasta file
```bash
./agptools.py assemble -o '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.fa' '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.agp'
```

# assemble centromere regions
make coord file for translocations with python script
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegionsContigsToRefseq' \
'../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getAllContigsInOldRefseq/AllCurrentContigs.fa'
```

```bash
show-coords -qcl '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/ContigsToRefseq.coord'
```

get orientation, order, filter on coverage, change current contigs if they interrrupt new contig with 100 changes
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/getOrderHarkessContigsInCentromere.py'
```

```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres' \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/juiste_volgorde_MPDIs_100verplaatsingen.fasta' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.fa'
```

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres.delta'  > '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres.filtered.delta' 
  ```

get mapping info
```bash
show-coords -qcl '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres.delta'  > '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres.coord' 
```

plotting:
```bash
mummerplot -p '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres'  '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/100Verplaatsingen/ContigsOrderedFilteredToRefseqV2WithGroupedCentromeres.filtered.delta' 
```

# get agp without gaps from flye contigs 
get agp without gaps from flye contigs
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/makeAgpFromCoord.py'
```

# append agp with unplaced contigs and with gaps
make agp with unplaced and gaps
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/findUnplaced.py'
```

reindex agp
```bash
python -m jcvi.formats.agp reindex '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenWithUnplaced.agp'
```

make fasta of agp  
```bash
python -m jcvi.formats.agp build '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenWithUnplaced.reindexed.agp' '../../../data/genome/02_deNovoAssembly/superScaffolds/getOrderdFilteredContigs/ContigsOrderedFilteredToRefseq/ContigsOrderedFilteredToRefseqContigs.fa' '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa'
```

# make file with only contigs in right order
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/getContigsInOrder.py'
```

# make chain file
```bash
perl '../../../data/genome/03_evaluation_genome/get_chain_file/crossmap_delta_to_chain.pl' --fwd_out '../../../data/genome/03_evaluation_genome/get_chain_file/forward.chain' --rev_out '../../../data/genome/03_evaluation_genome/get_chain_file/reverse.chain' '../../../data/genome/03_evaluation_genome/sequence_old_genome/NewRefToOldRef.delta'
```

# getting bed file gaps
getting bed file gaps
```bash
python '../../../data/genome/03_evaluation_genome/get_bed_file_gaps/get_bed_gaps.py' '../../../data/genome/02_deNovoAssembly/superScaffolds/findUnplaced/100verplaatsingenMetUnplaced.fa'
```
 
# liftover
liftover of genes
```bash
./liftOver -gff  '../../../gff3/gff_annotv0.2.edit.gff3' '../../../data/genome/03_evaluation_genome/get_chain_file/forward.chain' '../../../data/genome/03_evaluation_genome/get_chain_file/liftOver/genomeEvy.gff3' '../../../data/genome/03_evaluation_genome/get_chain_file/liftOver/unmapped'
```
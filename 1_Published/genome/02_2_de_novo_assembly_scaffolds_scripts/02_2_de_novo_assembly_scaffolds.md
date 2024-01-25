 # 1. running SSPACE for scaffolds
 running SSPACE:
 ```bash
 time sspace -l '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/configFileSSPACE' -s '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'  -T 5 -p -b '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/sspaceResultsV3' -v 2>&1 | tee '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/sspaceResultsV3error.txt'  
 ```

 # 2. headers omzetten voor mummer
 headers omzetten naar >scaffoldX:
 ```bash
 awk -v num=1 '/^>/ {$0 = ">scaffold" num; num++} 1' '../../../data/02_deNovoAssembly/scaffolds/SSPACE/v3/standard_output.final.scaffolds.fasta' > '../../../02_deNovoAssembly/data/02_deNovoAssembly/scaffolds/SSPACE/v3/standard_output.final.scaffolds_andere_nummering.fasta'
 ```
 
 # 3. making AGP file
 converting contigs and scaffolds to AGP (not my script):
 ```bash
 perl '../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/scaffoldsAndContigsToAGP.pl' -f '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/standard_output.final.scaffolds_andere_nummering.fasta' > '../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/flyeSspace.agp'
 ```
 
# 4. get coord file for order scaffolds

mapping harkess contigs to sspace scaffolds:
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/NewContigsV2standard_output' \
'../../../data/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/standard_output.final.scaffolds_andere_nummering.fasta' -l 200 -c 100 \
'../../../data/02_deNovoAssembly/superScaffolds/mappenOudeContigsOpNieuweContigs/merkersChromosomenHuidigGenoom.fasta' 2>&1 | tee '../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/terminalOutputMummerV2nieuwOpOudScaffoldsError.txt'
```

making filtered delta file harkess contigs to sspace scaffolds:
```bash
delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/NewContigsV2standard_output.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/NewContigsV2standard_output.filtered.delta'
```

making coord file harkess contigs to sspace scaffolds:
```bash
show-coords -qcl '../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/NewContigsV2standard_output.filtered.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/newScaffoldOrContigsAgianstRefSeq/NewContigsV2standard_output.coord'
```

# 4. ordering and orienting scaffolds
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/OrderedScaffoldsPythonScript/getOrderedSSPACEScaffolds.py'
```

 # 4. run pyscaf
run pyscaf
```bash
python pyScaf.py -f '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'  -t 7 --log '../../../data/pyscaf/ContigsToPyscaffRefseq/Pyscaf.output' -r '../../refseq/AsparagusCHR_V1.1.fa' -o '../../../data/pyscaf/ContigsToPyscaffRefseq/output_bestanden.fasta' --dotplot png --overlap 0.0
```

get genome stats
```bash
perl '../../../data/pyscaf/assemblathon2-analysis/assemblathon_stats.pl' '../../../data/pyscaf/pyScaf/scaffolds.fa'
```

 # 5. run MeDuSa
get other reference sequences 
```bash
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz #rice
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz # zea mays
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/ananas_comosus/dna/Ananas_comosus.F153.dna.toplevel.fa.gz #ananas comosus

gunzip '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/Ananas_comosus.F153.dna.toplevel.fa.gz'
gunzip '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz'
gunzip '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz'
```

headers aanpassen zodat alles na de spatie weg is
```bash
sed 's/ .*//' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/Ananas_comosus.F153.dna.toplevel.fa.gz' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/right_headers/Ananas_right_headers.fa'
sed 's/ .*//' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/right_headers/Oryza_sativa_right_headers.fa'
sed 's/ .*//' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/right_headers/Zea_mays_right_headers.fa'
```

in MeDuSa script: netcon_mummer.py change 'cpickle' import to 'pickle' line 5 

running MeDuSa
```bash
cp '../../refseq/AsparagusCHR_V1.1.fa' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/right_headers/AsparagusCHR_V1.1.fa'

medusa -i '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' -f  '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/draftsFolder/right_headers' -d -v 
```

# 6. ordering MeDuSa scaffolds 
mapping MeDuSa scaffolds to harkess refseq
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToCurrentRefseq/medusaToOld' \
'../../../refseq/AsparagusCHR_V1.1.fa' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fastaScaffold.fasta'
``` 

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToCurrentRefseq/medusaToOld.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToCurrentRefseq/medusaToOld.filtered.delta'
```

making coord file of medusa scaffolds to harkess refseq
```bash
show-coords -qcl '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToCurrentRefseq/medusaToOld.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/medusaToCurrentRefseq/medusaToOld.coord'
```

making agp with ordered MeDuSa scaffolds
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/getOrderedContigsMedusa.py' 
```

making fasta of agp
```bash
./agptools.py assemble -o '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/orderedOrientedScaffolds.fa' '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fastaScaffold.fasta' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/orderedOrientedScaffolds.agp'
```

# 7. evaluate assembly MeDuSa
mapping fasta against harkess reference
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq' \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/juiste_volgorde_MPDIs_100verplaatsingen.fasta' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/orderedOrientedScaffolds.fa'
```

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq.filtered.delta'
  ```

get coord file 
```bash
show-coords -qcl '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq.filtered.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq.coord'
```

plotten mummer:
```bash
mummerplot -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/getOrderOrientationScaffoldsCurrentGenome/medusaOrderedToCurrentRefseq.filtered.delta'
```

# get fasta of flye contigs in medusa scaffolds
make agp of flye contigs in medusa scaffolds 
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/makeAGPMedusa.py'
```
make fasta from agp
```bash
python '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/makeFastaFromAGP.py'
```

plotten fasta to harkess reference in contigs
```bash
nucmer -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigsOrderedToCurrentRef' \
'../../../data/genome/02_deNovoAssembly/superScaffolds/getCentromereRegions/juiste_volgorde_MPDIs_100verplaatsingen.fasta' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigs.fa'
```

filter delta file
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigsOrderedToCurrentRef.delta' > '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigsOrderedToCurrentRef.filtered.delta'
```


plotting harkess reference in contigs to medusa contigs
```bash
mummerplot -p '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigsOrderedToCurrentRef' '../../../data/genome/02_deNovoAssembly/superScaffolds/medusa/makeAGP/medusaContigsOrderedToCurrentRef.filtered.delta' 
```
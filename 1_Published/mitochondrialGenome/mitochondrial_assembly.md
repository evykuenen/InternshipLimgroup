# 1. downloading mitochondrial genome NCBI
   NCBI Reference Sequence: NC_053642.1
```bash
wget https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=NC_053642.1
```

# 2. flye assembly of new PacBio reads om contigs te krijgen
Getting contigs flye:
```bash
time '../../../data/genome/02_deNovoAssembly/contigs/tools/Flye/bin/' --pacbio-hifi '../../../data/PacBio/WGS/DH88_14/HiFi_data/REF88_DNA_m84074_230730_135048_s4.hifi_reads.bam' \
    --genome-size 1.3g \
    --threads 128 \
    --no-alt-contigs \
    --scaffold \
    --out-dir '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput'
```

# 3. SSPACE assembly of contigs for AGP to get scaffolds
running SSPACE:
```bash
time sspace -l configFileSSPACE -s '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'  -T 5 -p -b sspaceResultsV3 -v 2>&1 | tee '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/sspaceResultsV3error.txt' 
```

headers to >scaffoldX because MUMmer does not accept any other format header:
```bash
awk -v num=1 '/^>/ {$0 = ">scaffold" num; num++} 1' '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/standard_output.final.scaffolds.fasta' > '../../../data/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/standard_output.final.scaffolds_andere_nummering.fasta'
```


# 4. perl script scaffoldsAndContigsToAGP.pl to make AGP from new contigs and new scaffolds
converting contigs and scaffolds to AGP:
```bash
perl '../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/scaffoldsAndContigsToAGP.pl' -f '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/genome/02_deNovoAssembly/scaffolds/SSPACE/v3GoodOutput/standard_output.final.scaffolds_andere_nummering.fasta' >  '../../../data/genome/02_deNovoAssembly/scaffolds/makingAGPfile/flyeSspace.agp'
```

# 5. MUMmer mapping all new contigs to mitochondrial reference genome
activate mummer
```bash
conda activate mummer4
```

mapping of contigs against mitochondrien:
```bash
nucmer -p allContigstomitochondrien \
'../../../data/extraDNArefseqs/mitochondrialDNAasparagus.fasta' -l 200 -c 100 \
'../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'
```

# 6. python script ./02_deNovoAssembly/superScaffolds/getContigOrderAndOrientationFinal.py from .coords file and agp with new contigs make new agp with right order and orientation.
```bash
python '../../../02_deNovoAssembly/mitochondrien/getContigOrderAndOrientationMitochondrien.py'
```

# 7. reindexing agp file with jcvi
running jcvi mitochondrial:
```bash
python -m jcvi.formats.agp reindex '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.agp'
```

# 8. AGPtools make fasta with right order and orientation from orderedContigsWithOrientationMitochondrienAllContigsNoDoubles.agp

getting scaffolds:
```bash
./agptools.py assemble -o '../../../data/mitochondrien/NoDoublesPercentageMatchFinalRef/fastaMappedContigsAllContigsNoDoublesMatchPercentage.fa' '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/mitochondrien/contigsToMitochondrien/orderedContigsWithOrientationMitochondrienAllContigsNoDoubles.agp'
```

# 9. plot scaffolds against mitochondrial genome
mapping of contigs that map against mitochondrial genome without contigs that are inside other contig:
```bash
nucmer -p mappedContigsAllContigsNoDoublesMatch \
'../../../data/extraDNArefseqs/mitochondrialDNAasparagus.fasta' -l 20 -c 65 \
'../../../data/mitochondrien/NoDoublesPercentageMatchFinalRef/fastaMappedContigsAllContigsNoDoublesMatchPercentage.fa'


```

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.delta' > '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.filtered.delta'
  ```

mapping info:
```bash
show-coords -qcl '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.delta' > '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.coord'
```

show tiling:
```bash
show-tiling '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.delta' > '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.delta.tiling'
```

plot mummer:
```bash
mummerplot -p mappedContigsAllContigsNoDoublesMatch '../../../data/mitochondrien/contigsToMitochondrien/contigsVsOldRef/mappedContigsAllContigsNoDoublesMatch.delta'
```

# 10. getting reference genome of mitochondrial contigs
```bash
python '../../../mitochondrien/getGapsInAGPMitochondrien.py'
```

# 11. reindex agp
running jcvi mitochondrial to reindex agp:
```bash
python -m jcvi.formats.agp reindex '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.agp' 
```

# 12. get new refseq fasta of mitochondrial genome
getting refseq:
```bash
./agptools.py assemble -o '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMito.fa' '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.reindexed.agp'
```

# 13. plot newrefseq against mitochondrial genome
activate mummer
```bash
conda activate mummer4
```

mapping of contigs that map against mitochondrial genome without contigs that are inside other contig:
```bash
nucmer -p NewRefSeqToOldRefseq \
'../../../data/extraDNArefseqs/mitochondrialDNAasparagus.fasta' -l 20 -c 65 \
'../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMito.fa'
```

filter delta file:
```bash
 delta-filter -r -1 -l 2000 -i 95 -o 10 '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.delta' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.filtered.delta'
  ```

mapping info:
```bash
show-coords -qcl '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.delta' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.coord'
```

show tiling:
```bash
show-tiling '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.delta' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.delta.tiling'
```
plot mummer:
```bash
mummerplot -p NewRefSeqToOldRefseq '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefSeqToOldRefseq.delta'
```
# 14 creating blast database:
```bash
conda activate blast
```

```bash
new_contigs='../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta'
makeblastdb -in $new_contigs -dbtype nucl -parse_seqids -out '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/blast_database/new_contigs'
``` 

# 15 getting sequence of refseq to find matching contigs
```bash
# fasta without \n or headers
cp '../../../data/extraDNArefseqs/mitochondrialDNAasparagus.fasta' '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/mitochondrialDNAasparagusCopy.fa'

awk '/^>/{next} {print}' '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefVsOldRef/mitochondrialDNAasparagusCopy.fa'  > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefVsOldRef/newRefSeqMitoCopywithoutHeader.fa' 

awk '/^>/{printf "\n%s\n", $0; next} {printf "%s", $0} END {print ""}' '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/NewRefVsOldRef/newRefSeqMitoCopywithoutHeader.fa' | tr -d '\n ' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMitoWithoutEnters.fa'

# Extract region from 483611 bp(last from newMitoRefseq) till 492062 bp (last of current Mito Refseq)
tail -c 8451 '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMitoWithoutEnters.fa' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/region1.fa'

# Extract region from 0 bp to 700 bp
head -c 700 '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMitoWithoutEnters.fa' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/region2.fa'

# Concatenate the two regions into a new file
echo ">Last_10kbp_mito" | cat  - '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/region1.fa' '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/region2.fa' > '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/unknown_last_10kbp.fasta'
```

blast on contigs to find last 10kbp of refseq
```bash
blastn -db '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/blast_database/new_contigs' -query '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/unknown_last_10kbp.fasta -out blast_results.out' -outfmt 6
```

mappen 10kbp tegen refseq
```bash
nucmer -p last10kbp \
'../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMito.fa' -l 10 -c 65 \
'../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/unknown_last_10kbp.fasta'
```

plot mummer:
```bash
mummerplot -p last10kbp '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/blast/last10kbp.delta'
```

# sequenced parts to contigs
mappen all contigs to primer product sequences
```bash
nucmer -p sequencedproductsToAllContigs \
'../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' -l 10 -c 65 \
'../../../data/mitochondrien/sequencedPartMitochondrien/mitoprimerproductp4p10.fasta'
```

plot mummer:
```bash
mummerplot -p sequencedproductsToAllContigs '../../../data/mitochondrien/sequencedPartMitochondrien/sequencedproductsToAllContigs.delta'
```
```bash
show-coords -qcl '../../../data/mitochondrien/sequencedPartMitochondrien/sequencedproductsToAllContigs.delta' > '../../../data/mitochondrien/sequencedPartMitochondrien/sequencedproductsToAllContigs.coord'
```

gap size changed because of pcr results
```bash
'../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.reindexed.agp'
```

running jcvi mitochondrial to reindex agp:
```bash
python -m jcvi.formats.agp reindex '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.agp' 
```

# 12. get new refseq fasta of mitochondrial genome
getting refseq:
```bash
./agptools.py assemble -o '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/newRefSeqMito.fa' '../../../data/genome/02_deNovoAssembly/contigs/Flye_results/results_flye_V1_GoodOutput/30-contigger/contigs.fasta' '../../../data/mitochondrien/contigsToMitochondrien/NewRefVsOldRef/orderedContigsWithOrientationMitochondrienRefSeqNoDoubles.reindexed.agp'
```

# 13. get statistics aboud mitochondrien from geseq results (geseqresults can be obtained by going to https://chlorobox.mpimp-golm.mpg.de/geseq.html )

```bash
python '../../../data/mitochondrien/GeseqResults/job-results-20231024143612/getStatisticsMitochondrien.py'
python '../../../data/mitochondrien/GeseqResults/job-results-20231024143612/findGenesMitochondrien.py'
```
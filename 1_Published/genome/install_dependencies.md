
# FastqC
## Quality check
```bash
conda install -c bioconda fastqc
```

# Flye
## Getting contigs
downloading flye:
```bash
git clone https://github.com/fenderglass/Flye
cd Flye
make
```

# MUMmer
## Plotting alignments
installing mummer enviorment:
```bash
conda create -n mummer4
```

activate mummer:
```bash
conda activate mummer4
```

installing mummer:
```bash
conda install -c bioconda mummer
```

# BUSCO DOCKER
## Getting BUSCO score
docker voor Busco: 
``` bash
docker pull chrishah/busco-docker
```


# SSPACE
## Getting scaffolds
installing sspace:
```bash
sudo apt install sspace
```

# Pyscaf
## Getting scaffolds
installing Pyscaf:
```bash
git clone https://github.com/lpryszcz/pyScaf.git
git clone https://github.com/ucdavis-bioinformatics/assemblathon2-analysis.git
conda install -c bioconda last
```

# medusa 
## Getting scaffolds
installing medusa in mummer4 env:
```bash
conda install -c bioconda medusa
```

# AGPtools
## Make from AGP fasta or fasta from AGP
installing agptools:
```bash
git clone https://github.com/esrice/agptools.git
cd agptools
pip install .
echo $PATH
/home/student/.vscode-server/bin/f1b07bd25dfad64b0167beb15359ae573aecd2cc/bin/remote-cli:/home/student/.local/bin:/home/student/mambaforge/bin:/home/student/mambaforge/condabin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games
(base) student@chronos:~/projects/asp-pan-genome-evy$ cp /home/student/projects/asp-pan-genome-evy/data/genome/02_deNovoAssembly/superScaffolds/tools/agptools/agp/agptools.py /home/student/.local/bin
```

# Seqkit
## validating steps
installing seqkit:
```bash
conda install -c bioconda seqkit
```

# JCVI
## Make from AGP fasta or fasta from AGP and find marker locations
installing jcvi:
```bash
conda install bioconda::jcvi
```

# TIDK
## find telomere locations
installing tidk:
```bash
conda install -c bioconda tidk
```

# tandem repeat finder
## find tandem repeats for centromere locations
```bash
wget https://github.com/Benson-Genomics-Lab/TRF/archive/refs/tags/v4.09.1.zip

./configure --prefix=$HOME/trf_install
make
make install
```

# trf2gff
## make gff3 file from trf file
installing trf2gff
```bash
pip install git+https://github.com/Adamtaranto/TRF2GFF.git
```

# bedtools
## make bed file from gff3 file
installing bedtools
```bash
conda install bioconda::bedtools
```

# samtools
## index fasta files
installing samtools
```bash
conda install bioconda::samtools
```

# juicer
## get hic files for heat maps
installing hic pipeline and testing if it runs
```bash
conda install -c conda-forge autouri
pip install caper
pip install caper --upgrade
git clone https://github.com/ENCODE-DCC/hic-pipeline
caper run hic.wdl -i tests/functional/json/test_hic.json --docker
pip install -r requirements-scripts.txt
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar
```

# phasegenomics
## evaluating quality of reads 
install phase genomics qc
```bash
git clone https://github.com/phasegenomics/hic_qc.git
cd hic_qc/
conda env create -n hic_qc --file env.yml
conda activate hic_qc
python setup.py install --user
```

# samblaster
## marking duplicates for phasegenomics
install samblaster
```bash
conda install -c bioconda samblaster
```

download liftover
```bash
wget -r -np -nH --cut-dirs=2 -R index.html* http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
```

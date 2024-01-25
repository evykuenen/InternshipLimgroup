label="testHi-C"
R1="../../../data/genome/04_Hi-c/fastq/HIC003_S2_L001_R1_001.fastq"
R2="../../../data/genome/04_Hi-c/fastq/HIC003_S2_L001_R2_001.fastq"
genome="../../../data/genome/04_Hi-c/reference_genome/chr19HomoSapiens.fasta"
cpus= 8 

# first align Hi-C reads to your genome
mkdir ${label}~phase
cd ${label}~phase

bwa index $genome
bwa mem -t $cpus -5SP $genome $R1 $R2 > '../../../data/genome/04_Hi-c/1_quality_control/test/results_phasegenomics/aligned.sam'

cat '../../../data/genome/04_Hi-c/1_quality_control/test/results_phasegenomics/aligned.sam' | samblaster > '../../../data/genome/04_Hi-c/1_quality_control/test/results_phasegenomics/tmp.sam'
samtools view -@ $cpus -S -h -b -F 2316 '../../../data/genome/04_Hi-c/1_quality_control/test/results_phasegenomics/tmp.sam' > '../../../data/genome/04_Hi-c/1_quality_control/test/results_phasegenomics/blaster.bam'

# Now generate QC report
python '../../../data/genome/04_Hi-c/1_quality_control/test/results_phasegenomics/blaster'
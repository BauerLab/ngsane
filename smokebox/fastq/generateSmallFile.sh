WGS= # not using this anymore (transcript does the same job)
TRANSCR=1

module load samtools
module load jdk
module load picard

######################################
#Transcript
######################################
if [ -n "$TRANSCR" ]; then
	samtools view -b /home/cmis/bau04c/Documents/datahome/catomocolo/RNAseqPolyA/lean/tophat/14p.tph.bam 16:27184646-27472388  > fastq/test.sam
	samtools flagstat fastq/test.sam
	java -jar ${PICARD_HOME}/SamToFastq.jar \
		INPUT=fastq/test.sam \
		FASTQ=fastq/Region/test_R1.fastq 
		SECOND_END_FASTQ=fastq/Region/test_R2.fastq \
		VALIDATION_STRINGENCY=SILENT

	wc -l fastq/Region/test_R1.fastq
	wc -l fastq/Region/test_R2.fastq

fi

######################################
#Whole Genome
######################################
if [ -n "$WGS" ]; then
	dmget ${DATASTORE}/Documents/datahome/ErrorCorrection/Human/Original/bwa/ERR091571.asd.bam
	samtools view -b -f 3 ${DATASTORE}/Documents/datahome/ErrorCorrection/Human/Original/bwa/ERR091571.asd.bam 16:27184646-27472388  > WGS/t.sam
	samtools flagstat WGS/t.sam
	samtools fixmate -r WGS/t.sam WGS/test.sam
	java -jar ${PICARD_HOME}/SamToFastq.jar \
		INPUT= WGS/test.sam \
		FASTQ= WGS/ERR091571_r1.fastq \
		SECOND_END_FASTQ= WGS/ERR091571_r2.fastq \
		VALIDATION_STRINGENCY=SILENT

	wc -l WGS/ERR091571_r1.fastq
	wc -l WGS/ERR091571_r2.fastq

	gzip WGS/ERR091571_r1.fastq WGS/ERR091571_r2.fastq

fi





######################################
#Transcript
######################################
#chipseq stuff from fabian
#./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq -n ChIPseq_CTCF_chr16 ~/research/integration/TFs/H1esc/bowtie/wgEncodeBroadHistoneH1hescCtcfStdRawDataRep1.asd.bam chr16:27184646-27472388
#./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq -n ChIPseq_H3k9me3_chr16 ~/research/integration/TFs/H1esc/bowtie/wgEncodeBroadHistoneH1hescH3k09me3StdRawDataRep1.asd.bam chr16:27184646-27472388
#./extractRegionFromBam.sh -o ~/research/Sandbox_ngd1/fastq/ChIPseq_input -n ChIPseq_Input_chr16 ~/research/integration/TFs/H1esc_control/bowtie/wgEncodeBroadHistoneH1hescControlStdRawData.asd.bam chr16:27184646-27472388
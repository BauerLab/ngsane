echo "hello world"
echo "input $1"
echo "input $2"

if [ "$2" = "1" ]
then
    echo "true"
else
    echo "false"
fi

bwa aln /clusterdata/denis/shared/resources/Homo_sapiens_assembly18.fasta /clusterdata/denis/shared/Disc1/Hannibal_FC30MEJAAXX_seqs/fastq/Pool1_index10_read1.fastq
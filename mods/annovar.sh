#!/bin/bash

# Script running annovar
# QC:
# author: Denis C. Bauer
# date: Sept.2011

echo ">>>>> Annotate variants"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0)

required:
  -k | --toolkit <path>     location of the NGSANE repository 
  -i | --input <file>       vcf file
  -o | --outdir <path>      output dir
  -r | --reference <file>   reference genome

options:
  -t | --threads <nr>       number of CPUs to use (default: 1)
"
exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
THREADS=1
ADDRECAL="" # additional commands for variant recalibrator

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the NGSANE repository
        -t | --threads )        shift; THREADS=$1 ;; # number of CPUs to use
        -i1 | --input1 )        shift; snps=$1 ;; # SNP file
        -i2 | --input2 )        shift; snps2=$1 ;; # SNP file
        -i3 | --input3 )        shift; indel=$1 ;; # INDEL file
        -o | --outdir )         shift; OUT=$1 ;; # output dir
        -r | --reference )      shift; FASTA=$1 ;; # reference genome

        -n | --name )           shift; NAME=$1 ;; # name
        -g | --refseq )         shift; REFSEQROD=$1 ;; # refseq genome
        -H | --hapmap )         shift; HAPMAPVCF=$1 ;; # hapmap
        -d | --dbsnp )          shift; DBSNPVCF=$1 ;; # dbsnp
        -K | --1kg )            shift; ONEKGVCF=$1 ;; # 1000genomes data
        -L | --region )         shift; SEQREG=$1 ;; # (optional) region of specific interest, e.g. targeted reseq
        --maxGaussians )        shift; ADDRECAL=$ADDRECAL" --maxGaussians "$1 ;; #(additional params for recal)
        --percentBadVariants )  shift; ADDRECAL=$ADDRECAL" --percentBadVariants "$1 ;; #(additional params for recal)
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done



#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

export PATH=$ANNOVAR:$PATH


for f in $snps $snps2 $indel; do

    # get basename of f
	n=${f##*/}

    echo "********* convert to annovar format"
    $ANNOVAR/convert2annovar.pl $f -format vcf4 -filter PASS -includeinfo >$OUT/${n/vcf/txt}

    echo "********* autofilter"
    $ANNOVAR/summarize_annovar.pl --buildver hg19 $OUT/${n/vcf/txt} $ANNOVAR/humandb/ --remove -outfile $OUT/${n/vcf/sum}

    echo "********* summarize over all samples"
    python ${NGSANE_BASE}/bin/formatAnnovar.py $OUT/${n/vcf/sum}.genome_summary.csv $( grep "CHROM" $f ) > $OUT/${n/vcf/genome.summary.csv}


    #clean
    rm $OUT/${n/vcf/sum}.*_summary.csv
    rm $OUT/${n/vcf/txt}
    rm $OUT/*"_filtered"
    rm $OUT/*"_dropped"
    

done



echo ">>>>> Annotate variants - FINISHED"
echo ">>>>> enddate "`date`


exit


    echo " get sample ids"
    SAMPLES=`grep CHR $f | gawk '{split($0,arr,"FORMAT"); print arr[2]}'`


    for s in $SAMPLES; do 
	echo $s
    
	echo " extract individual samples"
	java -jar $GATKJAR/GenomeAnalysisTK.jar \
	    -T SelectVariants \
	    -R $FASTA \
	    -B:variant,VCF $f \
	    -sn $s \
	    -ef \
	    -o $OUT/${n/vcf/$s.vcf}


	echo " convert to annovar format"
	$ANNOVAR/convert2annovar.pl $OUT/${n/vcf/$s.vcf} -format vcf4 -filter PASS -includeinfo >$OUT/${n/vcf/$s.txt}
     

	echo " autofilter"
	$ANNOVAR/summarize_annovar.pl --buildver hg19 $OUT/${n/vcf/$s.txt} $ANNOVAR/humandb/ --remove -outfile $OUT/${n/vcf/$s.sum}

     RESULTS=$RESULTS" "$OUT/${n/vcf/$s.sum}.genome_summary.csv
     $ANNOVAR/auto_annovar.pl --buildver hg19 -model recessive $OUT/${n/vcf/$s.txt}  $ANNOVAR/humandb/

    done

    echo "python ${NGSANE_BASE}/bin/joinAnnovar.py $RESULTS > $OUT/${n/vcf/genome.summary.csv}"

    echo " summarize over all samples"
    python ${NGSANE_BASE}/bin/joinAnnovar.py $RESULTS > $OUT/${n/vcf/genome.summary.csv}




for i in gene knowngene ; do
    $ANNOVAR/annotate_variation.pl --buildver hg19 -geneanno -dbtype $i $OUT/${n/vcf/txt} $ANNOVAR/humandb/ >$OUT/${n/vcf/}.$i.txt
done

for i in band tfbs mirna mirnatarget segdup mce44way evofold dgv omimGene gwascatalog avsift ljb_pp2 ljb_mt ljb_phylop ljb_lrt; do
    $ANNOVAR/annotate_variation.pl --buildver hg19 -regionanno -dbtype $i $OUT/${n/vcf/txt} $ANNOVAR/humandb/ >$OUT/${n/vcf/}.$i.txt
done

for i in avsift ljb_pp2 ljb_mt ljb_phylop ljb_lrt; do
    $ANNOVAR/annotate_variation.pl --buildver hg19 -filteranno -dbtype $i $OUT/${n/vcf/txt} $ANNOVAR/humandb/ >$OUT/${n/vcf/}.$i.txt
done

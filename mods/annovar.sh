#!/bin/bash -e

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
        -f | --file )           shift; INPUTFILE=$1 ;;  # input file
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir 
        -n | --name )           shift; NAME=$1 ;; # name
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file                                                  
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done



#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

for MODULE in $MODULE_ANNOVAR; do module load $MODULE; done  # save way to load modules that itself load other modules
export PATH=$PATH_ANNOVAR:$PATH
module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)


echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"

f=$INPUTFILE

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
fi
    
echo -e "\n********* $CHECKPOINT\n"



# get basename of f
n=${f##*/}

################################################################################
CHECKPOINT="convert to annovar format"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	#head -n 1000 $f >test

    command="convert2annovar.pl $f --format vcf4old --filter PASS --includeinfo >$OUTDIR/${n/vcf/txt}"
    echo $command && eval $command
    
    # mark checkpoint
    if [ -e $OUTDIR/${n/vcf/txt} ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="autofilter"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

    #summarize_annovar.pl --buildver $LIBRARY $OUTDIR/${n/vcf/txt} $DATABASE --remove -outfile $OUTDIR/${n/vcf/sum}
	command="table_annovar.pl $OUTDIR/${n/vcf/txt} $DATABASE --buildver $LIBRARY --protocol refGene,knownGene,ensGene,wgEncodeGencodeManualV4,gerp++elem,phastConsElements46way,genomicSuperDups,esp6500si_all,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,cg46,cosmic64,snp129,snp132,snp138,avsift,ljb2_all --operation g,g,g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f --outfile $OUTDIR/${n/vcf/sum} --remove --otherinfo"
	#command="table_annovar.pl $OUTDIR/${n/vcf/txt} $DATABASE --buildver $LIBRARY --protocol refGene,genomicSuperDups,cosmic64 --operation g,r,f --outfile $OUTDIR/${n/vcf/sum} --remove --otherinfo"
	echo $command && eval $command

    # mark checkpoint
    if [ -e $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.txt" ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

################################################################################
CHECKPOINT="summarize over all samples"

if [[ -n "$RECOVERFROM" ]] && [[ $(grep -P "^\*{9} $CHECKPOINT" $RECOVERFROM | wc -l ) -gt 0 ]] ; then
    echo "::::::::: passed $CHECKPOINT"
else 

	# add the genotype information to the sample and make it then an csv file (cannot do it in the previous step because the genotype info is tab separated...)
    REP=$(grep "CHROM" $f) # | cut -f 10- | sed 's/\t/,/g')
    sed "s/Otherinfo/$REP/" $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.txt" |	gawk -F \\t -v OFS="\",\"" '{$1=$1; print "\""$0"\""}' > $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.csv"

# mark checkpoint
    if [ -e $OUTDIR/${n/vcf/sum}.$LIBRARY"_multianno.csv" ];then echo -e "\n********* $CHECKPOINT\n"; unset RECOVERFROM; else echo "[ERROR] checkpoint failed: $CHECKPOINT"; exit 1; fi
fi

#clean
rm $OUTDIR/${n/vcf/txt}

echo ">>>>> Annotate variants - FINISHED"
echo ">>>>> enddate "`date`



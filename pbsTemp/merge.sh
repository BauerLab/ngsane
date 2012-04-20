#!/bin/bash

# merges bam or vcf files, run directly or qsub 
# note merged bam files comtain aligned reads only.
# author: Denis C. Bauer
# date: Dec.2010

######################################
# merge all bam files from a sample
#####################################

export PERL5LIB=${PERL5LIB}:/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/perl
export PATH=$PATH:/clusterdata/hiseq_apps/bin/freeze001/tabix-0.2.3

#INPUTS
HISEQINF=$1   # location of the HiSeqInf repository
FILES=$2      # file with paths to bamfiles
OUT=$3        # output dir
NAME=$4       # output name
TASK=$5       # either "bam" or "vcf"
QOUT=$6       # either empty or a directory for the qsub output
WAIT=$7       # wait

PRIORITY="-l hp=TRUE"

#PROGRAMS
. $HISEQINF/pbsTemp/header.sh

VCFTOOLS=/clusterdata/hiseq_apps/bin/devel/vcftools_0.1.5/bin


echo ">>>>> call merge $TASK"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> merge.sh $HISEQINF $FILES $OUT $NAME $TASK $QOUT $WAIT $PRIORITY"

if [ $TASK = "bam" ]; then
    echo "merging bam "$( less $FILES )

    #delete old merged files
    if [ -e $OUT/$NAME ]; then
	rm $OUT/$NAME
    fi

    TMP=$OUT/${NAME/.bam/mergeBams.tmp}

    echo "********* merge bam"
    echo "java -Xmx10g -jar $PICARD/MergeSamFiles.jar OUTPUT=$OUT/$NAME \\" >$TMP
    for f in $( less $FILES ); do echo "INPUT=$f \\" >> $TMP; done
    echo "USE_THREADING=true" >> $TMP
    echo "$SAMTOOLS index $OUT/$NAME" >> $TMP
    echo "$SAMTOOLS flagstat $OUT/$NAME > $OUT/$NAME.stats" >> $TMP
    echo "rm $TMP" >>$TMP
    chmod -u=rwx $TMP

    # if qsub 
    if [ -n "$QOUT" ]; then

	if [ -e $QOUT/${NAME/bam/mergeBams}.out ]; then rm $QOUT/${NAME/bam/mergeBams}.out; fi

	qsub $PRIORITY -b y -cwd -j y -o $QOUT/${NAME/bam/mergeBams}.out -N ${NAME/bam/mergeBams} -hold_jid=$WAIT \
	    -pe mpich 2 $TMP
    # if normal commandline running
    else
	$TMP
    fi

elif [ $TASK = "vcf" ]; then

    echo "********* merge vcf"
    if [ -e merge$TASK.tmp ]; then rm merge$TASK.tmp; fi
    if [ -e $OUT/mergeRun.tmp ]; then rm $OUT/mergeRun.tmp; fi
    
    # zip
    echo "********* zip"
    for f in $( less $FILES ); do
	n=$(echo $f|sed 's/\///g')
	#echo "set path = ( $SAMUTILS $path )" >$OUT/zip$n.tmp
	#echo "export PATH=$PATH:$SAMUTILS">$OUT/zip$n.tmp
	echo "$SAMUTILS/bgzip -c $f >$f.gz" >> $OUT/zip$n.tmp
	echo "$SAMUTILS/tabix -p vcf -f $f.gz" >> $OUT/zip$n.tmp
	echo "rm $OUT/zip$n.tmp" >>$OUT/zip$n.tmp
	chmod -u=rwx $OUT/zip$n.tmp
	# qsub or run directly
	if [ -n "$QOUT" ]; then
	    qsub $PRIORITY -b y -cwd -j y -o $QOUT/zip$n.out -N zip$n -hold_jid $WAIT $OUT/zip$n.tmp
	else
	    $OUT/zip$n.tmp
	fi
	echo "$f.gz" >> merge$TASK.tmp
    done

    
    # merge by qsub or run directly
    echo "********* merge"
    #echo "export PERL5LIB=$PERL5LIB:/clusterdata/hiseq_apps/bin/freeze001/VCFtools_0.1.3.2/perl">$OUT/mergeRun.tmp
    #echo "set path = ( $SAMUTILS $path )" >$OUT/mergeRun.tmp
    echo "export PATH=$PATH:$SAMUTILS">$OUT/mergeRun.tmp
    echo "$VCFTOOLS/vcf-merge "$(less merge$TASK.tmp)" --silent >$OUT/$NAME">>$OUT/mergeRun.tmp
    echo "rm $OUT/mergeRun.tmp" >>$OUT/mergeRun.tmp
    chmod -u=rwx $OUT/mergeRun.tmp
    if [ -n "$QOUT" ]; then
	if [ -e $QOUT/merge$TASK ]; then rm $QOUT/merge$TASK; fi
	qsub $PRIORITY -b y -cwd -j y -o $QOUT/merge$TASK -N merge$TASK -hold_jid zip* \
	    $OUT/mergeRun.tmp
    else
        $VCFTOOLS/vcf-merge $(cat merge$TASK.tmp) --silent >$OUT/$NAME
    fi

    
    rm merge$TASK.tmp

fi
#  
# Navigates through an Illumina Runs folder and performs bclConversion
# as well as Gerald transformation (and Eland aligning). It expects to
# find a config.txt in the root of the Run folder, specifying what to
# do during the Gerald run
#
# INPUT: 
# runs folder, e.g. /illumina/Runs/101105_SN417_0109_B80C1MABXX
# date, e.g. 26-11-2010 (optional) if the bclConversion was on a
#       different day
#
# author: Denis C. Bauer
# date: Nov.2010


function usage {
    echo "MakeFastQ.sh <Run dir> <check files> <run bclConverter> <run Gerald> <Log dir> [<BaseCall date>]"
    exit
}


if [ ! $# -gt 3 ]; then usage ; fi

#DEFAULTS
fastq=0
align=0
bclConv=0
olb=0
gerald=0
check=0
INIBASECALLS=""
OVERWRITE=""
THREADS=24
MULTIP="1"

#PRIORITY="-l hp=TRUE"

MYPATH="/usr/lib64:/usr/lib:/usr/bin:$PATH:$OLB:$CASAVA"

#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -w | --workDir )        shift; WD=$1 ;; # working dir
        -c | --check )          shift; check=$1 ;; # check all files are there
	-f | --fastq )          shift; fastq=$1 ;; #fastq conversion
	-s | --singleplex )     MULTIP="" ;; #not multiplexed
        -b | --bclConv )        shift; bclConv=$1 ;; # bclConversion
        -g | --gerald )         shift; gerald=$1 ;; # gerald
        -o | --olb )            shift; olb=$1 ;; # olb
	-l | --log )            shift; LOG=$1 ;; # logdir
	-d | --data )           shift; BCLDATE=$1 ;; # certain data to pass on to gerald
	-O | --overwrite )      shift; OVERWRITE="1";;
	--basecalls )           shift; INIBASECALLS=$1 ;; # certain folder for the to pass on to gerald
        -h | --help )           usage ;;
        * )                     usage
    esac
    shift
done



QSUB="/opt/gridengine/bin/lx26-amd64/qsub"
OLB="/clusterdata/hiseq_apps/bin/freeze001/OLB-1.9.3/bin"
#CASAVA="/clusterdata/hiseq_apps/bin/freeze001/CASAVA_1.7.1a5/bin/"
#CASAVA="/clusterdata/hiseq_apps/bin/freeze001/CASAVA_v1.8.0a16-build/bin/"
#CASAVA="/clusterdata/hiseq_apps/bin/freeze001/CASAVA_v1.8.0-build/bin/"
CASAVA="/clusterdata/hiseq_apps/bin/freeze001/CASAVA_v1.8.1/build/bin/"
SCRIPTS="/clusterdata/hiseq_apps/hiSeqInf/"

ID=`whoami`
#DATE=`date "+%d-%m-%Y"`
DATE=`date "+%Y-%m-%d"`
APP=$DATE"_"$ID
# if part of the analysis was done on another day
if [ -n "$BCLDATE" ]; then APP=$BCLDATE"_"$ID; fi
# test if log dir is there
if [ ! -d $WD/$LOG ]; then mkdir $WD/$LOG; fi

n=`basename $WD`

echo `date`
echo "$ID runs MakeFastQ in $WD check:$check fastq:$fastq align:$align bclConv:$bclConv gerald:$gerald and olb:$olb $OVERWRITE"

############################################
#   check all files are there
############################################
if [ "$check" -eq "1" ]; then
    $SCRIPTS/pbsTemp/checkAllHiseqFiles.sh $WD > $WD/$LOG/fileCheck.out

    ERROR=`grep -c "Error:" $WD/$LOG/fileCheck.out`
    if [ ! $ERROR -eq "0" ]; then
	echo "Error: there are files missing and I cannot recover them see $WD/$LOG/fileCheck.out"
	exit -1
    fi

    MISSING=`grep -c "Missing:" $WD/$LOG/fileCheck.out`
    if [ ! $MISSING -eq "0" ]; then
	olb="1"
	echo "Warning: Missig files we need to deviate from normal pipeline olb:$olb"
    fi
fi


############################################
#   Basecalling
############################################
if [ "$bclConv" -eq "1" ]
then
    echo "do conversion"

    BCLPREPOUT=$WD/$LOG/bclPrep.out
    BCLOUT=$WD/$LOG/bcl.out

    if [ -n "$OVERWRITE" ]; then
	#overwrite
	echo "overwrite"
	if [ -e $BCLPREPOUT ]; then rm $BCLPREPOUT; fi
	if [ -e $BCLOUT ]; then rm $BCLOUT; fi
    else
	# new number of dir and log files
	if [ -d $WD/Data/Intensities/BaseCalls_$APP ]; then
	    let COUNTER=$(ls $WD/Data/Intensities/BaseCalls_$APP* | grep ":" | wc -l)+1
	    COUNTER="."$COUNTER
	    BASECALLOUT=$WD/Data/Intensities/BaseCalls_$APP$COUNTER
	    BCLPREPOUT=$BCLPREPOUT$COUNTER
	    BCLOUT=$BCLOUT$COUNTER
	fi
    fi

    BASECALLOUT=$WD/Data/Intensities/BaseCalls_$APP$COUNTER
    echo $BASECALLOUT

    $OLB/setupBclToQseq.py -b $WD/Data/Intensities/BaseCalls/ \
	-o $BASECALLOUT >$BCLPREPOUT 2>&1

    cd $BASECALLOUT
#    if [ "$olb" -eq "1" ]; then
#	python $SCRIPTS/bin/modifyConfig.py $WD/Data/Intensities/BaseCalls_$APP/config.xml --set $WD/$LOG/fileCheck.out
#    fi

#    $QSUB -l hp=TRUE -j y -o $BCLOUT -cwd -N bclConv"_"$n -pe mpich 16 -b y \
#	make -j 16

    $QSUB -l hp=TRUE -j y -o $BCLOUT -cwd -N bclConv"_"$n -pe make 16 -b y -v $MYPATH \
	qmake -inherit -recursive --


    cd $WD
fi




############################################
#   Fastq conversion and de multiplexing
############################################
if [ "$fastq" -eq "1" ]
then
    echo "do fastq conversion"

    FASTQPREPOUT=$WD/$LOG/fastqPrep.out
    FASTQOUT=$WD/$LOG/fastq.out
    UNALIGNEDOUT=$WD/Data/Unaligned_$APP$COUNTER

    if [ -n "$OVERWRITE" ]; then
	#overwrite
	echo "overwrite"
	if [ -e $FASTQPREPOUT ]; then rm $FASTQPREPOUT; fi
	if [ -e $FASTQOUT ]; then rm $FASTQOUT; fi
    else
	# new number of dir and log files
	if [ -d $WD/Data/Unaligned_$APP ]; then
	    let COUNTER=$(ls $WD/Data/Unaligned_$APP* | grep ":" | wc -l)+1
	    COUNTER="."$COUNTER
	    UNALIGNEDOUT=$WD/Data/Unaligned_$APP$COUNTER
	    FASTQPREPOUT=$FASTQPREPOUT$COUNTER
	    FASTQOUT=$FASTQOUT$COUNTER
	fi
    fi

    echo $UNALIGNEDOUT

    SAMPLESHEET="--sample-sheet $WD/SampleSheet.csv"

    # it if is multiplexed use the special mask with the two ns after the index
    if [ -n "$MULTIP" ]; then
	#MASK=" --use-bases-mask y100,I6nn,Y100"
        MASK=" --use-bases-mask y101,I6n,Y101"
	echo "$MASK"
    fi

    #for getting the full sequence irrespective of sample sheet info and multiplexing
    #SAMPLESHEET=""
    #MASK="--use-bases-mask n*,Y6n,n*"

    echo "$SAMPLESHEET"
    echo "$MASK"
    
     $CASAVA/configureBclToFastq.pl --input-dir $WD/Data/Intensities/BaseCalls/ $MASK \
	-output-dir $UNALIGNEDOUT $SAMPLESHEET >$FASTQPREPOUT 2>&1

    cd $UNALIGNEDOUT

    $QSUB $PRIORITY -j y -o $FASTQOUT -cwd -N "fastq_"$n -pe mpich 16 -b y -l vf=1G -v $MYPATH \
	make -j 16

    #$QSUB -l hp=TRUE -j y -o $FASTQOUT -cwd -N "fastq_"$n -pe make 16 -b y -v $MYPATH \
#	qmake -cwd -v $MYPATH -inherit -- all


    cd $WD

    echo "echo finished > $WD/$LOG/done.txt" > $WD/$LOG/writefinMsg$n.tmp
    echo "rm $WD/$LOG/writefinMsg$n.tmp" >> $WD/$LOG/writefinMsg$n.tmp
    chmod -u=wrx $WD/$LOG/writefinMsg$n.tmp

    $QSUB $PRIORITY -j y -o "doneFlag" -cwd -N finished"_"$n -hold_jid "fastq_"$n -b y -l vf=1G \
	$WD/$LOG/writefinMsg$n.tmp


fi



############################################
#   Alignment
############################################
if [ "$align" -eq "1" ]; then
    echo "make aligment"

    ALIGNPREPOUT=$WD/$LOG/alignPrep.out
    ALIGNOUT=$WD/$LOG/align.out

    if [ -e $WD/$LOG/done.txt ]; then rm $WD/$LOG/done.txt; fi

    if [ -n "$OVERWRITE" ]; then
	#overwrite
	echo "overwrite"
	if [ -e $ALIGNPREPOUT ]; then rm $ALIGNPREPOUT; fi
	if [ -e $ALIGNOUT ]; then rm $ALIGNQOUT; fi
	ALIGNEDOUT=$WD/Data/Aligned_$APP$COUNTER
    else
	# new number of dir and log files
	if [ -d $WD/Data/Aligned_$APP ]; then
	    let COUNTER=$(ls $WD/Data/Aligned_$APP* | grep ":" | wc -l)+1
	    COUNTER="."$COUNTER
	    ALIGNEDOUT=$WD/Data/Aligned_$APP$COUNTER
	    ALIGNPREPOUT=$ALIGNPREPOUT$COUNTER
	    ALIGNOUT=$ALIGNOUT$COUNTER
	fi
    fi

    echo $ALIGNEDOUT

    # if not the offline basecaller had to be used to create the bcl conversion
    if [ "$olb" -eq "1" ]; then 
	UNALIGNED=$WD/Data/Intensities/`ls $WD/Data/Intensities/ | grep "Bustard" | tail -n1`
    else
	UNALIGNED=$WD/Data/Intensities/`ls $WD/Data/* | grep "Unaligned_" | tail -n1`
    fi
    
    if [ -n "$INIBASECALLS" ]; then UNALIGNED=$INIBASECALLS; fi

    echo $UNALIGNED

    # wait for olb and bcl conf
    HOLD=" -hold_jid fastq_"$n
    if [ $olb -eq "1" ]; then
	HOLD=$HOLD",move_"$n",olb_"$n
    fi

    $CASAVA/configureAlignment.pl $WD/config.txt --EXPT_DIR $UNALIGNED \
	--OUT_DIR $ALIGNOUT --make >$ALIGNPREPOUT 2>&1 

    cd $ALIGNOUT
    $QSUB $PRIORITY -j y -o $ALIGNOUT -cwd -N "align_"$n $HOLD \
	-pe mpich $THREADS -b y \
	make -j $THREADS ALIGN=yes
    cd $WD

    echo "echo finished > $WD/$LOG/done.txt" > $WD/$LOG/writefinMsg$n.tmp
    echo "rm $WD/$LOG/writefinMsg$n.tmp" >> $WD/$LOG/writefinMsg$n.tmp
    chmod -u=wrx $WD/$LOG/writefinMsg$n.tmp

    $QSUB $PRIORITY -j y -o $GERALDOUT -cwd -N finished"_"$n -hold_jid "*"gerald"_"$n -b y \
	$WD/$LOG/writefinMsg$n.tmp
fi





############################################
#   OLB --tiles=s_3_0044
############################################
if [ "$olb" -eq "1" ]
then

    if [ -e $WD/$LOG/fileCheck.out ]; then
	RUN=`tail -n 1 $WD/$LOG/fileCheck.out | cut -b 3-`
	echo "do OLB for $RUN"
	LIMIT="--tiles=$RUN"
    else
	echo "do OLB for all"
	LIMIT="--tiles=s_1,s_2,s_3,s_4,s_5,s_6,s_7,s_8"
    fi

    #let COUNTER=$(ls MyLogs/olb.*.out | wc -l)+1
    #OLBOUT=$WD/$LOG/"olb."$COUNTER".out"
    #OLBPREPOUT=$WD/$LOG/"olbPrep."$COUNTER".out"
    #echo $OLBOUT

    if [ -e $WD/$LOG/olbPrep.out ]; then rm $WD/$LOG/olbPrep.out; fi
    if [ -e $WD/$LOG/olb.out ]; then rm $WD/$LOG/olb.out; fi

    $OLB/bustard.py --CIF $LIMIT $WD/Data/Intensities/ --with-qseq --make >$OLBPREPOUT 2>&1
    BUSTARDDIR=`ls $WD/Data/Intensities/ | grep "Bustard" | tail -n1`

    echo $BUSTARDDIR

    cd $WD/Data/Intensities/$BUSTARDDIR
   
     #$QSUB -l hp=TRUE -j y -o $OLBOUT -cwd -N olb"_"$n -pe mpich 96 -b y \
    #    make -j 96

     $QSUB -l hp=TRUE -j y -o $OLBOUT -cwd -N olb"_"$n -pe make 1-$THREADS -b y -v $MYPATH \
	  qmake -inherit -recursive --

    cd $WD

    if [ -e $WD/$LOG/fileCheck.out ]; then
	$QSUB -l hp=TRUE -j y -o $OLBOUT -cwd -N move"_"$n -b y -hold_jid olb"_"$n\
	    cp $WD/Data/Intensities/$BUSTARDDIR/*_qseq.txt $WD/Data/Intensities/BaseCalls_$APP/
    fi

fi


############################################
#   Gerald
############################################
if [ "$gerald" -gt "0" ]; then
    echo "run Gerald"

    GERALDPREPOUT=$WD/$LOG/geraldPrep.out
    GERALDOUT=$WD/$LOG/gerald.out

    if [ -e $GERALDPREPOUT ]; then rm $GERALDPREPOUT; fi
    if [ -e $GERALDOUT ]; then rm $GERALDOUT; fi
    if [ -e $WD/$LOG/done.txt ]; then rm $WD/$LOG/done.txt; fi
  
    # if not the offline basecaller had to be used to create the bcl conversion
    if [ "$olb" -eq "1" ]; then 
	BASECALLS=$WD/Data/Intensities/`ls $WD/Data/Intensities/ | grep "Bustard" | tail -n1`
    else
	BASECALLS=$WD/Data/Intensities/`ls $WD/Data/Intensities/ | grep "BaseCalls_" | tail -n1`
    fi
    
    #BASECALLS=$WD/Data/Intensities/Bustard1.9.0_01-04-2011_denis.2/BaseCalls_$APP
    if [ -n "$INIBASECALLS" ]; then BASECALLS=INIBASECALLS; fi

    echo $BASECALLS

    # wait for olb and bcl conf
    HOLD=" -hold_jid bclConv_"$n
    if [ $olb -eq "1" ]; then
	HOLD=$HOLD",move_"$n",olb_"$n
    fi

    if [ -e $WD/SampleSheet.csv ]; then
	echo "demultiplexing"
	
	$CASAVA/demultiplex.pl --input-dir $BASECALLS \
	    --output-dir $BASECALLS/Demult_$APP \
	    --sample-sheet $WD/SampleSheet.csv \
	    --alignment-config $WD/config.template.txt >$GERALDPREPOUT 2>&1
	echo $BASECALLS/Demult_$APP
	cd $BASECALLS/Demult_$APP
	$QSUB -l hp=TRUE -j y -o $GERALDOUT -cwd -N "dem_gerald_"$n $HOLD \
	    -pe mpich $THREADS -b y \
	    make -j $THREADS ALIGN=yes
	#$QSUB -l hp=TRUE -j y -o $GERALDOUT -cwd -N "dem_gerald_"$n $HOLD \
	#    -pe make 1-$THREADS -b y -v $MYPATH \
	#    qmake -inherit -recursive -- ALIGN=yes

	cd $WD

    else
	echo "normal gerald"
	$CASAVA/GERALD.pl $WD/config.txt --EXPT_DIR $BASECALLS \
	    --make >$GERALDPREPOUT 2>&1
	GERALDDIR=`ls $BASECALLS | grep "GERALD" | tail -n1`
	cd $BASECALLS/$GERALDDIR
	echo $BASECALLS/$GERALDDIR
	$QSUB -l hp=TRUE -j y -o $GERALDOUT -cwd -N gerald"_"$n $HOLD \
	    -pe mpich $THREADS -b y \
	    make -j $THREADS
	#$QSUB -l hp=TRUE -j y -o $GERALDOUT -cwd -N "norm_gerald_"$n $HOLD \
	#    -pe mpich $THREADS -b y -v $MYPATH \
	#    qmake -inherit -recursive --
	cd $WD

    fi

    echo "echo finished > $WD/$LOG/done.txt" > $WD/$LOG/writefinMsg$n.tmp
    echo "rm $WD/$LOG/writefinMsg$n.tmp" >> $WD/$LOG/writefinMsg$n.tmp
    chmod -u=wrx $WD/$LOG/writefinMsg$n.tmp

    $QSUB -l hp=TRUE -j y -o $GERALDOUT -cwd -N finished"_"$n -hold_jid "*"gerald"_"$n -b y \
	$WD/$LOG/writefinMsg$n.tmp
fi

#!/bin/bash

# Overall trigger pipeline
# author: Denis C. Bauer
# date: Nov.2010

function usage {
echo -e "usage: $(basename $0) CONFIG [TASK]

Script interpreting the config.txt file in the project directory and submitting
tasks to the queue.

required:
  CONFIG     config.txt file specifying what needs to be done and
               where the resources are located

options for TASK:
  empty      start dry-run: make dirs, delete old files, print what will be done
  armed      submit tasks to the queue
  verify     check the pbs logfiles for errors 
  html       check the pbs logfiles for errors and and make summary HTML page
  clean      clean up
"
exit
}

if [ ! $# -gt 0 ]; then usage ; fi

CONFIG=$1
ADDITIONALTASK=$2

# get all the specs defined in the config  (note both configs are necessary)
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

#PRIORITY="-l hp=TRUE"
#PRIORITY="-p 0.9"

############################################
#   Verify
############################################
if [ -n "$ADDITIONALTASK" ]; then
    if [ "$ADDITIONALTASK" = "verify" ]; then
	echo ">>>>>>>>>> $ADDITIONALTASK"
	if [ -n "$RUNMAPPINGBWA" ]; then $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/bwa.sh $QOUT/$TASKBWA; fi
	if [ -n "$recalibrateQualScore" ]; then $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/reCalAln.sh $QOUT/$TASKRCA; fi
	if [ -n "$GATKcallIndels" ]; then $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/gatkIndel.sh $QOUT/$TASKIND; fi
	if [ -n "$GATKcallSNPS" ]; then $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/gatkSNPs.sh $QOUT/$TASKSNP; fi
	if [ -n "$RUNTOPHATCUFF" ]; then $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/tophatcuff.sh $QOUT/$TASKTOPHAT; fi
	if [ -n "$RUNCUFFDIFF" ]; then $HISEQINF/pbsTemp/QC.sh $HISEQINF/pbsTemp/cuffdiff.sh $QOUT/$TASKCUFFDIFF; fi
	exit
    elif [ "$ADDITIONALTASK" = "html" ]; then
	echo ">>>>>>>>>> $ADDITIONALTASK"
	$HISEQINF/pbsTemp/makeSummary.sh $HISEQINF $CONFIG
	exit
    elif [ "$ADDITIONALTASK" = "clean" ]; then
	echo ">>>>>>>>>> $ADDITIONALTASK"
	exit
    elif [ "$ADDITIONALTASK" = "armed" ]; then
	echo ">>>>>>>>>> $ADDITIONALTASK"
	ARMED="--armed"
    elif [ "$ADDITIONALTASK" = "keep" ]; then
        echo ">>>>>>>>>> $ADDITIONALTASK"
        ARMED="--keep"
    elif [ "$ADDITIONALTASK" = "direct" ]; then
        echo ">>>>>>>>>> $ADDITIONALTASK"
        ARMED="--direct"
    elif [ "$ADDITIONALTASK" = "first" ]; then
        echo ">>>>>>>>>> $ADDITIONALTASK"
        ARMED="--first --armed"
    elif [ "$ADDITIONALTASK" = "postonly" ]; then
        echo ">>>>>>>>>> $ADDITIONALTASK"
        ARMED="--postonly"
    else
	echo -e "don't understand "$ADDITIONALTASK" \nI understand only \"verify\" and \"clean\""
	exit -1
    fi
fi

# ensure out directory is there 
for dir in ${DIR[@]}; do
  if [ ! -d $OUT/$dir ]; then mkdir $OUT/$dir; fi
done

if [ ! -d $QOUT ]; then mkdir $QOUT; fi


############################################
#   FastQC summary of fastq files
#   QC:~/hiSeqInf/pbsTemp/bwaQC.sh qout/bwa/
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bwa/*.$ASD.bam
############################################

if [ -n "$fastQC" ]
then

    echo "********* $TASKFASTQC"
    MAX=16
    
    if [ ! -d $QOUT/$TASKFASTQC ]; then mkdir $QOUT/$TASKFASTQC; fi
    if [ ! -d $OUT/runStats ]; then mkdir $OUT/runStats; fi
    if [ -d $OUT/runStats/$TASKFASTQC ]; then 
	rm -r $OUT/runStats/$TASKFASTQC/
    fi
    mkdir $OUT/runStats/$TASKFASTQC
    
    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')
#    for d in ${DIR[@]}; do
#	FILES=$FILES" "$( ls $OUT/fastq/$d/*$FASTQ )
#    done

#    echo $FILES

#    CPUS=`echo $FILES | wc -w`
#    if [ "$CPUS" -gt "$MAX" ]; then echo "reduce to $MAX CPUs"; CPUS=$MAX; fi
     
    # run fastQC
    if [ -n "$ARMED" ]; then
	qsub -o $QOUT/$TASKFASTQC/$NAME.out -N $TASKFASTQC"_"$NAME -v CONFIG=$CONFIG  $HISEQINF/pbsTemp/fastQC.sh
#	$BINQSUB -j oe -o $QOUT/$TASKFASTQC/$NAME.out -w $(pwd) -N $TASKFASTQC"_"$NAME -l nodes=2:ppn=8 -l vmem=20G -l walltime=2:00:00 \
#	    -command "$FASTQC --nogroup -t $CPUS --outdir $OUT/runStats/$TASKFASTQC `echo $FILES`"
    fi

fi

############################################
#   FASTXTOOLKIT based trimming to a certain length
############################################

if [ -n "$RUNTRIMMING" ]
then

    TASKTRIMMING="fastxtrim"
    echo "********* $TASKTRIMMING"

    for dir in ${DIR[@]}; do

      if [ ! -d $QOUT/$TASKTRIMMING ]; then mkdir $QOUT/$TASKTRIMMING; fi
      if [ ! -d $OUT/fastq/$dir$TRIMLENGTH ]; then mkdir $OUT/fastq/$dir$TRIMLENGTH; fi

      for f in $( ls $SOURCE/fastq/$dir/*$FASTQ ); do
	  n=`basename $f`
	  NAME=$dir"_"$n
	  echo $NAME

	  if [ -e $QOUT/$TASKTRIMMING/$NAME.out ]; then rm $QOUT/$TASKTRIMMING/$NAME.out; fi

          # run trimmer -z for gz output
	  if [ -n "$ARMED" ]; then
	      if [ $FASTQ=="fastq.gz" ]; then
		  echo "gunzip -c $f >${f/.gz/}" >$f.tmp
		  echo "rm $f.tmp" >>$f.tmp
		  chmod -u=rwx $f.tmp
		  qsub $PRIORITY -j y -o $QOUT/$TASKTRIMMING/$NAME.out -cwd -b y -N "unzip_"$NAME \
		      $f.tmp
		  HOLD="-hold_jid unzip_"$NAME
		  f=${f/.gz/}
		  n=${n/.gz/}
	      fi
	      qsub $PRIORITY -j y -o $QOUT/$TASKTRIMMING/$NAME.out -cwd -b y -N $TASKTRIMMING"_"$NAME $HOLD\
		  $FASTXTK/fastx_trimmer -f 1 -l $TRIMLENGTH -i $f -o $OUT/fastq/$dir$TRIMLENGTH/$n -Q 33

	      if [ -n "$ZIPAGAIN" ]; then
		  qsub $PRIORITY -j y -o $QOUT/$TASKTRIMMING/$NAME.out -cwd -b y -N "zip_"$NAME -hold_jid $TASKTRIMMING"_"$NAME $HOLD\
	               gzip $OUT/fastq/$dir$TRIMLENGTH/$n
	      fi

	  fi
      done
  done

fi

############################################
#   FASTXTOOLKIT based trimming to a certain length
############################################

if [ -n "$RUNCUSTOMPLEX" ]
then
    TASKCUSTOMPLEX="cdmult"
    echo "********* $TASKCUSTOMPLEX"
    for dir in ${DIR[@]}; do

      if [ ! -d $QOUT/$TASKCUSTOMPLEX ]; then mkdir $QOUT/$TASKCUSTOMPLEX; fi
      if [ ! -d $OUT/fastq/$dir$TASKCUSTOMPLEX ]; then mkdir $OUT/fastq/$dir$TASKCUSTOMPLEX; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ ); do
	  n=`basename $f`
	  p=${n/"_"$READONE.$FASTQ/}
	  NAME=$dir"_"$n
	  echo $NAME

          # run trimmer -z for gz output
	  if [ -n "$ARMED" ]; then

	      if [ -e $QOUT/$TASKCUSTOMPLEX/$NAME.out ]; then rm $QOUT/$TASKCUSTOMPLEX/$NAME.out; fi

	      qsub $PRIORITY -j y -o $QOUT/$TASKCUSTOMPLEX/$NAME.out -cwd -b y -N $TASKCUSTOMPLEX"_"$NAME \
		  $HISEQINF/pbsTemp/customplex.sh -k $CONFIG -f $f -b $CUSTOMBARCODE -p $p -o $OUT/fastq/$dir$TASKCUSTOMPLEX

	  fi
      done
   done	      
fi


############################################
# downsampling
############################################
if [ -n "$RUNDOWNSAMPLING" ]; then

    echo "********** Downsampling"
    let ABB=$READNUMBER/100000
    TASKDOWNSAMPLE="sample"$ABB"K"
    for dir in ${DIR[@]}; do

      if [ ! -d $QOUT/$TASKDOWNSAMPLE ]; then mkdir $QOUT/$TASKDOWNSAMPLE; fi
      if [ ! -d $OUT/$dir/$TASKDOWNSAMPLE ]; then mkdir $OUT/$dir/$TASKDOWNSAMPLE; fi


      for f in $( ls $OUT/$dir/$TASKRCA/*$ASR.bam ); do
          n=`basename $f`
          NAME=$dir"_"$n
          echo $dir"/"$TASKRCA"/"$NAME" -> "$dir/$TASKDOWNSAMPLE

          if [ -e $QOUT/$TASKDOWNSAMPLE/$NAME.out ]; then rm $QOUT/$TASKDOWNSAMPLE/$NAME.out; fi

          if [ -n "$ARMED" ]; then
	  		 qsub $PRIORITY -j y -o $QOUT/$TASKDOWNSAMPLE/$NAME.out -cwd -b y \
	  		 -N $TASKDOWNSAMPLE"_"$NAME -l vf=4G \
			$HISEQINF/pbsTemp/downsample.sh -k $HISEQINF -i $f -o $OUT/$dir/$TASKDOWNSAMPLE \
			-s $READNUMBER -r $FASTA

		  fi
      done
  done

fi


############################################
#   Mapping using BWA
#   QC:~/hiSeqInf/pbsTemp/bwaQC.sh qout/bwa/
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bwa/*.$ASD.bam
############################################

if [ -n "$RUNMAPPINGBWA" ]
then

    echo "********* $TASKBWA"

    if [ ! -d $QOUT/$TASKBWA ]; then mkdir $QOUT/$TASKBWA; fi

    for dir in ${DIR[@]}
      do
      
      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKBWA ]; then mkdir $OUT/$dir/$TASKBWA; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ); do
      #for f in $(ls $SOURCE/fastq/$dir/ERA000206_hg_$READONE.$FASTQ ); do
	n=`basename $f`
	name=${n/'_'$READONE.$FASTQ/}
	echo ">>>>>"$dir$n
		
	#Sumit (ca. 30 min on AV 1297205 reads) note ${dir/_seqs/} is just a prefix for the
	# readgroup name -l h_rt=03:00:00  -pe mpich 10
	if [ -n "$ARMED" ]; then

            # remove old pbs output
            if [ -e $QOUT/$TASKBWA/$dir'_'$name.out ]; then rm $QOUT/$TASKBWA/$dir'_'$name.*; fi

	    $BINQSUB -j oe -o $QOUT/$TASKBWA/$dir'_'$name'.out' -w $(pwd) -l $NODES_BWA \
		-l vmem=$MEMORY_BWA"G" -N $TASKBWA'_'$dir'_'$name -l walltime=$WALLTIME_BWA \
		-command "$HISEQINF/pbsTemp/bwa.sh $BWAADDPARM -k $CONFIG -t $CPU_BWA -m $(expr $MEMORY_BWA - 1 ) -f $f -r $FASTA \
                -o $OUT/$dir/$TASKBWA --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi $dir \
                --fastqName $FASTQ -R $SEQREG"
	fi
      done
    done

fi

if [ -n "$RUNMAPPINGBWA2" ]; then
    $HISEQINF/pbsTemp/pbsTemp.sh $ARMED -k $CONFIG -t $TASKBWA -o fastq -e "_"$READONE.$FASTQ -n $NODES_BWA -m $MEMORY_BWA"G" -w $WALLTIME_BWA \
        --command "$HISEQINF/pbsTemp/bwa.sh $BWAADDPARM -k $CONFIG -t $CPU_BWA -m $(expr $MEMORY_BWA - 1 ) -f <FILE> -r $FASTA \
                -o $OUT/<DIR>/$TASKBWA --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR> \
                --fastqName $FASTQ -R $SEQREG"
fi


############################################
#   Mapping using BOWTIE
############################################

if [ -n "$RUNMAPPINGBOWTIE" ]
then
    echo "********* $TASKBOWTIE"
    CPUS=$CPU_BOWTIE

    if [ ! -d $QOUT/$TASKBOWTIE ]; then mkdir $QOUT/$TASKBOWTIE; fi

    for dir in ${DIR[@]}; do

      if [ ! -d $OUT/$dir/$TASKBOWTIE ]; then mkdir $OUT/$dir/$TASKBOWTIE; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ); do
        n=`basename $f`
        name=${n/'_'$READONE.$FASTQ/}
        echo ">>>>>"$dir$n
        # remove old pbs output
	if [ -e $QOUT/$TASKBOWTIE/$dir'_'$name.out ]; then rm $QOUT/$TASKBOWTIE/$dir'_'$name.*; fi
        if [ -n "$ARMED" ]; then
           $BINQSUB -j oe -o $QOUT/$TASKBOWTIE/$dir'_'$name'.out' -w $(pwd) -l $NODES_BOWTIE \
                -l vmem=$MEMORY_BOWTIE"G" -N $TASKBOWTIE'_'$dir'_'$name -l walltime=$WALLTIME_BOWTIE \
                -command "$HISEQINF/pbsTemp/bowtie2.sh $BWAADDPARM -k $CONFIG -t $CPU_BOWTIE -m $(expr $MEMORY_BOWTIE - 1 ) -f $f -r $FASTA \
                -o $OUT/$dir/$TASKBOWTIE --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi $dir \
                --fastqName $FASTQ"
        fi
      done
    done

fi



if [ -n "$RUNMAPPINGBOWTIE2" ]; then
    $HISEQINF/pbsTemp/pbsTemp.sh $ARMED -k $CONFIG -t $TASKBOWTIE -o fastq -e "_"$READONE.$FASTQ -n $NODES_BOWTIE -m $MEMORY_BOWTIE"G" -w $WALLTIME_BOWTIE \
	--command "$HISEQINF/pbsTemp/bowtie2.sh $BOWTIEADDPARM -k $CONFIG -t $CPU_BOWTIE -m $(expr $MEMORY_BOWTIE - 1 ) -f <FILE> -r $FASTA -o $OUT/<DIR>/$TASKBOWTIE \
        --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR> --fastqName <NAME>"
fi


############################################
# combine flowcell lenes files for recalibration
#
############################################

if [ -n "$mergeBWAbams" ]; then
    # make tmp files for first step combining
    ./combined/mergeguide/README
    
    # combine them in lane bam files
    for e in $( ls $OUT/combined/mergeguide/lanes/ ); do
	merge.sh $HISEQINF $OUT/combined/mergeguide/lanes/$e $OUT/combined/bwa/ ${e/tmp/$ASD.bam} bam qout/merged/
    done
fi


########
# run tophat,
# already upgraded to take the config rather than the toolkit CAREFULL to move $OUT ->$OUTDIR 
# already cleanup in armed only
########
if [ -n "$RUNTOPHATCUFF" ]; then
    
    echo "********* $TASKTOPHAT"
    
    #CPUS_TOPHAT=24
 
    # ensure directory is there
    if [ ! -d $QOUT/$TASKTOPHAT ]; then mkdir $QOUT/$TASKTOPHAT; fi
    for dir in ${DIR[@]}
      do
       
      if [ ! -d $OUT/$dir/$TASKTOPHAT ]; then mkdir $OUT/$dir/$TASKTOPHAT; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	name=${n/'_'$READONE.$FASTQ/}
	echo $name
	echo ">>>>>"$dir"/"$TASKTOPHAT"/"$n" ("$TASKCUFF")"
	
        #submit
	if [ -n "$ARMED" ]; then

            #cleanup old qouts
	    if [ -e $QOUT/$TASKTOPHAT/$dir'_'$name.out ]; then rm $QOUT/$TASKTOPHAT/$dir'_'$name.out; fi

	    $BINQSUB -j oe -o $QOUT/$TASKTOPHAT/$dir'_'$name.out -w $(pwd) -l walltime=$WALLTIME_TOPHAT \
		-N $TASKTOPHAT"_"$dir"_"$name -l $NODES_TOPHAT -l vmem=$MEMORY_TOPHAT"G" \
		-command "$HISEQINF/pbsTemp/tophatcuff.sh $TOPHATADDPARM -k $CONFIG -r $FASTA -f $f \
		-t $CPU_TOPHAT -o $OUT/$dir/$TASKTOPHAT/$name/ -a $REFSEQGTF"
	    #exit
	fi
      done


    done
fi


if [ -n "$RUNTOPHATCUFF2" ]; then
    $HISEQINF/pbsTemp/pbsTemp.sh $ARMED -k $CONFIG -t $TASKTOPHAT -o fastq -e "_"$READONE.$FASTQ -n $NODES_TOPHAT -m $MEMORY_TOPHAT"G" -w $WALLTIME_TOPHAT \
        --command "$HISEQINF/pbsTemp/tophatcuff.sh $TOPHATADDPARM -k $CONFIG -f <FILE> \
         -t $CPU_TOPHAT -o $OUT/<DIR>/$TASKTOPHAT/<NAME> "

fi




############################################
#   recalibrate quality scores OLD
#   http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
#   http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
#   http://picard.sourceforge.net/command-line-overview.shtml#FixMateInformation
#   full pipe: http://www.broadinstitute.org/gsa/wiki/index.php/Whole_genome,_deep_coverage
# IN:$SOURCE/$dir/fastq/*$READONE.fastq
# OUT: $OUT/$dir/reCal/*.$$ASR.bam
############################################
if [ -n "$recalibrateQualScore" ]; then

    echo "********* $TASKRCA"

    # generates the .dict file
    #java -Xmx4g -jar /home/Software/picard-tools-1.22/CreateSequenceDictionary.jar R= $FASTA O= ${$FASTA/.fasta/.dict}
    #sort dbSNP rod file according to fasta
    #/home/Software/Sting/perl/sortByRef.pl --k 2 $GATKSUP/tmp.dbsnp.txt $FASTA.fai > $DBROD

    if [ ! -d $QOUT/$TASKRCA ]; then mkdir $QOUT/$TASKRCA; fi
    #ensure dirs are there...
    if [ ! -d $OUT/combined/$TASKRCA ]; then mkdir $OUT/combined/$TASKRCA; fi

    for f in $( ls $OUT/combined/bwa/*$ASD.bam )
#    for f in $( ls combined/bwa/Nina_FC3056JAAXX_l6.$ASD.bam )
	do
	n=`basename $f`
	name=${n/.bam/}
	echo ">>>>>"$n


	# wait on pipeline steps
	HOLD="-hold_jid *mergeBams*"
	# remove old pbs output
	if [ -e $QOUT/$TASKRCA/$name.out ]; then rm $QOUT/$TASKRCA/$name.*; fi

	#Sumit (ca. 80 min on AV 1297205 reads) -l h_rt=20:00:00
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKRCA/$name'.out' -cwd -b y -l vf=20G \
		-l mem_free=20G -l h_vmem=20G -N $TASKRCA'_'$name $HOLD\
		$HISEQINF/pbsTemp/reCalAln.sh $HISEQINF $f $FASTA $DBROD \
		$OUT/combined/$TASKRCA $SEQREG
	fi
	
    done

fi




############################################ 
#   recalibrate quality scores
#   http://www.broadinstitute.org/gsa/wiki/index.php/Base_quality_score_recalibration
#   http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels
#   http://picard.sourceforge.net/command-line-overview.shtml#FixMateInformation
#   full pipe: http://www.broadinstitute.org/gsa/wiki/index.php/Whole_genome,_deep_coverage
# IN:$SOURCE/$dir/fastq/*$READONE.fastq
# OUT: $OUT/$dir/reCal/*.$$ASR.bam
############################################
if [ -n "$RUNREALRECAL" ]; then

    echo "********* $TASKRCA"
    CPU=16
    exit

    # generates the .dict file
    #java -Xmx4g -jar /home/Software/picard-tools-1.22/CreateSequenceDictionary.jar R= $FASTA O= ${$FASTA/.fasta/.dict}
    #sort dbSNP rod file according to fasta
    #/home/Software/Sting/perl/sortByRef.pl --k 2 $GATKSUP/tmp.dbsnp.txt $FASTA.fai > $DBROD

    if [ ! -d $QOUT/$TASKRCA ]; then mkdir $QOUT/$TASKRCA; fi

    for dir in ${DIR[@]}
      do

     #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKRCA ]; then mkdir $OUT/$dir/$TASKRCA; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	n2=${n/'_'$READONE.$FASTQ/.$ASD.bam}
	name=${n/'_'$READONE.$FASTQ/}
	echo ">>>>>"$dir$n2

	# wait on pipeline steps
	if [ -n "$RUNMAPPINGBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi
	# remove old pbs output
	if [ -e $QOUT/$TASKRCA/$dir'_'$name.out ]; then rm $QOUT/$TASKRCA/$dir'_'$name.*; fi
			
	#Sumit (ca. 80 min on AV 1297205 reads) -l h_rt=20:00:00
	if [ -n "$ARMED" ]; then
	    $BINQSUB -j oe -o $QOUT/$TASKRCA/$dir'_'$name'.out' -w $(pwd) -l $NODES_RECAL \
		-l vmem=$MEMORY_RECAL'G' -N $TASKRCA'_'$dir'_'$name -l walltime=$WALLTIME_RECAL \
		-command "$HISEQINF/pbsTemp/reCalAln.sh -k $HISEQINF -f $OUT/$dir/$TASKBWA/$n2 -r $FASTA -d $DBROD \
		-o $OUT/$dir/$TASKRCA -t $CPU_RECAL $RECALADDPARAM"
	fi
	
      done
    done

fi


if [ -n "$RUNREALRECAL2" ]; then

    $HISEQINF/pbsTemp/pbsTemp.sh $ARMED -r -k $CONFIG -t $TASKRCA -o $TASKBWA/ -e .$ASD.bam \
        -n $NODES_RECAL -m $MEMORY_RECAL"G" -w $WALLTIME_RECAL \
        --command "$HISEQINF/pbsTemp/reCalAln.sh $RECALADDPARAM -k $CONFIG -f <FILE> -r $FASTA -d $DBROD -o $OUT/<DIR>/$TASKRCA -t $CPU_RECAL"


fi




############################################
# combine all into one bamfile
#
############################################

if [ -n "$mergeReCalbams" ]; then
    # make tmp files for first step combining
    ls combined/reCalAln/*.bam >combined/mergeguide/combineAll.txt
    
    # combine them in lane bam files
    $HISEQINF/pbsTemp/merge.sh $HISEQINF $OUT/combined/mergeguide/combineAll.txt $OUT/combined/ DISC1_all.bam bam qout/merged/
fi



############################################
# DepthOfCoverage
# expects to be run fom <dir>/<TASKRCA>/<name>.<ASR>.bam
# e.g. Run/reCalAln/name.ashrr.bam
# change that by setting TASKRCA=TASKBWA and ASR=ASD
############################################

if [ -n "$DEPTHOFCOVERAGE" ]
then
    CPUS=24

    echo "********* $TASKDOC"

    if [ ! -d $QOUT/$TASKDOC ]; then mkdir $QOUT/$TASKDOC; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/'_'$READONE.$FASTQ/.$ASR.bam}
	name=${n/'_'$READONE.$FASTQ/}
	echo ">>>>>"$dir/$TASKRCA/$n2

	if [ ! -d $OUT/$dir/$TASKDOC ]; then mkdir $OUT/$dir/$TASKDOC; fi

	# remove old pbs output
	if [ -e $QOUT/$TASKDOC/$dir'_'$name'.out' ]; then rm $QOUT/$TASKDOC/$dir'_'$name'.out'; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$$RUNREALRECAL" ]; then HOLD="-hold_jid "$TASKRCA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKDOC/$dir'_'$name'.out' -cwd -b y -pe mpich $CPUS \
		-l mem_free=11G -l h_vmem=11G -l vf=500K -N $TASKDOC'_'$dir'_'$name $HOLD\
		$HISEQINF/pbsTemp/gatkDOC.sh -k $HISEQINF -f $OUT/$dir/$TASKRCA/$n2 -r $FASTA \
		-o $OUT/$dir/$TASKDOC -t $CPUS $DOCADDPARAM
	fi

      done
    done

fi

if [ -n "$DEPTHOFCOVERAGE2" ]; then

    $HISEQINF/pbsTemp/pbsTemp.sh $ARMED -r -k $CONFIG -t $TASKDOC -o $TASKRCA/ -e .$ASR.bam \
	-n $NODES_GATKDOC -m $MEMORY_GATKDOC"G" -w $WALLTIME_GATKDOC \
	--command "$HISEQINF/pbsTemp/gatkDOC.sh $DOCADDPARAM -k $HISEQINF -f <FILE> -r $FASTA -o $OUT/<DIR>/$TASKDOC -t $CPU_GATKDOC"

fi



############################################
# downsample
############################################

if [ -n "$DOWNSAMPLE" ]
then

    echo "********* $TASKDOWN"

    if [ ! -d $QOUT/$TASKDOWN ]; then mkdir $QOUT/$TASKDOWN; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/'_'$READONE.$FASTQ/.$ASD.bam}
	name=${n/'_'$READONE.$FASTQ/}
	echo ">>>>>"$dir$n2

	if [ ! -d $OUT/$dir/$TASKDOWN ]; then mkdir $OUT/$dir/$TASKDOWN; fi

	# remove old pbs output
	if [ -e $QOUT/$TASKDOWN/$dir"_"$n2"0.out" ]; then rm $QOUT/$TASKDOWN/$dir"_"$n2.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$mappingBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    #downsample
	    python $HISEQINF/bin/downsample.py -i $OUT/$dir/bwa/$n2 -o $OUT/$dir/$TASKDOWN/ -t downsample/$dir --region $SEQREG -w 500 -s 500 -q $QOUT/$TASKDOWN/$dir

	fi

      done
    done

# echo $( ls Hannibal_FC30MEJAAXX/downsample/*500d.bam ) > mergeBamfiles500d.tmp
# qsub -b y -cwd -j y -o qout/$COMBDIR/merged.asd.500d -pe mpich 2 /clusterdata/hiseq_apps/hiSeqInf/pbsTemp/merge.sh $HISEQINF mergeBamfiles500d.tmp merged merged.$ASD.500d.bam bam

fi


############################################
# combine dindel vcfs
############################################

if [ -n "$DINDELcombineVCF" ]
then

    echo "********* $TASKDINS combine"

    if [ ! -d $QOUT/$TASKDINS ]; then mkdir $QOUT/$TASKDINS; fi
    if [ -e $OUT/dindelVCFmerge.tmp ]; then rm $OUT/dindelVCFmerge.tmp; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	n2=${n/'_'$READONE.$FASTQ/.ashrr.bam.dindel.VCF}
	echo "$OUT/$dir/dindelS/$n2" >> dindelVCFmerge.tmp
      done
    done

    # remove old pbs output
    if [ -e $QOUT/$TASKDINS/dindelVCFmerge.out ]; then rm $QOUT/$TASKDINS/dindelVCFmerge.out; fi

    #Submit
    if [ -n "$ARMED" ]; then
	qsub $PRIORITY -j y -o $QOUT/$TASKDINS/dindelVCFmerge.out -cwd -b y -l h_vmem=12G -N dindelVCFmerge \
	    $HISEQINF/pbsTemp/merge.sh $HISEQINF dindelVCFmerge.tmp $OUT/genotype dindelSeparate.vcf vcf
	    
    fi

fi


############################################
# combine dindel vcfs
############################################

if [ -n "$DINDELcombine" ]
then

    echo "********* $TASKDINS combine"

    if [ ! -d $QOUT/$TASKDINS ]; then mkdir $QOUT/$TASKDINS; fi
    if [ -e $OUT/dindelVCFmerge.tmp ]; then rm $OUT/dindelVCFmerge.tmp; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	n2=${n/'_'$READONE.fastq/.ashrr.bam.dindel.VCF}
	#qsub -j y -o $QOUT/$TASKDINS/reg$dir"_"$n2.out -cwd -b y -N reg$dir"_"$n2.out \
	    python $DINDELHOME/dindel-1.01-python/convertVCFToDindel.py -i $OUT/$dir/dindelS/$n2 \
	        -o $OUT/$dir/dindelS/$n2.reg -r $FASTA
	echo "$OUT/$dir/dindelS/$n2.reg" >> dindelVCFmerge.tmp
      done
    done

    # remove old pbs output
    if [ -e $QOUT/$TASKDINS/dindelVCFmerge.out ]; then rm $QOUT/$TASKDINS/dindelVCFmerge.out; fi

    #Submit
    if [ -n "$ARMED" ]; then
	#echo "cat $( cat dindelVCFmerge.tmp) | sort -k 1,1 -k 2,2n -u > $OUT/genotype/dindelSeparate.cand" >qsubmerge.tmp
	#echo "rm qsubmerge.tmp" >>qsubmerge.tmp
	#chmod -u=rwx qsubmerge.tmp
	#qsub -j y -o $QOUT/$TASKDINS/dindelVCFmerge.out -cwd -b y -N dindelVCFmerge -hold_jid reg* \
	#    ./qsubmerge.tmp
	    cat $( cat dindelVCFmerge.tmp) | sort -k 1,1 -k 2,2n -u > $COMBDIR/dindelSeparate.variants.txt
    fi

fi


############################################
#   call indels with DINDEL on all individuals combined
#   QC: ~/hiSeqInf/pbsTemp/dindelQC.sh qout/dindel/
############################################

if [ -n "$DINDELcallIndelsCombined" ]
then

    echo "********* $TASKDINC"

    if [ ! -d $QOUT/$TASKDINC ]; then mkdir $QOUT/$TASKDINC; fi
    if [ ! -d $COMBDIR/$TASKDINC ]; then mkdir $COMBDIR/$TASKDINC; fi
 
    for dir in ${DIR[@]}; do
	CANDIDATES=$CANDIDATES" "`echo $(ls $dir/dindel/*dindel.VCF)`
    done
    CANDIDATES=${CANDIDATES// /,}

    #Submit
    if [ -n "$ARMED" ]; then
	    $HISEQINF/pbsTemp/dindelV2.sh $HISEQINF $CANDIDATES $FASTA $COMBDIR/$TASKDINC/

    fi
fi



############################################
#   call indels with GATK
############################################

#TODO run them together with gatkIndelV2.sh

if [ -n "$GATKcallIndelsSeperate" ]
then

    echo "********* $TASKIND"
    if [ ! -d $QOUT/$TASKIND ]; then mkdir $QOUT/$TASKIND; fi
    if [ -e $QOUT/$TASKIND/ids.txt ]; then rm $QOUT/$TASKIND/ids.txt; fi
    if [ -e $QOUT/$TASKIND/sum.out ]; then rm $QOUT/$TASKIND/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKIND ]; then mkdir $OUT/$dir/$TASKIND; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq )
	do
	n=`basename $f`
	n2=${n/'_'$READONE.fastq/.$ASR.bam}
	name=${n/'_'$READONE.fastq/}
	echo ">>>>>"$dir$n2

	# remove old pbs output
	if [ -e $QOUT/$TASKIND/$dir"_"$name.out ]; then rm $QOUT/$TASKIND/$dir"_"$name.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$RUNMAPPINGBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi
	if [ -n "$recalibrateQualScore" ]; then HOLD="-hold_jid "$TASKRCA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKIND/$dir'_'$name'.out' -cwd -b y \
		-l h_vmem=12G -N $TASKIND'_'$dir'_'$name $HOLD\
		$HISEQINF/pbsTemp/gatkIndel.sh $HISEQINF $OUT/$dir/$TASKRCA/$n2 $FASTA $DBROD \
		$REFSEQROD $OUT/$dir/$TASKIND
	fi


      done
    done


    #combine into one file
#    qsub -j y -o $QOUT/gatkInd/final.out -cwd -b y -hold_jid `cat $QOUT/gatkInd/ids.txt`\
#	echo "merge"
#		merge-vcf `ls $QOUT/dindel/*.dindel.VCF` > genotyping/indels.dindel.vcf
#		merge-vcf A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

fi


############################################
#   call indels with GATK
############################################

#TODO run them together with gatkIndelV2.sh

if [ -n "$GATKcallIndelsCombined" ]
then

    echo "********* $TASKIND"

    if [ ! -d $QOUT/$TASKIND ]; then mkdir $QOUT/$TASKIND; fi
    if [ -e $QOUT/$TASKIND/ids.txt ]; then rm $QOUT/$TASKIND/ids.txt; fi
    if [ ! -d $OUT/genotype ]; then mkdir $OUT/genotype; fi

    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')

    for dir in ${DIR[@]};do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq );do
	n=`basename $f`
	echo $OUT/$dir/$TASKRCA/${n/'_'$READONE.fastq/.$ASR.bam} >> $TASKIND"bamfiles.tmp"
      done
    done

    # remove old pbs output
    if [ -e $QOUT/$TASKIND/$NAME.out ]; then rm $QOUT/$TASKIND/$NAME.out ; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    if [ -n "$recalibrateQualScore" ]; then 
	HOLD="-hold_jid "
	for dir in ${DIR[@]};do
	    HOLD=$HOLD","$TASKRCA"_"$dir"*"
        done
	HOLD=${HOLD/,/}
    fi

    #Submit
    if [ -n "$ARMED" ]; then
	qsub -j y -o $QOUT/$TASKIND/$NAME.out -cwd -b y \
	    -l mem_free=20G -l h_vmem=20G -N $TASKIND"_"$NAME.out $HOLD\
	    $HISEQINF/pbsTemp/gatkIndelV2.sh $HISEQINF $TASKIND"bamfiles.tmp" $FASTA $DBROD \
	    $REFSEQROD $OUT/genotype $SEQREG
    fi

    
fi


############################################
#   call indels with GATK -- call one vcf file over all folders
############################################

if [ -n "$RUNVARCALLS" ]
then

    echo "********* $TASKVAR"
    CPUS=24

    if [ ! -d $QOUT/$TASKVAR ]; then mkdir $QOUT/$TASKVAR; fi
    if [ ! -d $OUT/$TASKVAR ]; then mkdir $OUT/$TASKVAR; fi

    RUNVARS=""

    # one per folder
    if [ -n "$VARPERFOLDER" ]; then
	for dir in ${DIR[@]};do
	    RUNVARS=$RUNVARS" "$dir
	    if [ ! -d $OUT/$TASKVAR/$dir ]; then mkdir $OUT/$TASKVAR/$dir; fi
            if [ -e  $OUT/$TASKVAR/$dir/$TASKVAR"bamfiles.tmp" ]; then rm  $OUT/$TASKVAR/$dir/$TASKVAR"bamfiles.tmp"; fi
	    for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );do
		n=`basename $f`
		echo $OUT/$dir/$TASKRCA/${n/'_'$READONE.$FASTQ/.$ASR.bam} >> $OUT/$TASKVAR/$dir/$TASKVAR"bamfiles.tmp"
	    done
	done
    # one over all folders
    else
	NAME=$(echo ${DIR[@]}|sed 's/ /_/g')
	RUNVARS=$NAME
	if [ ! -d $OUT/$TASKVAR/$NAME ]; then mkdir $OUT/$TASKVAR/$NAME; fi
        if [ -e  $OUT/$TASKVAR/$NAME/$TASKVAR"bamfiles.tmp" ]; then rm  $OUT/$TASKVAR/$NAME/$TASKVAR"bamfiles.tmp"; fi
	for dir in ${DIR[@]};do
	    for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );do
		n=`basename $f`
		echo $OUT/$dir/$TASKRCA/${n/'_'$READONE.$FASTQ/.$ASR.bam} >> $OUT/$TASKVAR/$NAME/$TASKVAR"bamfiles.tmp"
	    done
	done
	
    fi

    for dir in $RUNVARS; do #$(ls -d $OUT/$TASKVAR/* ); do
	NAME=`basename $dir`
        # remove old pbs output
	if [ -e $QOUT/$TASKVAR/$NAME.out ]; then rm $QOUT/$TASKVAR/$NAME.out ; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    #if [ -n "$recalibrateQualScore" ]; then 
#	HOLD="-hold_jid "
#	for dir in ${DIR[@]};do
#	    HOLD=$HOLD","$TASKRCA"_"$dir"*"
 #       done
#	HOLD=${HOLD/,/}
 #   fi

        #Submit #-pe mpich $CPUS 
	if [ -n "$ARMED" ]; then
	    $BINQSUB -j oe -o $QOUT/$TASKVAR/$NAME.out -w $(pwd) -l $NODES_VAR $HOLD \
		-l vmem=$MEMORY_VAR"G" -N $TASKVAR'_'$NAME -l walltime=$WALLTIME_VAR \
		-command "$HISEQINF/pbsTemp/gatkSNPs.sh -k $CONFIG -i $OUT/$TASKVAR/$NAME/$TASKVAR'bamfiles.tmp' -t $CPU_VAR \
		-r $FASTA -d $DBROD -o $OUT/$TASKVAR/$NAME -n $NAME -H $HAPMAPVCF -K $ONEKGVCF $VARADDPARAM"
	fi
    done

    
fi


########
# run differental expression detection
########
if [ -n "$RUNANNOTATION" ]; then

    # ensure directory is there
    if [ ! -d $QOUT/$TASKANNOVAR ]; then mkdir $QOUT/$TASKANNOVAR; fi
    if [ ! -d $OUT/$TASKANNOVAR ]; then mkdir $OUT/$TASKANNOVAR; fi

    
#    for dir in $( ls -d $TASKVAR/* ); do
     for d in ${DIR[@]}; do	

	dir=$OUT/$TASKVAR/$d

	n=`basename $dir`
	echo $n

	namesnp=$OUT/$TASKVAR/$n/$n.filter.snps.vcf
	namesnp2=$OUT/$TASKVAR/$n/$n.recalfilt.snps.vcf
	nameindel=$OUT/$TASKVAR/$n/$n.filter.indel.vcf

        #cleanup old qouts
        if [ -e $QOUT/$TASKANNOVAR/$n.out ]; then rm $QOUT/$TASKANNOVAR/$n.out; fi
	#ensure subfolder is there
	if [ ! -d $OUT/$TASKANNOVAR/$n ]; then mkdir $OUT/$TASKANNOVAR/$n; fi


        #check if this is part of the pipe and jobsubmission needs to wait
        if [ -n "$TASKVAR" ]; then HOLD="-hold_jid "$TASKVAR"_"$n; fi
	#echo "-hold_jid "$TASKVAR"_"$n
	HOLD="-hold_jid "$TASKVAR"_"$n

        #submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASKANNOVAR/$n.out \
		-N $TASKANNOVAR'_'$n $HOLD\
		$HISEQINF/pbsTemp/annovar.sh -k $HISEQINF -i1 $namesnp -i2 $namesnp2 -i3 $nameindel -r $FASTA \
		-o $OUT/$TASKANNOVAR/$n 

	fi
   done

fi


############################################
#   call snps with GATK
############################################

if [ -n "$GATKcallSNPS" ]
then

    echo "********* $TASKSNP"

    if [ ! -d $QOUT/$TASKSNP ]; then mkdir $QOUT/$TASKSNP; fi
    if [ -e $QOUT/$TASKSNP/ids.txt ]; then rm $QOUT/$TASKSNP/ids.txt; fi
    if [ -e $QOUT/$TASKSNP/sum.out ]; then rm $QOUT/$TASKSNP/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKSNP ]; then mkdir $OUT/$dir/$TASKSNP; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	n2=${n/'_'$READONE.$FASTQ/.$ASR.bam}
	name=${n/'_'$READONE.$FASTQ/}
	echo ">>>>>"$dir$n2

	# remove old pbs output
	if [ -e $QOUT/$TASKSNP/$dir"_"$name.out ]; then rm $QOUT/$TASKSNP/$dir"_"$name.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$GATKcallIndelsSeperate" ] || [ -n "$GATKcallIndelsCombined" ]
	    then HOLD="-hold_jid "$TASKIND"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKSNP/$dir'_'$name'.out' -cwd -b y \
		-l h_vmem=12G -N $TASKSNP'_'$dir'_'$name $HOLD\
		$HISEQINF/pbsTemp/gatkSNPs.sh $CONFIG $OUT/$dir/$TASKRCA/$n2 $FASTA $DBVCF \
		$REFSEQROD $OUT/$dir/$TASKSNP $OUT/$dir/$TASKIND | cut -d "" -f 2 >>$QOUT/$TASKSNP/ids.txt
	fi


      done
    done


    #combine into one file
#    qsub -j y -o $QOUT/gatkInd/final.out -cwd -b y -hold_jid `cat $QOUT/gatkInd/ids.txt`\
#	echo "merge"
#		merge-vcf `ls $QOUT/dindel/*.dindel.VCF` > genotyping/indels.dindel.vcf
#		merge-vcf A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

fi




########
# run differental expression detection
########
if [ -n "$RUNCUFFDIFF" ]; then

    CPUS=24

    # ensure directory is there
    if [ ! -d $QOUT/$TASKCUFFDIFF ]; then mkdir $QOUT/$TASKCUFFDIFF; fi
    if [ ! -d $OUT/$TASKCUFFDIFF ]; then mkdir $OUT/$TASKCUFFDIFF; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	name=$name${n/'_'$READONE.$FASTQ/}","
      done
    done
    name=$(echo $name | sed 's/\(.*\),/\1/')
    n=${name/,/_u_}

    #cleanup old qouts
    if [ -e $QOUT/$TASKCUFFDIFF/run.out ]; then rm $QOUT/$TASKCUFFDIFF/run.out; fi

    #check if this is part of the pipe and jobsubmission needs to wait
    if [ -n "$RUNTOPHATCUFF" ]; then HOLD="-hold_jid "$TASKTOPHAT"_"$dir"_"$name; fi

    #submit
    if [ -n "$ARMED" ]; then
      qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASKCUFFDIFF/$n.out \
     	    -N $TASKCUFFDIFF'_'$n -pe mpich $CPUS $HOLD\
            $HISEQINF/pbsTemp/cuffdiff.sh -k $HISEQINF -b $name -r $FASTA \
	    -t $CPUS -o $OUT/$TASKCUFFDIFF/$n -a $REFSEQGTF

	fi

fi


############################################
#   Mapping using rrbsmap
#   
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/rrbs/*.rrbsd.bam
############################################

if [ -n "$RUNMAPPINGRRBS" ]
then

    echo "********* $TASKRRBS"
    CPUS=32

    if [ ! -d $QOUT/$TASKRRBS ]; then mkdir $QOUT/$TASKRRBS; fi

    for dir in ${DIR[@]}
      do
      
      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKRRBS ]; then mkdir $OUT/$dir/$TASKRRBS; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	name=${n/'_'$READONE.$FASTQ/}
	echo ">>>>>"$dir$n
		
	# remove old pbs output
	if [ -e $QOUT/$TASKRRBS/$dir'_'$name.out ]; then rm $QOUT/$TASKRRBS/$dir'_'$name.*; fi

	if [ -n "$ARMED" ]; then
	   qsub $PRIORITY -j y -o $QOUT/$TASKRRBS/$dir'_'$name'.out' -cwd -b y -pe mpich $CPUS \
		-l h_vmem=12G -N $TASKRRBS'_'$dir'_'$name \
		$HISEQINF/pbsTemp/rrbsmap.sh $RRBSMAPADDPARM -k $HISEQINF -t $CPUS -f $f -r $FASTA -o $OUT/$dir/$TASKRRBS \
		--rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi $dir --rgpu $FLOWCELL --fastqName $FASTQ 

	fi

      done
    done

fi

############################################ 
# IN */bwa/*.bam
# OUT */bwa/*.ann
############################################ 
if [ -n "$RUNANNOTATINGBAM" ]; then
    $HISEQINF/pbsTemp/pbsTemp.sh --nodir -r $ARMED -k $CONFIG -t $TASKBAMANN --origin $TASKBWA -e .bam -n $NODES_BAMANN \
	-m $MEMORY_BAMANN'G' -w $WALLTIME_BAMANN \
        --command "$HISEQINF/pbsTemp/annotateBam.sh -k $CONFIG -f <FILE>"
fi


############################################
# IN: */bwa/*.bam
# OUT: */bwa_var/*.clean.vcf
############################################

if [ -n "$RUNSAMVAR" ]; then
   $HISEQINF/pbsTemp/pbsTemp.sh -r $ARMED -k $CONFIG -t $TASKBWA-$TASKSAMVAR --origin $TASKBWA -e .$ASD.bam -n $NODES_SAMVAR \
        -m $MEMORY_SAMVAR'G' -w $WALLTIME_SAMVAR \
        --command "$HISEQINF/pbsTemp/samSNPs.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBWA-$TASKSAMVAR" \
	--postcommand "$HISEQINF/pbsTemp/samSNPscollect.sh -k $CONFIG -f <FILE> -o $OUT/variant/$TASKBWA-$TASKSAMVAR-<DIR>"
fi
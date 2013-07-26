#!/bin/bash

# Overall trigger pipeline
# author: Denis C. Bauer
# date: Nov.2010

function usage {
echo -e "usage: $(basename $0) CONFIG [TASK]

Script interpreting the CONFIG file in the project directory and submitting
tasks to the queue.

required:
  CONFIG     config.txt file specifying what needs to be done and
               where the resources are located

options for TASK:
  empty      start dry-run: make dirs, delete old files, print what will be done
  fetchdata  get data from remote server (via smbclient)
  pushresult puts results to remote server (via smbclient)
  armed      submit tasks to the queue
  direct     run task directly (e.g. on node after qrsh)
  postonly   run only the postanalysis step of a TASK
  html       check the cluster logfiles for errors and and make summary HTML page
"
exit
}

if [ ! $# -gt 0 ]; then usage ; fi

CONFIG=$1
ADDITIONALTASK=$2

# get all the specs defined in the config  (note both configs are necessary)
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

############################################
#   Verify
############################################
if [ -n "$ADDITIONALTASK" ]; then
    if [ "$ADDITIONALTASK" = "verify" ]; then
	echo ">>>>>>>>>> $ADDITIONALTASK"
	if [ -n "$RUNMAPPINGBWA" ]; then ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/bwa.sh $QOUT/$TASKBWA; fi
	if [ -n "$recalibrateQualScore" ]; then ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/reCalAln.sh $QOUT/$TASKRCA; fi
	if [ -n "$GATKcallIndels" ]; then ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/gatkIndel.sh $QOUT/$TASKIND; fi
	if [ -n "$GATKcallSNPS" ]; then ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/gatkSNPs.sh $QOUT/$TASKSNP; fi
	if [ -n "$RUNTOPHATCUFF" ]; then ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/tophatcuff.sh $QOUT/$TASKTOPHAT; fi
	if [ -n "$RUNCUFFDIFF" ]; then ${NGSANE_BASE}/mods/QC.sh ${NGSANE_BASE}/mods/cuffdiff.sh $QOUT/$TASKCUFFDIFF; fi
	exit
	elif [ "$ADDITIONALTASK" = "fetchdata" ]; then
	    echo ">>>>>>>>>> $ADDITIONALTASK"
	    ${NGSANE_BASE}/mods/fetchRawDataFromServer.sh -k $CONFIG
	    exit
	elif [ "$ADDITIONALTASK" = "pushresult" ]; then
	    echo ">>>>>>>>>> $ADDITIONALTASK"
	    ${NGSANE_BASE}/mods/pushResultToServer.sh -k $CONFIG
	    exit
    elif [ "$ADDITIONALTASK" = "html" ]; then
	   echo ">>>>>>>>>> $ADDITIONALTASK"
	   ${NGSANE_BASE}/mods/makeSummary.sh ${NGSANE_BASE} $CONFIG
	    exit
    elif [ "$ADDITIONALTASK" = "armed" ]; then
	    echo ">>>>>>>>>> $ADDITIONALTASK"
	    ARMED="--armed"
            echo -n "Double check! Then type 'safetyoff' and hit enter to launch the job: "
            read safetyoff
            if [ "$safetyoff" != "safetyoff" ];then
                echo "Holstering..."
                exit 0
            else
                echo "... take cover!"
            fi
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
	echo -e "[ERROR] don't understand $ADDITIONALTASK"
	exit -1
    fi
fi

# test if source data is defined
echo "${DIR[@]}"
if [[ -z "${DIR[@]}" ]]; then
  echo "[ERROR] no input directories specified (DIR)."
  exit 1
fi

# ensure out directory is there 
for dir in ${DIR[@]}; do
  if [ ! -d $OUT/$dir ]; then mkdir -p $OUT/$dir; fi
done

if [ ! -d $QOUT ]; then mkdir -p $QOUT; fi
if [ ! -d $TMP ]; then mkdir -p $TMP; fi

############################################
#   FastQC summary of fastq files
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $OUT/runstats/fastQC/*
############################################

if [ -n "$RUNFASTQC" ]; then
    if [ -z "$TASKFASTQC" ] || [ -z "$NODES_FASTQC" ] || [ -z "$CPU_FASTQC" ] || [ -z "$MEMORY_FASTQC" ] || [ -z "$WALLTIME_FASTQC" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKFASTQC -i fastq -e $READONE.$FASTQ -n $NODES_FASTQC \
	-c $CPU_FASTQC -m $MEMORY_FASTQC"G" -w $WALLTIME_FASTQC \
	--postcommand "${NGSANE_BASE}/mods/fastQC.sh -k $CONFIG" 
fi

############################################
#   CUTADAPT remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_cutadapt/*read1.fastq
############################################

if [ -n "$RUNCUTADAPT" ]; then
    if [ -z "$TASKCUTADAPT" ] || [ -z "$NODES_CUTADAPT" ] || [ -z "$CPU_CUTADAPT" ] || [ -z "$MEMORY_CUTADAPT" ] || [ -z "$WALLTIME_CUTADAPT" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKCUTADAPT -i fastq -e $READONE.$FASTQ -n $NODES_CUTADAPT \
	-c $CPU_CUTADAPT -m $MEMORY_CUTADAPT"G" -w $WALLTIME_CUTADAPT \
	--command "${NGSANE_BASE}/mods/cutadapt.sh -k $CONFIG -f <FILE>" 
fi

############################################
#   TRIMMOMATIC remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_trimmomatic/*read1.fastq
############################################

if [ -n "$RUNTRIMMOMATIC" ]; then
    $QSUB $ARMED -d -k $CONFIG -t $TASKTRIMMOMATIC -i fastq -e $READONE.$FASTQ -n $NODES_TRIMMOMATIC \
        -c $CPU_TRIMMOMATIC -m $MEMORY_TRIMMOMATIC"G" -w $WALLTIME_TRIMMOMATIC \
        --command "$NGSANE_BASE/mods/trimmomatic.sh -k $CONFIG -f <FILE>"
fi

############################################
#   TRIMGALORE remove contaminants
#
# IN : $SOURCE/fastq/$dir/*read1.fastq
# OUT: $SOURCE/fastq/$dir_trimgalore/*read1.fastq
############################################ 
if [ -n "$RUNTRIMGALORE" ]; then
    if [ -z "$TASKTRIMGALORE" ] || [ -z "$NODES_TRIMGALORE" ] || [ -z "$CPU_TRIMGALORE" ] || [ -z "$MEMORY_TRIMGALORE" ] || [ -z "$WALLTIME_TRIMGALORE" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -d -k $CONFIG -t $TASKTRIMGALORE -i fastq -e $READONE.$FASTQ -n $NODES_TRIMGALORE \
        -c $CPU_TRIMGALORE -m $MEMORY_TRIMGALORE"G" -w $WALLTIME_TRIMGALORE \
        --command "${NGSANE_BASE}/mods/trimgalore.sh -k $CONFIG -f <FILE>"
fi

############################################ 
# IN : */bwa/*.bam
# OUT: */bwa/*.ann
############################################ 
if [ -n "$RUNANNOTATINGBAM" ]; then
    if [ -z "$TASKBWA" ] || [ -z "$NODES_BAMANN" ] || [ -z "$CPU_BAMANN" ] || [ -z "$MEMORY_BAMANN" ] || [ -z "$WALLTIME_BAMANN" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB --nodir -r $ARMED -k $CONFIG -t $TASKBAMANN -i $TASKBWA -e .bam \
    -n $NODES_BAMANN -c $CPU_BAMANN -m $MEMORY_BAMANN'G' -w $WALLTIME_BAMANN \
        --command "${NGSANE_BASE}/mods/annotateBam.sh -k $CONFIG -f <FILE>"
fi


############################################
# IN: */bwa/*.bam
# OUT: */bwa_var/*.clean.vcf
############################################

if [ -n "$RUNSAMVAR" ]; then
    if [ -z "$TASKBWA" ] || [ -z "$NODES_SAMVAR" ] || [ -z "$CPU_SAMVAR" ] || [ -z "$MEMORY_SAMVAR" ] || [ -z "$WALLTIME_SAMVAR" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB -r $ARMED -k $CONFIG -t $TASKBWA-$TASKSAMVAR -i $TASKBWA -e .$ASD.bam \
       -n $NODES_SAMVAR -c $CPU_SAMVAR -m $MEMORY_SAMVAR'G' -w $WALLTIME_SAMVAR \
       --command "${NGSANE_BASE}/mods/samSNPs.sh -k $CONFIG -f <FILE> -o $OUT/<DIR>/$TASKBWA-$TASKSAMVAR" \
	   --postcommand "${NGSANE_BASE}/mods/samSNPscollect.sh -k $CONFIG -f <FILE> -o $OUT/variant/$TASKBWA-$TASKSAMVAR-<DIR>"
fi

############################################
#   Mapping using HiCUP
############################################

if [ -n "$RUNHICUP" ]; then
    if [ -z "$TASKHICUP" ] || [ -z "$NODES_HICUP" ] || [ -z "$CPU_HICUP" ] || [ -z "$MEMORY_HICUP" ] || [ -z "$WALLTIME_HICUP" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASKHICUP -i fastq -e $READONE.$FASTQ -n $NODES_HICUP -c $CPU_HICUP \
    	-m $MEMORY_HICUP"G" -w $WALLTIME_HICUP \
        --command "${NGSANE_BASE}/mods/hicup.sh $HICUPADDPARM -k $CONFIG -t $CPU_HICUP -m $(expr $MEMORY_HICUP - 1 ) -f <FILE> -r $FASTA --digest '$HICUP_RENZYMES' -o $OUT/<DIR>/$TASKHICUP --fastqName <NAME>"
fi

############################################
#  Assessing HiC data with hiclib
#
# IN: $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/hiclib/*.hdf5
############################################

if [ -n "$RUNHICLIB" ]; then
    if [ -z "$TASKHICLIB" ] || [ -z "$NODES_HICLIB" ] || [ -z "$CPU_HICLIB" ] || [ -z "$MEMORY_HICLIB" ] || [ -z "$WALLTIME_HICLIB" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKHICLIB -i fastq -e $READONE.$FASTQ \
    	-n $NODES_HICLIB -c $CPU_HICLIB -m $MEMORY_HICLIB"G" -w $WALLTIME_HICLIB \
    	--postnodes $NODES_HICLIB_POSTCOMMAND --postcpu $CPU_HICLIB_POSTCOMMAND \
        --command "${NGSANE_BASE}/mods/hiclibMapping.sh $HICLIBADDPARM -k $CONFIG --threads $CPU_HICLIB --fastq <FILE> --enzymes '$HICLIB_RENZYMES' --outdir $OUT/<DIR>/$TASKHICLIB --fastqName <NAME>" \
        --postcommand "${NGSANE_BASE}/mods/hiclibCorrelate.sh $HICLIBADDPARM -f <FILE> -k $CONFIG --outdir $OUT/hiclib/$TASKHICLIB-<DIR>"

fi


############################################
#   Mapping using BWA
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bwa/*.$ASD.bam
############################################

if [ -n "$RUNMAPPINGBWA2" ]; then
    if [ -z "$TASKBWA" ] || [ -z "$NODES_BWA" ] || [ -z "$CPU_BWA" ] || [ -z "$MEMORY_BWA" ] || [ -z "$WALLTIME_BWA" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKBWA -i fastq -e $READONE.$FASTQ -n $NODES_BWA -c $CPU_BWA -m $MEMORY_BWA"G" -w $WALLTIME_BWA \
        --command "${NGSANE_BASE}/mods/bwa.sh $BWAADDPARM -k $CONFIG -t $CPU_BWA -m $(expr $MEMORY_BWA - 1 ) -f <FILE> -r $FASTA \
                -o $OUT/<DIR>/$TASKBWA --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR> \
                --fastqName $FASTQ -R $SEQREG"
fi

############################################
#   Mapping using Bowtie v1
#
# IN:$SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie/*.bam
############################################

if [ -n "$RUNMAPPINGBOWTIE" ]; then
    if [ -z "$TASKBOWTIE" ] || [ -z "$NODES_BOWTIE" ] || [ -z "$CPU_BOWTIE" ] || [ -z "$MEMORY_BOWTIE" ] || [ -z "$WALLTIME_BOWTIE" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi

    $QSUB $ARMED -k $CONFIG -t $TASKBOWTIE -i fastq -e $READONE.$FASTQ -n $NODES_BOWTIE -c $CPU_BOWTIE -m $MEMORY_BOWTIE"G" -w $WALLTIME_BOWTIE \
        --command "${NGSANE_BASE}/mods/bowtie.sh $BOWTIEADDPARM -k $CONFIG -t $CPU_BOWTIE -m $(expr $MEMORY_BOWTIE - 1 ) -f <FILE> -r $FASTA \
        -o $OUT/<DIR>/$TASKBOWTIE --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR> \
        --fastqName $FASTQ -R $SEQREG"
fi



############################################
#  Mapping with bowtie v2
#
# IN: $SOURCE/$dir/fastq/*read1.fastq
# OUT: $OUT/$dir/bowtie/*.bam
############################################
if [ -n "$RUNMAPPINGBOWTIE2" ]; then
    if [ -z "$TASKBOWTIE" ] || [ -z "$NODES_BOWTIE" ] || [ -z "$CPU_BOWTIE" ] || [ -z "$MEMORY_BOWTIE" ] || [ -z "$WALLTIME_BOWTIE" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -k $CONFIG -t $TASKBOWTIE -i fastq -e $READONE.$FASTQ -n $NODES_BOWTIE -c $CPU_BOWTIE -m $MEMORY_BOWTIE"G" -w $WALLTIME_BOWTIE \
	--command "${NGSANE_BASE}/mods/bowtie2.sh $BOWTIEADDPARM -k $CONFIG -t $CPU_BOWTIE -m $(expr $MEMORY_BOWTIE - 1 ) -f <FILE> -r $FASTA -o $OUT/<DIR>/$TASKBOWTIE \
        --rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi <DIR> --fastqName <NAME>"
fi

############################################
#  Analysis with homer
#
# IN: $SOURCE/$dir/bowtie/*.bam
# OUT: $OUT/$dir/homerhic/
############################################
if [ -n "$RUNHOMERHIC" ]; then
    if [ -z "$TASKHOMERHIC" ] || [ -z "$NODES_HOMERHIC" ] || [ -z "$CPU_HOMERHIC" ] || [ -z "$MEMORY_HOMERHIC" ] || [ -z "$WALLTIME_HOMERHIC" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
    
    $QSUB $ARMED -r -k $CONFIG -t $TASKHOMERHIC -i $TASKBWA -e $READONE.$ASD.bam -n $NODES_HOMERHIC -c $CPU_HOMERHIC -m $MEMORY_HOMERHIC"G" -w $WALLTIME_HOMERHIC \
	--command "${NGSANE_BASE}/mods/hicHomer.sh -k $CONFIG -t $CPU_HOMERHIC -f <FILE> -o $OUT/<DIR>/$TASKHOMERHIC"
fi

############################################
#  Creating normalized (wig) files with wiggler
#
# IN: $SOURCE/<DIR>/bowtie/*.bam
# OUT: $OUT/$dir/wiggler/
############################################
if [ -n "$RUNWIGGLER" ]; then
    if [ -z "$TASKWIGGLER" ] || [ -z "$NODES_WIGGLER" ] || [ -z "$CPU_WIGGLER" ] || [ -z "$MEMORY_WIGGLER" ] || [ -z "$WALLTIME_WIGGLER" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi

    $QSUB $ARMED -r -k $CONFIG -t $TASKWIGGLER -i $TASKBWA -e .$ASD.bam -n $NODES_WIGGLER -c $CPU_WIGGLER -m $MEMORY_WIGGLER"G" -w $WALLTIME_WIGGLER \
        --postcommand "${NGSANE_BASE}/mods/wiggler.sh -k $CONFIG -f <FILE> -m $MEMORY_WIGGLER -o $OUT/<DIR>/$TASKWIGGLER" 
fi

############################################
# downsampling
############################################
if [ -n "$RUNDOWNSAMPLING" ]; then

    echo "********** Downsampling"
    let ABB=$READNUMBER/100000
    TASKDOWNSAMPLE="sample"$ABB"K"
    for dir in ${DIR[@]}; do

      if [ ! -d $QOUT/$TASKDOWNSAMPLE ]; then mkdir -p $QOUT/$TASKDOWNSAMPLE; fi
      if [ ! -d $OUT/$dir/$TASKDOWNSAMPLE ]; then mkdir -p $OUT/$dir/$TASKDOWNSAMPLE; fi


      for f in $( ls $OUT/$dir/$TASKRCA/*$ASR.bam ); do
          n=`basename $f`
          NAME=$dir"_"$n
          echo $dir"/"$TASKRCA"/"$NAME" -> "$dir/$TASKDOWNSAMPLE

          if [ -e $QOUT/$TASKDOWNSAMPLE/$NAME.out ]; then rm $QOUT/$TASKDOWNSAMPLE/$NAME.out; fi

          if [ -n "$ARMED" ]; then
	  		 qsub $PRIORITY -j y -o $QOUT/$TASKDOWNSAMPLE/$NAME.out -cwd -b y \
	  		 -N $TASKDOWNSAMPLE"_"$NAME -l vf=4G \
			${NGSANE_BASE}/mods/downsample.sh -k ${NGSANE_BASE} -i $f -o $OUT/$dir/$TASKDOWNSAMPLE \
			-s $READNUMBER -r $FASTA

		  fi
      done
  done

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
	merge.sh ${NGSANE_BASE} $OUT/combined/mergeguide/lanes/$e $OUT/combined/bwa/ ${e/tmp/$ASD.bam} bam qout/merged/
    done
fi


########
# run tophat,
# already upgraded to take the config rather than the toolkit CAREFULL to move $OUT ->$OUTDIR 
# already cleanup in armed only
########
if [ -n "$RUNTOPHATCUFF" ]; then
    if [ -z "$TASKTOPHAT" ] || [ -z "$NODES_TOPHAT" ] || [ -z "$CPU_TOPHAT" ] || [ -z "$MEMORY_TOPHAT" ] || [ -z "$WALLTIME_TOPHAT" ]; then echo "[ERROR] Server misconfigured"; exit 1; fi
        
    echo "********* $TASKTOPHAT"
    
    #CPUS_TOPHAT=24
 
    # ensure directory is there
    if [ ! -d $QOUT/$TASKTOPHAT ]; then mkdir -p $QOUT/$TASKTOPHAT; fi
    for dir in ${DIR[@]}
      do
       
      if [ ! -d $OUT/$dir/$TASKTOPHAT ]; then mkdir -p $OUT/$dir/$TASKTOPHAT; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	name=${n/%$READONE.$FASTQ/}
	echo $name
	echo ">>>>>"$dir"/"$TASKTOPHAT"/"$n" ("$TASKCUFF")"
	
        #submit
	if [ -n "$ARMED" ]; then

            #cleanup old qouts
	    if [ -e $QOUT/$TASKTOPHAT/$dir'_'$name.out ]; then rm $QOUT/$TASKTOPHAT/$dir'_'$name.out; fi

	    $BINQSUB -j oe -o $QOUT/$TASKTOPHAT/$dir'_'$name.out -w $(pwd) -l walltime=$WALLTIME_TOPHAT \
		-N $TASKTOPHAT"_"$dir"_"$name -l $NODES_TOPHAT -l vmem=$MEMORY_TOPHAT"G" \
		-command "${NGSANE_BASE}/mods/tophatcuff.sh $TOPHATADDPARM -k $CONFIG -r $FASTA -f $f \
		-t $CPU_TOPHAT -o $OUT/$dir/$TASKTOPHAT/$name/ -a $REFSEQGTF"
	    #exit
	fi
      done


    done
fi


if [ -n "$RUNTOPHATCUFF2" ]; then
  $QSUB $ARMED -k $CONFIG -t $TASKTOPHAT -i fastq -e $READONE.$FASTQ -n $NODES_TOPHAT -c $CPU_TOPHAT -m $MEMORY_TOPHAT"G" -w $WALLTIME_TOPHAT \
        --command "${NGSANE_BASE}/mods/tophatcuff.sh $TOPHATADDPARM -k $CONFIG -f <FILE> \
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

    if [ ! -d $QOUT/$TASKRCA ]; then mkdir -p $QOUT/$TASKRCA; fi
    #ensure dirs are there...
    if [ ! -d $OUT/combined/$TASKRCA ]; then mkdir -p $OUT/combined/$TASKRCA; fi

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
		${NGSANE_BASE}/mods/reCalAln.sh ${NGSANE_BASE} $f $FASTA $DBROD \
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

    if [ ! -d $QOUT/$TASKRCA ]; then mkdir -p $QOUT/$TASKRCA; fi

    for dir in ${DIR[@]}
      do

     #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKRCA ]; then mkdir -p $OUT/$dir/$TASKRCA; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASD.bam}
	name=${n/%$READONE.$FASTQ/}
	echo ">>>>>"$dir$n2

	# wait on pipeline steps
	if [ -n "$RUNMAPPINGBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi
	# remove old pbs output
	if [ -e $QOUT/$TASKRCA/$dir'_'$name.out ]; then rm $QOUT/$TASKRCA/$dir'_'$name.*; fi
			
	#Sumit (ca. 80 min on AV 1297205 reads) -l h_rt=20:00:00
	if [ -n "$ARMED" ]; then
	    $BINQSUB -j oe -o $QOUT/$TASKRCA/$dir'_'$name'.out' -w $(pwd) -l $NODES_RECAL \
		-l vmem=$MEMORY_RECAL'G' -N $TASKRCA'_'$dir'_'$name -l walltime=$WALLTIME_RECAL \
		-command "${NGSANE_BASE}/mods/reCalAln.sh -k ${NGSANE_BASE} -f $OUT/$dir/$TASKBWA/$n2 -r $FASTA -d $DBROD \
		-o $OUT/$dir/$TASKRCA -t $CPU_RECAL $RECALADDPARAM"
	fi
	
      done
    done

fi


if [ -n "$RUNREALRECAL2" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $TASKRCA -i $TASKBWA/ -e .$ASD.bam \
        -n $NODES_RECAL -c $CPU_RECAL -m $MEMORY_RECAL"G" -w $WALLTIME_RECAL \
        --command "${NGSANE_BASE}/mods/reCalAln.sh $RECALADDPARAM -k $CONFIG -f <FILE> -r $FASTA -d $DBROD -o $OUT/<DIR>/$TASKRCA -t $CPU_RECAL"

fi

if [ -n "$RUNREALRECAL3" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $TASKRCA -i $TASKBWA/ -e .$ASD.bam \
        -n $NODES_RECAL -c $CPU_RECAL -m $MEMORY_RECAL"G" -w $WALLTIME_RECAL \
        --command "${NGSANE_BASE}/mods/reCalAln2.sh $RECALADDPARAM -k $CONFIG -f <FILE> -r $FASTA -d $DBROD -o $OUT/<DIR>/$TASKRCA -t $CPU_RECAL"

fi


############################################
# combine all into one bamfile
#
############################################

if [ -n "$mergeReCalbams" ]; then
    # make tmp files for first step combining
    ls combined/reCalAln/*.bam >combined/mergeguide/combineAll.txt
    
    # combine them in lane bam files
    ${NGSANE_BASE}/mods/merge.sh ${NGSANE_BASE} $OUT/combined/mergeguide/combineAll.txt $OUT/combined/ DISC1_all.bam bam qout/merged/
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

    if [ ! -d $QOUT/$TASKDOC ]; then mkdir -p $QOUT/$TASKDOC; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
	echo ">>>>>"$dir/$TASKRCA/$n2

	if [ ! -d $OUT/$dir/$TASKDOC ]; then mkdir -p $OUT/$dir/$TASKDOC; fi

	# remove old pbs output
	if [ -e $QOUT/$TASKDOC/$dir'_'$name'.out' ]; then rm $QOUT/$TASKDOC/$dir'_'$name'.out'; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$$RUNREALRECAL" ]; then HOLD="-hold_jid "$TASKRCA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -j y -o $QOUT/$TASKDOC/$dir'_'$name'.out' -cwd -b y -pe mpich $CPUS \
		-l mem_free=11G -l h_vmem=11G -l vf=500K -N $TASKDOC'_'$dir'_'$name $HOLD\
		${NGSANE_BASE}/mods/gatkDOC.sh -k ${NGSANE_BASE} -f $OUT/$dir/$TASKRCA/$n2 -r $FASTA \
		-o $OUT/$dir/$TASKDOC -t $CPUS $DOCADDPARAM
	fi

      done
    done

fi

if [ -n "$DEPTHOFCOVERAGE2" ]; then

    $QSUB $ARMED -r -k $CONFIG -t $TASKDOC -i $TASKRCA/ -e .$ASR.bam \
	-n $NODES_GATKDOC -c $CPU_GATKDOC -m $MEMORY_GATKDOC"G" -w $WALLTIME_GATKDOC \
	--command "${NGSANE_BASE}/mods/gatkDOC.sh $DOCADDPARAM -k ${NGSANE_BASE} -f <FILE> -r $FASTA -o $OUT/<DIR>/$TASKDOC -t $CPU_GATKDOC"

fi



############################################
# downsample
############################################

if [ -n "$DOWNSAMPLE" ]
then

    echo "********* $TASKDOWN"

    if [ ! -d $QOUT/$TASKDOWN ]; then mkdir -p $QOUT/$TASKDOWN; fi

    for dir in ${DIR[@]}
      do

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASD.bam}
	name=${n/%$READONE.$FASTQ/}
	echo ">>>>>"$dir$n2

	if [ ! -d $OUT/$dir/$TASKDOWN ]; then mkdir -p $OUT/$dir/$TASKDOWN; fi

	# remove old pbs output
	if [ -e $QOUT/$TASKDOWN/$dir"_"$n2"0.out" ]; then rm $QOUT/$TASKDOWN/$dir"_"$n2.*; fi

	#check if this is part of the pipe and jobsubmission needs to wait
	if [ -n "$mappingBWA" ]; then HOLD="-hold_jid "$TASKBWA"_"$dir"_"$name; fi

	#Submit
	if [ -n "$ARMED" ]; then
	    #downsample
	    python ${NGSANE_BASE}/tools/downsample.py -i $OUT/$dir/bwa/$n2 -o $OUT/$dir/$TASKDOWN/ -t downsample/$dir --region $SEQREG -w 500 -s 500 -q $QOUT/$TASKDOWN/$dir

	fi

      done
    done

# echo $( ls Hannibal_FC30MEJAAXX/downsample/*500d.bam ) > mergeBamfiles500d.tmp
# qsub -b y -cwd -j y -o qout/$COMBDIR/merged.asd.500d -pe mpich 2 /clusterdata/hiseq_apps/hiSeqInf/mods/merge.sh ${NGSANE_BASE} mergeBamfiles500d.tmp merged merged.$ASD.500d.bam bam

fi


############################################
# combine dindel vcfs
############################################

if [ -n "$DINDELcombineVCF" ]
then

    echo "********* $TASKDINS combine"

    if [ ! -d $QOUT/$TASKDINS ]; then mkdir -p $QOUT/$TASKDINS; fi
    if [ -e $OUT/dindelVCFmerge.tmp ]; then rm $OUT/dindelVCFmerge.tmp; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.ashrr.bam.dindel.VCF}
	echo "$OUT/$dir/dindelS/$n2" >> dindelVCFmerge.tmp
      done
    done

    # remove old pbs output
    if [ -e $QOUT/$TASKDINS/dindelVCFmerge.out ]; then rm $QOUT/$TASKDINS/dindelVCFmerge.out; fi

    #Submit
    if [ -n "$ARMED" ]; then
	qsub $PRIORITY -j y -o $QOUT/$TASKDINS/dindelVCFmerge.out -cwd -b y -l h_vmem=12G -N dindelVCFmerge \
	    ${NGSANE_BASE}/mods/merge.sh ${NGSANE_BASE} dindelVCFmerge.tmp $OUT/genotype dindelSeparate.vcf vcf
	    
    fi

fi


############################################
# combine dindel vcfs
############################################

if [ -n "$DINDELcombine" ]
then

    echo "********* $TASKDINS combine"

    if [ ! -d $QOUT/$TASKDINS ]; then mkdir -p $QOUT/$TASKDINS; fi
    if [ -e $OUT/dindelVCFmerge.tmp ]; then rm $OUT/dindelVCFmerge.tmp; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.ashrr.bam.dindel.VCF}
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
#   QC: ~/hiSeqInf/mods/dindelQC.sh qout/dindel/
############################################

if [ -n "$DINDELcallIndelsCombined" ]
then

    echo "********* $TASKDINC"

    if [ ! -d $QOUT/$TASKDINC ]; then mkdir -p $QOUT/$TASKDINC; fi
    if [ ! -d $COMBDIR/$TASKDINC ]; then mkdir -p $COMBDIR/$TASKDINC; fi
 
    for dir in ${DIR[@]}; do
	CANDIDATES=$CANDIDATES" "`echo $(ls $dir/dindel/*dindel.VCF)`
    done
    CANDIDATES=${CANDIDATES// /,}

    #Submit
    if [ -n "$ARMED" ]; then
	    ${NGSANE_BASE}/mods/dindelV2.sh ${NGSANE_BASE} $CANDIDATES $FASTA $COMBDIR/$TASKDINC/

    fi
fi



############################################
#   call indels with GATK
############################################

#TODO run them together with gatkIndelV2.sh

if [ -n "$GATKcallIndelsSeperate" ]
then

    echo "********* $TASKIND"
    if [ ! -d $QOUT/$TASKIND ]; then mkdir -p $QOUT/$TASKIND; fi
    if [ -e $QOUT/$TASKIND/ids.txt ]; then rm $QOUT/$TASKIND/ids.txt; fi
    if [ -e $QOUT/$TASKIND/sum.out ]; then rm $QOUT/$TASKIND/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKIND ]; then mkdir -p $OUT/$dir/$TASKIND; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq )
	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
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
		${NGSANE_BASE}/mods/gatkIndel.sh ${NGSANE_BASE} $OUT/$dir/$TASKRCA/$n2 $FASTA $DBROD \
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

    if [ ! -d $QOUT/$TASKIND ]; then mkdir -p $QOUT/$TASKIND; fi
    if [ -e $QOUT/$TASKIND/ids.txt ]; then rm $QOUT/$TASKIND/ids.txt; fi
    if [ ! -d $OUT/genotype ]; then mkdir -p $OUT/genotype; fi

    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')

    for dir in ${DIR[@]};do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.fastq );do
	n=`basename $f`
	echo $OUT/$dir/$TASKRCA/${n/%$READONE.$FASTQ/.$ASR.bam} >> $TASKIND"bamfiles.tmp"
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
	    ${NGSANE_BASE}/mods/gatkIndelV2.sh ${NGSANE_BASE} $TASKIND"bamfiles.tmp" $FASTA $DBROD \
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

    if [ ! -d $QOUT/$TASKVAR ]; then mkdir -p $QOUT/$TASKVAR; fi
    if [ ! -d $OUT/$TASKVAR ]; then mkdir -p $OUT/$TASKVAR; fi

    RUNVARS=""

    # one per folder
    if [ -n "$VARPERFOLDER" ]; then
	for dir in ${DIR[@]};do
	    RUNVARS=$RUNVARS" "$dir
	    if [ ! -d $OUT/$TASKVAR/$dir ]; then mkdir -p $OUT/$TASKVAR/$dir; fi
            if [ -e  $OUT/$TASKVAR/$dir/$TASKVAR"bamfiles.tmp" ]; then rm  $OUT/$TASKVAR/$dir/$TASKVAR"bamfiles.tmp"; fi
	    for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );do
		n=`basename $f`
		echo $OUT/$dir/$TASKRCA/${n/%$READONE.$FASTQ/.$ASR.bam} >> $OUT/$TASKVAR/$dir/$TASKVAR"bamfiles.tmp"
	    done
	done
    # one over all folders
    else
	NAME=$(echo ${DIR[@]}|sed 's/ /_/g')
	RUNVARS=$NAME
	if [ ! -d $OUT/$TASKVAR/$NAME ]; then mkdir -p $OUT/$TASKVAR/$NAME; fi
        if [ -e  $OUT/$TASKVAR/$NAME/$TASKVAR"bamfiles.tmp" ]; then rm  $OUT/$TASKVAR/$NAME/$TASKVAR"bamfiles.tmp"; fi
	for dir in ${DIR[@]};do
	    for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );do
		n=`basename $f`
		echo $OUT/$dir/$TASKRCA/${n/%$READONE.$FASTQ/.$ASR.bam} >> $OUT/$TASKVAR/$NAME/$TASKVAR"bamfiles.tmp"
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
	    #$BINQSUB -j oe -o $QOUT/$TASKVAR/$NAME.out -w $(pwd) -l $NODES_VAR $HOLD \
		#-l vmem=$MEMORY_VAR"G" -N $TASKVAR'_'$NAME -l walltime=$WALLTIME_VAR \
		#-command "${NGSANE_BASE}/mods/gatkSNPs.sh -k $CONFIG -i $OUT/$TASKVAR/$NAME/$TASKVAR'bamfiles.tmp' -t $CPU_VAR \
		#-r $FASTA -d $DBROD -o $OUT/$TASKVAR/$NAME -n $NAME -H $HAPMAPVCF -K $ONEKGVCF $VARADDPARAM"
    	    
	    $BINQSUB -a "$QSUBEXTRA" -k $CONFIG -m $MEMORY_VAR"G" -n $NODES_VAR -w $WALLTIME_VAR \
	    	-j $TASKVAR'_'$NAME -o $QOUT/$TASKVAR/$NAME.out \
			--command "${NGSANE_BASE}/mods/gatkSNPs.sh -k $CONFIG \
			-i $OUT/$TASKVAR/$NAME/$TASKVAR'bamfiles.tmp' -t $CPU_VAR \
			-r $FASTA -d $DBROD -o $OUT/$TASKVAR/$NAME -n $NAME \
			-H $HAPMAPVCF -K $ONEKGVCF $VARADDPARAM"


	fi
    done

    
fi

if [ -n "$RUNVARCALLS3" ]; then
    NAME=$(echo ${DIR[@]}|sed 's/ /_/g')
    $QSUB $ARMED -d -r -k $CONFIG -t $TASKVAR -i $TASKRCA/  -e .$ASR.bam -n $NODES_VAR \
        -c $CPU_VAR -m $MEMORY_VAR"G" -w $WALLTIME_VAR \
        --postcommand "${NGSANE_BASE}/mods/gatkSNPs2.sh -k $CONFIG \
                        -i <FILE> -t $CPU_VAR \
                        -r $FASTA -d $DBROD -o $OUT/$TASKVAR/$NAME -n $NAME \
                        -H $HAPMAPVCF -K $ONEKGVCF $VARADDPARAM"
fi




########
# run differental expression detection
########
if [ -n "$RUNANNOTATION" ]; then

    # ensure directory is there
    if [ ! -d $QOUT/$TASKANNOVAR ]; then mkdir -p $QOUT/$TASKANNOVAR; fi
    if [ ! -d $OUT/$TASKANNOVAR ]; then mkdir -p $OUT/$TASKANNOVAR; fi

    
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
	if [ ! -d $OUT/$TASKANNOVAR/$n ]; then mkdir -p $OUT/$TASKANNOVAR/$n; fi


        #check if this is part of the pipe and jobsubmission needs to wait
        if [ -n "$TASKVAR" ]; then HOLD="-hold_jid "$TASKVAR"_"$n; fi
	#echo "-hold_jid "$TASKVAR"_"$n
	HOLD="-hold_jid "$TASKVAR"_"$n

        #submit
	if [ -n "$ARMED" ]; then
	    qsub $PRIORITY -b y -cwd -j y -o $QOUT/$TASKANNOVAR/$n.out \
		-N $TASKANNOVAR'_'$n $HOLD\
		${NGSANE_BASE}/mods/annovar.sh -k ${NGSANE_BASE} -i1 $namesnp -i2 $namesnp2 -i3 $nameindel -r $FASTA \
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

    if [ ! -d $QOUT/$TASKSNP ]; then mkdir -p $QOUT/$TASKSNP; fi
    if [ -e $QOUT/$TASKSNP/ids.txt ]; then rm $QOUT/$TASKSNP/ids.txt; fi
    if [ -e $QOUT/$TASKSNP/sum.out ]; then rm $QOUT/$TASKSNP/sum.out; fi


    for dir in ${DIR[@]}
      do

      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKSNP ]; then mkdir -p $OUT/$dir/$TASKSNP; fi

      

      #for f in $( ls $SOURCE/$dir/aln2/*$ASR.bam )
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	n2=${n/%$READONE.$FASTQ/.$ASR.bam}
	name=${n/%$READONE.$FASTQ/}
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
		${NGSANE_BASE}/mods/gatkSNPs.sh $CONFIG $OUT/$dir/$TASKRCA/$n2 $FASTA $DBVCF \
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
    if [ ! -d $QOUT/$TASKCUFFDIFF ]; then mkdir -p $QOUT/$TASKCUFFDIFF; fi
    if [ ! -d $OUT/$TASKCUFFDIFF ]; then mkdir -p $OUT/$TASKCUFFDIFF; fi

    for dir in ${DIR[@]}; do
      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ );	do
	n=`basename $f`
	name=$name${n/%$READONE.$FASTQ/}","
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
            ${NGSANE_BASE}/mods/cuffdiff.sh -k ${NGSANE_BASE} -b $name -r $FASTA \
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

    if [ ! -d $QOUT/$TASKRRBS ]; then mkdir -p $QOUT/$TASKRRBS; fi

    for dir in ${DIR[@]}
      do
      
      #ensure dirs are there...
      if [ ! -d $OUT/$dir/$TASKRRBS ]; then mkdir -p $OUT/$dir/$TASKRRBS; fi

      for f in $( ls $SOURCE/fastq/$dir/*$READONE.$FASTQ )
	do
	n=`basename $f`
	name=${n/%$READONE.$FASTQ/}
	echo ">>>>>"$dir$n
		
	# remove old pbs output
	if [ -e $QOUT/$TASKRRBS/$dir'_'$name.out ]; then rm $QOUT/$TASKRRBS/$dir'_'$name.*; fi

	if [ -n "$ARMED" ]; then
	   qsub $PRIORITY -j y -o $QOUT/$TASKRRBS/$dir'_'$name'.out' -cwd -b y -pe mpich $CPUS \
		-l h_vmem=12G -N $TASKRRBS'_'$dir'_'$name \
		${NGSANE_BASE}/mods/rrbsmap.sh $RRBSMAPADDPARM -k ${NGSANE_BASE} -t $CPUS -f $f -r $FASTA -o $OUT/$dir/$TASKRRBS \
		--rgid $EXPID --rglb $LIBRARY --rgpl $PLATFORM --rgsi $dir --rgpu $FLOWCELL --fastqName $FASTQ 

	fi

      done
    done

fi


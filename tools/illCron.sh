
. /etc/profile.d/sge-binaries.sh

#PROGRAMS
. $HISEQINF/conf/header.sh
. $CONFIG

FORCE=$1

RUNS="/illumina/Runs/"
EMAIL="z.zhang5@uq.edu.au"
LOG="MyLogs"
SCRIPTS=$HISEQINF/mods

for r in $( ls $RUNS )
do

  # make mylog
  if [ ! -d $RUNS/$r/$LOG ]; then mkdir $RUNS/$r/$LOG; fi

  # find out how many "reads" the machine deposits
  if [ -e $RUNS/$r/RunInfo.xml ]; then 
      num=$( grep -w 'Read' $RUNS/$r/RunInfo.xml | wc -l )
  fi

  if [ -e $RUNS/$r/$LOG/ignore ]; then echo "-------> $r ignored"; continue; fi

  # look for HiSeq run completion flag
  if [ -e $RUNS/$r/Basecalling_Netcopy_complete_*$num.txt ]; then

      # check that this folder is not already analyzed
      filan=`ls $RUNS/$r/$LOG/fastqPrep.out* 2> /dev/null | wc -l`
      if [ $filan -eq 0 ] ; then

	  # if run is prepared start it
	  if [ -e $RUNS/$r/SampleSheet.csv ]; then
	      echo "-------> $r starting run"
	      if [ $num -lt 3 ]; then MULT="--singleplex"; echo "single plex"; fi
	      echo "$SCRIPTS/MakeFastQ.sh -w $RUNS/$r --fastq 1 -l $LOG $MULT" > $RUNS/$r/$LOG/call.txt
	      echo "starting analysis $r" | mail -s "[HiSeq] $r start run" $EMAIL
	      $SCRIPTS/MakeFastQ.sh -w $RUNS/$r --fastq 1 -l $LOG $MULT >>$RUNS/$r/$LOG/call.txt
	      
	      
	  # if unprepared and not already send email
	  else
	      echo "-------> $r action required"
	      if [ ! -e $RUNS/$r/$LOG/actionReq.txt ]; then
		  echo "action required $r finished but no SampleSheet.txt file there" | mail -s "[HiSeq] $r action required" $EMAIL
		  echo "run finished "`date` > $RUNS/$r/$LOG/actionReq.txt
	      fi
	  fi
      else
	  # if makeFastq run command wrote the finished flag check all files for errors
	  if [ -e  $RUNS/$r/$LOG/done.txt ]; then
	      errors=`grep "rror" $RUNS/$r/$LOG/*.out`
	      # if there is an error write error flag and send mail
	      if [ -n "$errors" ]; then
		  if [ ! -e $RUNS/$r/$LOG/errors.txt ]; then
		      echo -e "$errors errors in $r \n$anno" | mail -s "[HiSeq] $r error" $EMAIL
		      echo "$errors errors in $r "`date` >$RUNS/$r/$LOG/errors.txt
		  fi
		  echo "-------> $r errors; $anno"
	      else
		  if [ ! -e $RUNS/$r/$LOG/doneMail.txt ]; then
		      echo -e "finished in $r \n$anno" | mail -s "[HiSeq] $r done" $EMAIL
		      echo "done in $r "`date` >$RUNS/$r/$LOG/doneMail.txt
		  fi
		  echo "-------> $r **ready"
	      fi
	  else
	      echo "-------> $r running analysis"
	  fi
      fi
  else
      echo "-------> $r incomplete (Hiseq still running)"
  fi
done


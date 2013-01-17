#
# General template for submitting a job 
#


#INPUTS
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;; # location of the HiSeqInf repository
	-o | --origin )         shift; ORIGIN=$1 ;; # subfile in $SOURCE
	-e | --fileending )     shift; ENDING=$1 ;; # select source files of a specific ending
	-t | --task )           shift; TASK=$1 ;; # what to do
	-n | --nodes )          shift; NODES=$1;;
	-m | --memory )         shift; MEMORY=$1;;
	-w | --walltime )       shift; WALLTIME=$1;;
	-c | --command )        shift; COMMAND=$1;;
	-r | --reverse )        REV="1";;
	-d | --nodir )          NODIR="1";;
	-a | --armed )          ARMED="armed";;
	--keep )                KEEP="keep";;
	--direct )              DIRECT="direct";;
	--first )               FIRST="first";;
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. $HISEQINF/pbsTemp/header.sh
. $CONFIG

echo "********* $TASK"

if [ ! -d $QOUT/$TASK ]; then mkdir $QOUT/$TASK; fi

## Select files in dir to run
if [ ! -e $QOUT/$TASK/runnow.tmp ]; then
    for dir in ${DIR[@]}; do
      
      #ensure dirs are there...
      if [ ! -n "$NODIR" ]; then
	  if [ ! -d $OUT/$dir/$TASK ]; then mkdir $OUT/$dir/$TASK; fi
      fi
      # print out 
      if [ -n "$REV" ]; then
	  for f in $( ls $SOURCE/$dir/$ORIGIN/*$ENDING); do
              echo $f >> $QOUT/$TASK/runnow.tmp
          done
      else
	  for f in $( ls $SOURCE/$ORIGIN/$dir/*$ENDING); do
	      echo $f >> $QOUT/$TASK/runnow.tmp
	  done
      fi
  done
fi


for i in $(cat $QOUT/$TASK/runnow.tmp); do

    n=$(basename $i) 
    # ending : fastq/xx or xx/bwa
    dir=$(dirname $i | gawk '{n=split($1,arr,"/"); print arr[n]}')
    if [ -n "$REV" ]; then dir=$(dirname $i | gawk '{n=split($1,arr,"/"); print arr[n-1]}'); fi
    name=${n/$ENDING/}
    echo ">>>>>"$dir"/"$name
                

    COMMAND2=${COMMAND//<FILE>/$i}
    COMMAND2=${COMMAND2//<DIR>/$dir}
    COMMAND2=${COMMAND2//<NAME>/$name}

    echo $COMMAND2


    if [ -n "$DIRECT" ]; then eval $COMMAND2; fi

    if [ -n "$ARMED" ]; then

        # remove old pbs output
	if [ -e $QOUT/$TASK/$dir'_'$name.out ]; then rm $QOUT/$TASK/$dir'_'$name.*; fi

        $BINQSUB $QSUBEXTRA -j oe -o $QOUT/$TASK/$dir'_'$name'.out' -w $(pwd) -l $NODES \
            -l vmem=$MEMORY -N $TASK'_'$dir'_'$name -l walltime=$WALLTIME \
            -command "$COMMAND2"

	if [ -n "$FIRST" ]; then exit; fi

    fi
done

if [ ! -n "$KEEP" ]; then  rm $QOUT/$TASK/runnow.tmp ; fi

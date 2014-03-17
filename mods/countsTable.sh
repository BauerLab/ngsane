
#!/bin/bash -e

# Template for new mods
# Work through the TODOs and replace all [TEMPLATE...] patterns
# author: Fabian Buske
# date: November 2013


echo ">>>>> [TEMPLATE purpose]"
echo ">>>>> startdate "`date`
echo ">>>>> hostname "`hostname`
echo ">>>>> job_name "$JOB_NAME
echo ">>>>> job_id "$JOB_ID
echo ">>>>> $(basename $0) $*"

function usage {
echo -e "usage: $(basename $0) -k NGSANE -f INPUTFILE -o OUTDIR [OPTIONS]"
exit
}

if [ ! $# -gt 3 ]; then usage ; fi

#INPUTS                                                                                                           
while [ "$1" != "" ]; do
    case $1 in
        -k | --toolkit )        shift; CONFIG=$1 ;;     # location of the NGSANE repository                       
        -f | --file )           shift; FILES=$1 ;;  # input file                                                       
        -o | --outdir )         shift; OUTDIR=$1 ;;     # output dir                                                     
        --recover-from )        shift; RECOVERFROM=$1 ;; # attempt to recover from log file
        -h | --help )           usage ;;
        * )                     echo "don't understand "$1
    esac
    shift
done

#PROGRAMS
. $CONFIG
. ${NGSANE_BASE}/conf/header.sh
. $CONFIG

################################################################################
CHECKPOINT="programs"

module list
echo "PATH=$PATH"

echo -e "--NGSANE      --\n" $(trigger.sh -v 2>&1)
## [TODO] test and output versions of software utilized in this mod 

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="parameters"





echo "[NOTE] Files: $FILES"
OLDFS=$IFS
IFS=","
DATASETS=""
for f in $FILES; do
    # get basename of f

    n=${f/%.$ASD.bam/}
    FILE=${n/$TASK_TOPHAT/$TASK_HTSEQCOUNT}
    # get directory
    d=$(dirname $f)
    d=${d##*/}


    # add to dataset
    if [ -n "$FILE" ]; then 
    	DATASETS="${DATASETS[@]} ${FILE[@]}"
    fi
done
IFS=$OLDFS

echo "$DATASETS"
exit 1




# delete old bam files unless attempting to recover
if [ -z "$RECOVERFROM" ]; then
    ## TODO remove primary result files from pervious runs
fi

## TODO remove comments if paired/single library preps should be detected based on the READ identifier patterns
#if [ "$INPUTFILE" != "${INPUTFILE/$READONE/$READTWO}" ] && [ -e ${INPUTFILE/$READONE/$READTWO} ]; then
#    PAIRED="1"
#else
#    PAIRED="0"
#fi

## TODO remove comments if compressed status of input files should be detected
#ZCAT="zcat"
#if [[ ${INPUTFILE##*.} != "gz" ]]; then ZCAT="cat"; fi

# unique temp folder that should be used to store temporary files
THISTMP=$TMP"/"$(whoami)"/"$(echo $OUTDIR | md5sum | cut -d' ' -f1)
mkdir -p $THISTMP

echo -e "\n********* $CHECKPOINT\n"
################################################################################
CHECKPOINT="recall files from tape"

if [ -n "$DMGET" ]; then
	dmget -a $INPUTFILE
    dmget -a $OUTDIR/*
    # TODO add additional resources that are required and may need recovery from tape
fi
    
echo -e "\n********* $CHECKPOINT\n"


THISTMP="${OUTDIR}/tmp/"   
	
################################################################################
CHECKPOINT="Make tables of counts."

echo "[NOTE] Make tables of counts."
     
for GTF in  "GTF" "GTF_masked"

  do
     
      for MODE in "union" "intersection-strict" "intersection-nonempty"

         do     
 
            for ATTR in "gene_id" "transcript_id"
    
               do 
   
               [ -f ${THISTMP}/files.txt ] &&  rm ${THISTMP}/files.txt
               touch ${THISTMP}/files.txt
               
                for FILE in  ${DATASETS}${GTF}.${MODE}.${ATTR}
                  do
                  [ -f $FILE ] && echo $FILE "Found" >> ${THISTMP}/files.txt || echo "Not found" >> ${THISTMP}/files.txt
                  done
    
                    if grep -q "Not found" ${THISTMP}/files.txt
                         
                         then
                         echo "[NOTE] ${GTF}.${MODE}.${ATTR} - not found, skipping."            
                         else
                         echo "[NOTE] ${GTF}.${MODE}.${ATTR} - found, making count tables."  
              
                            [ -f ${THISTMP}/joinedfile.txt ] && rm ${THISTMP}/joinedfile.txt
							
								for i in ${FILES}${GTF}.${MODE}.${ATTR}

									do

										 if [ ! -f ${THISTMP}/joinedfile.txt ] 
										
										   then
										
										      cat ${i} > ${THISTMP}/joinedfile.txt 

										    else

											   cut -f 2 ${i} | paste ${THISTMP}/joinedfile.txt -  > ${THISTMP}/tmp.txt 
											
											   mv ${THISTMP}/tmp.txt ${THISTMP}/joinedfile.txt
											   
										fi
									done
	
							echo ${FILES}${GTF}.${MODE}.${ATTR} | sed 's/ /,/g' | sed "s/\/${GTF}.${MODE}.${ATTR}//g" > ${THISTMP}/tmp.txt
							
							awk '{print "gene," $0;}' ${THISTMP}/tmp.txt > ${THISTMP}/out.csv
	
	                        cat ${THISTMP}/joinedfile.txt | sed 's/\t/,/g' >> ${THISTMP}/out.csv

	                        mv ${THISTMP}/out.csv ${OUTDIR}/${GTF}.${MODE}.${ATTR}.csv
    
 
                          fi
 
            done
 
        done

	done

################################################################################
CHECKPOINT="cleanup."


  [ -f ${THISTMP}/joinedfile.txt ] && rm ${THISTMP}/joinedfile.txt
  
  [ -f ${THISTMP}/joinedfile.txt ] && rm ${THISTMP}/out.csv
  
  [ -f ${THISTMP}/joinedfile.txt ] && rm ${THISTMP}/joinedfile.txt
  
  [ -f ${THISTMP}/joinedfile.txt ] && rm ${THISTMP}/files.txt
  


echo -e "\n********* $CHECKPOINT\n"
################################################################################

echo ">>>>> alignment with TopHat - FINISHED"
echo ">>>>> enddate "`date`





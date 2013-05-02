#!/bin/bash

DIR=$1

missing=""


# check "_pos.txt" are there
for l in $(seq 8); do
    for t in `seq 1 8; seq 21 28; seq 41 48; seq 61 68`; do
	if [ ! -e $DIR/Data/Intensities/s"_"$l"_"`printf %04i $t`"_pos.txt" ]; then
	    echo "Error: "$DIR/Data/Intensities/s"_"$l"_"`printf %04i $t`"_pos.txt"
	fi
    done
done

# check ".filter" are there
for l in $(seq 8); do
    for t in `seq 1 8; seq 21 28; seq 41 48; seq 61 68`; do
	if [ ! -e $DIR/Data/Intensities/BaseCalls/s"_"$l"_"`printf %04i $t`".filter" ]; then
	    echo "Error: "$DIR/Data/Intensities/BaseCalls/s"_"$l"_"`printf %04i $t`".filter"
	fi
    done
done

missing=""

#find out how many cycles are there suppose to be
cycles=`grep "Cycles Incorporation" $DIR/RunInfo.xml | cut -d "\"" -f 2`

# check ".stats" and ".bcl"
for c in $(seq $cycles); do
    for l in $(seq 8); do
	for t in `seq 1 8; seq 21 28; seq 41 48; seq 61 68`; do
	    DIRL=L`printf %03i $l`
	    DCYC=C$c".1"
	    FILE=$DIR/Data/Intensities/BaseCalls/$DIRL/$DCYC/s"_"$l"_"$t
	    if [ ! -e $FILE.bcl ]; then
		echo "Missing: "$FILE.bcl
		missing=$missing","s"_"$l"_"`printf %04i $t`
	    fi
	    if [ ! -e $FILE.stats ]; then
		echo "Missing: "$FILE.stats
		missing=$missing","s"_"$l"_"`printf %04i $t`
	    fi
	done
    done
done

echo ">$missing"
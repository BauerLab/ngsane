PHR=$1
DIR=$2

echo `date`

for f in $( ls $DIR ); do
    r=`grep -c "$PHR" $DIR$f `
    if [ "$r" = "0" ]; then
	echo $f
    fi
done
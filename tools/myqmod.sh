
function usage {
echo -e "usage: $(basename $0) <start> <end> [TASK]

Hold and release jobs

required:
  start       first job to hold
  end         last job to hold

options for TASK:
  -h          hold
  -r          release
  -d          delete
"
exit
}

if [ ! $# -gt 0 ]; then usage ; fi


if [[ "$3" = "-h" ]]; then
echo "Hold jobs ? y|n"
read input
if [ "$input" != "y" ];then exit 0; fi
for ((i=$1; i<=$2;i++)); do
    echo "qhold $i"
    qhold $i
done
fi

if [[ "$3" = "-r" ]]; then
echo "Release jobs ? y|n"
read input
if [ "$input" != "y" ];then exit 0; fi
for ((i=$1; i<=$2;i++)); do
    echo "qrls $i"
    qrls $i
done
fi

if [[ "$3" = "-d" ]]; then
echo "Delete jobs ? y|n"
read input
if [ "$input" != "y" ];then exit 0; fi
for ((i=$1; i<=$2;i++)); do
    echo "qdel $i"
    qrls $i
done
fi
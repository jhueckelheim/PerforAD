KERNEL_INDEX=$1
MODE=$2
NUM_THREADS=$3

for i in {1..10}
do
    N=$((10**$i))
   eval /home/cs-kukr1/stencil_ad/PerforAD/generated/runner.sh $KERNEL_INDEX $MODE $N $NUM_THREADS
done




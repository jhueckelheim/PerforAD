COMMAND="$@"

for i in {0..8}
do
    THREADS=$((2**$i))
    C="OMP_NUM_THREADS=$THREADS $COMMAND"
    for j in {0..4}
    do
    echo $C
    eval $C
done
done

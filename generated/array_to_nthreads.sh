KERNEL_INDEX=$1
MODE=$2

NUM_THREADS=$((2**$SLURM_ARRAY_TASK_ID))

/home/cs-kukr1/stencil_ad/PerforAD/generated/vary-n.sh $KERNEL_INDEX $MODE $NUM_THREADS

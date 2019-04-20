KERNEL_INDEX=$1
MODE=$2
DOMAIN_SIZE=$3
NUM_THREADS=$4

export KMP_AFFINITY=scatter

KERNEL_LIST=(burgers1d driver1dwave driver3dwave)

KERNEL_NAME=${KERNEL_LIST[$KERNEL_INDEX]}
BASE_PATH=/home/cs-kukr1/stencil_ad/PerforAD/generated/
COMMAND="OMP_NUM_THREADS=${NUM_THREADS} ${BASE_PATH}${KERNEL_NAME} $DOMAIN_SIZE $MODE >> /home/cs-kukr1/stencil_ad/PerforAD/generated/results.txt 2>&1"


echo $COMMAND >> /home/cs-kukr1/stencil_ad/PerforAD/generated/results.txt

eval $COMMAND >> /home/cs-kukr1/stencil_ad/PerforAD/generated/results.txt
echo -e "\n***\n" >> /home/cs-kukr1/stencil_ad/PerforAD/generated/results.txt

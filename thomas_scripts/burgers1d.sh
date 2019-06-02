#!/bin/bash -l
# Batch script to run an OpenMP threaded job on Legion with the upgraded
# software stack under SGE.
# 1. Force bash as the executing shell.
#$ -S /bin/bash
# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=6:0:0
# 3. Request 1 gigabyte of RAM for each core/thread (must be an integer)
#$ -l mem=5G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=15G
# 5. Set the name of the job.
#$ -N stenciladburgers1d
# 6. Select 12 threads (the most possible on Legion).
#$ -pe smp 12
# 7. Set the working directory to somewhere in your scratch space.  This is
# a necessary step with the upgraded software stack as compute nodes cannot
# write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/mmm0334/Scratch/output
#$ -A Imperial_ESE
# 8. Run the application.

export OMP_PLACES=cores
export OMP_PROC_BIND=close
/shared/ucl/apps/numactl/2.0.12/bin/numactl -H
/shared/ucl/apps/numactl/2.0.12/bin/numactl -s
/shared/ucl/apps/numactl/2.0.12/bin/numactl --cpunodebind=0 --membind=0 /home/mmm0334/stencil_ad/PerforAD/generated/vary-threads.sh /home/mmm0334/stencil_ad/PerforAD/generated/burgers1d 100000000

#!/bin/bash
#SBATCH --job-name=pflrw_d3e2_L1206_t1_N64_EdS_CLPT_mKPunc
#SBATCH --partition sciama2.q
#SBATCH --time=170:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=15
#SBATCH --output=/users/munozr/mycodes/pflrw_thorns/slurm_out/%x-%j.out
#SBATCH --error=/users/munozr/mycodes/pflrw_thorns/slurm_out/%x-%j.err
#SBATCH --exclude=node110,node162,node163

now=$(date)
echo "Start time : $now"
echo "Job ID : $SLURM_JOB_ID"
echo "Job name : $SLURM_JOB_NAME"

echo "=========================================="
echo "Loading modules"
module purge
module load system/intel64
module load intel_comp/2019.2
module load openmpi/4.0.1

echo "=========================================="
echo "Environment:"
echo "---- Nodes"
echo "Nodes assigned : $SLURM_JOB_NODELIST"
echo "Number of nodes allocated : $SLURM_NNODES"
echo "Memory per node : $SLURM_MEM_PER_NODE"
echo "---- Tasks"
echo "Number of tasks : $SLURM_NTASKS"
echo "Number of tasks requested per node : $SLURM_NTASKS_PER_NODE"
echo "Number of tasks requested per core : $SLURM_NTASKS_PER_CORE"
echo "Number of tasks to be initiated on each node : $SLURM_TASKS_PER_NODE"
echo "---- CPUs"
echo "Number of CPUs per task : $SLURM_CPUS_PER_TASK"
echo "Number of CPUs on the allocated node : $SLURM_CPUS_ON_NODE"
echo "Count of processors available to the job on this node : $SLURM_JOB_CPUS_PER_NODE"
echo "Memory per CPU : $SLURM_MEM_PER_CPU"

echo "=========================================="
echo "Setting up directories"

echo "Create simulationg directory"
JOBPATH=/mnt/lustre2/ET_sims/$SLURM_JOB_NAME
mkdir -p $JOBPATH

# While output-i directory already exists, update SIMDIR
i=0
SIMDIR="${JOBPATH}/output-$(printf "%04d" $i)"
while [ -d "$SIMDIR" ]
do
  echo "$SIMDIR already exists"
  i=$((i + 1))
  SIMDIR="${JOBPATH}/output-$(printf "%04d" $i)"
done
echo "Creating $SIMDIR"
mkdir $SIMDIR
cd $SIMDIR

# If it's a restart, copy checkpoint file
if [ "$i" -gt 0 ]
then
  echo "Copy checkpoint file"
  mkdir $SIMDIR/$SLURM_JOB_NAME
  PREVSIM="${JOBPATH}/output-$(printf "%04d" $((i - 1)))"
  cp $PREVSIM/$SLURM_JOB_NAME/checkpoint.chkpt.it_699909.file_*.h5 $SIMDIR/$SLURM_JOB_NAME/
fi

echo "Copy parameter file"
cp /users/munozr/mycodes/pflrw_thorns/par/$SLURM_JOB_NAME.par $SIMDIR/$SLURM_JOB_NAME.par

echo "Create stdout and stderr soft links"
STDNAME=${SLURM_JOB_NAME}-${SLURM_JOB_ID}
ln -s /users/munozr/mycodes/pflrw_thorns/slurm_out/$STDNAME.out $SIMDIR/$STDNAME.out
ln -s /users/munozr/mycodes/pflrw_thorns/slurm_out/$STDNAME.err $SIMDIR/$STDNAME.err

echo "=========================================="
echo "Run simulation"
mpirun -n $SLURM_NNODES /users/munozr/code/2020_05/Cactus/exe/cactus_sim $SIMDIR/$SLURM_JOB_NAME.par

now=$(date)
echo "End time : $now"    
    

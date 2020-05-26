#PBS -lselect=1:ncpus=32:mem=50gb
#PBS -lwalltime=24:00:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
module load mpi
mkdir $PBS_O_WORKDIR/output/potentialData
mpiexec -n 32 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/potentialTest -s 64 -p $PBS_O_WORKDIR/output/potentialData -i $PBS_O_WORKDIR/singleMonopoleData -n 4 -m 8
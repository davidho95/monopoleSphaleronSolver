#PBS -lselect=1:ncpus=4:mem=1gb
#PBS -lwalltime=00:30:00

export LD_LIBRARY_PATH=/apps/gcc/8.2.0/lib64:$LD_LIBRARY_PATH
module load mpi
mkdir $PBS_O_WORKDIR/singleMonopoleDataTest
mpiexec -n 4 /rds/general/user/dlh16/home/monopoleSphaleronSolver/bin/generateSingleMonopole -s 16 -p $PBS_O_WORKDIR/singleMonopoleDataTest -n 2 -m 2
cp -r $PBS_O_WORKDIR/singleMonopoleDataTest /rds/general/user/dlh16/home/monopoleSphaleronSolver/output
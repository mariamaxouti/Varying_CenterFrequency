 #!/bin/sh
 #PBS -l walltime=10:00:00
 #PBS -l select=1:ncpus=64:mem=2000gb

module load matlab/R2020b

cd $PBS_O_WORKDIR

## execution of program
matlab -nosplash -nodisplay -nodesktop -singleCompThread -logfile myTestRun_output.log -r "Iter_TR"

mkdir $PBS_O_WORKDIR/$PBS_JOBID
cp * $PBS_O_WORKDIR/$PBS_JOBID



#!/bin/bash -l
#SBATCH -N 2         #Use 2 nodes
#SBATCH -t 00:3:00  #Set 30 minute time limit
#SBATCH -q debug   #Submit to the regular QOS
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes
source /global/homes/t/tyounkin/atomIPS/atom-install-cori/GITR/env.edison.sh
source /global/homes/t/tyounkin/atomIPS/atom-install-cori/GITR/ftridyn/env.cori.sh
srun -n 8 ~/code/python2.7.14/bin/python ftMPI.py 

export OMPI_ROOT=/project/projectdirs/atom/users/elwasif/ompi/install_4.0
export PATH=$OMPI_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$OMPI_ROOT/lib64:$OMPI_ROOT/lib:$LD_LIBRARY_PATH
export MANPATH=$OMPI_ROOT/share/man:$MANPATH

srun -n $SLURM_NNODES --ntasks-per-node=1 hostname > .node_names.$$
for n in $(cat .node_names.$$); do echo "$n slots=24" >> .hostfile.$$ ; done
orte-dvm --hostfile ./.hostfile.$$ &

# MDNafion
All code related to molecular dynamics simulations of Nafion and p-TPN1 composites

***How to run in parallel with MPI***

```bash
mpirun -n Np lmp_mpi -in in.file
```


where `Np` is the number of MPI processes to be used and `in.file` is the name of LAMMPS input file.

`lmp_mpi` is the LAMMPS exectuable. The path to executable needs to either be full or be added to 
user's path as

```bash
export PATH=/home/user/lammps/builds/executable:$PATH
```

Path can be anything, elaborated for better demonstration.


***List of packages to install with LAMMPS***

QEQ

MISC

KSPACE

MOLECULE

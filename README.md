# APM
INF560_Project_APM

Different versions are kept in different branches with corresponded name. Compilation is done by make command.


Usage: 
seq branch: ./apm approximation_factor dna_database pattern1 pattern2 ...
MPI branch: salloc -N 1 -n 4 mpirun ./apm approximation_factor dna_database pattern1 pattern2 (n>=2)
openMP branch: salloc -n 4 mpirun ./apm approximation_factor dna_database pattern1 pattern2 
hybrid branch: salloc -N 1 -n 4 mpirun ./apm approximation_factor dna_database pattern1 pattern2 (n>=2)

Different versions are kept in different branches with corresponded name.

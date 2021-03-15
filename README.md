# APM
INF560_Project_APM

Different versions are kept in different branches with corresponded name. Compilation is done by make command.


## Usage: 

### seq branch: 
./apm approximation_factor dna_database pattern1 pattern2 ...

### MPI branch: 
salloc -N 1 -n 4 mpirun ./apm approximation_factor dna_database pattern1 pattern2 

salloc -N 1 -n 9 mpirun ./apm 2 ./dna/chr9.fa AAAA AAAA AAAA AAAA AAAA AAAA AAAA AAAA

### openMP branch: 
salloc -n 1 ./apm approximation_factor dna_database pattern1 pattern2 

salloc -n 1 ./apm 2 ./dna/chr9.fa AAAA AAAA AAAA AAAA AAAA AAAA AAAA AAAA
### hybrid branch: 
salloc -N 1 -n 4 mpirun ./apm approximation_factor dna_database pattern1 pattern2 

salloc -N 1 -n 9 mpirun ./apm 2 ./dna/chr9.fa AAAA AAAA AAAA AAAA AAAA AAAA AAAA AAAA

### gpu branch: 
./test approximation_factor dna_database pattern1 pattern2

./test 2 ./dna/chr9.fa AAAA AAAA AAAA AAAA AAAA AAAA AAAA AAAA


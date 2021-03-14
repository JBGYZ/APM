#include <cuda_runtime.h>
#include <cuda.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>

#include <iostream>
#include <memory>
#include <cassert>

int __device__ min3(int a, int b, int c) {
    return ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)));
}


int __device__ levenshtein_d(char *s1, char *s2, int len, int *column) {
    unsigned int x, y, lastdiag, olddiag;

    for (y = 1; y <= len; y++) {
        column[y] = y;
    }
    for (x = 1; x <= len; x++) {
        column[0] = x;
        lastdiag = x - 1;
        for (y = 1; y <= len; y++) {
            olddiag = column[y];
            column[y] = min3(
                    column[y] + 1,
                    column[y - 1] + 1,
                    lastdiag + (s1[y - 1] == s2[x - 1] ? 0 : 1)
            );
            lastdiag = olddiag;

        }
    }
    return (column[len]);
}

void __global__ kernelwrapper(int* d_n_matches, char * d_buf, char * d_pattern, int i, int size_pattern, int offset, int n_bytes, int approx_factor){

    /* Traverse the input data up to the end of the file */
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int distance = 0 ;
    int size ;

    size = size_pattern ;
    int* columns = (int *) malloc((size_pattern + 1) * sizeof(int));
    while (j < n_bytes) {
        if (n_bytes - j < size_pattern ){
            size = n_bytes - j ;
        }

        distance = levenshtein_d(d_pattern + offset, &d_buf[j], size, columns ) ;
        if ( distance <= approx_factor) {
            atomicAdd(&d_n_matches[i], 1);
        }

        j += stride;
    }
    free(columns);

}



#define APM_DEBUG 0

char * 
read_input_file( char * filename, int * size )
{
    char * buf ;
    off_t fsize;
    int fd = 0 ;
    int n_bytes = 1 ;

    /* Open the text file */
    fd = open( filename, O_RDONLY ) ;
    if ( fd == -1 ) 
    {
        fprintf( stderr, "Unable to open the text file <%s>\n", filename ) ;
        return NULL ;
    }


    /* Get the number of characters in the textfile */
    fsize = lseek(fd, 0, SEEK_END);
    if ( fsize == -1 )
    {
        fprintf( stderr, "Unable to lseek to the end\n" ) ;
        return NULL ;
    }

#if APM_DEBUG
    printf( "File length: %lld\n", fsize ) ;
#endif

    /* Go back to the beginning of the input file */
    if ( lseek(fd, 0, SEEK_SET) == -1 ) 
    {
        fprintf( stderr, "Unable to lseek to start\n" ) ;
        return NULL ;
    }

    /* Allocate data to copy the target text */
    buf = (char *)malloc( fsize * sizeof ( char ) ) ;
    if ( buf == NULL ) 
    {
        fprintf( stderr, "Unable to allocate %lld byte(s) for main array\n",
                fsize ) ;
        return NULL ;
    }

    n_bytes = read( fd, buf, fsize ) ;
    if ( n_bytes != fsize ) 
    {
        fprintf( stderr, 
                "Unable to copy %lld byte(s) from text file (%d byte(s) copied)\n",
                fsize, n_bytes) ;
        return NULL ;
    }

#if APM_DEBUG
    printf( "Number of read bytes: %d\n", n_bytes ) ;
#endif

    *size = n_bytes ;


    close( fd ) ;


    return buf ;
}

int 
main( int argc, char ** argv )
{
  char ** pattern ;
  char * filename ;
  int approx_factor = 0 ;
  int nb_patterns = 0 ;
  int i, j ;
  char * buf ;
  struct timeval t1, t2, t3;
  double duration ;
  int n_bytes ;
  int * n_matches ;

  /* Check number of arguments */
  if ( argc < 4 ) 
  {
    printf( "Usage: %s approximation_factor "
            "dna_database pattern1 pattern2 ...\n", 
            argv[0] ) ;
    return 1 ;
  }

  /* Get the distance factor */
  approx_factor = atoi( argv[1] ) ;

  /* Grab the filename containing the target text */
  filename = argv[2] ;

  /* Get the number of patterns that the user wants to search for */
  nb_patterns = argc - 3 ;

  /* Fill the pattern array */
  pattern = (char **)malloc( nb_patterns * sizeof( char * ) ) ;
  if ( pattern == NULL ) 
  {
      fprintf( stderr, 
              "Unable to allocate array of pattern of size %d\n", 
              nb_patterns ) ;
      return 1 ;
  }

  /* Grab the patterns */
  for ( i = 0 ; i < nb_patterns ; i++ ) 
  {
      int l ;

      l = strlen(argv[i+3]) ;
      if ( l <= 0 ) 
      {
          fprintf( stderr, "Error while parsing argument %d\n", i+3 ) ;
          return 1 ;
      }

      pattern[i] = (char *)malloc( (l+1) * sizeof( char ) ) ;
      if ( pattern[i] == NULL ) 
      {
          fprintf( stderr, "Unable to allocate string of size %d\n", l ) ;
          return 1 ;
      }

      strncpy( pattern[i], argv[i+3], (l+1) ) ;
  }


  printf( "Approximate Pattern Mathing: "
          "looking for %d pattern(s) in file %s w/ distance of %d\n", 
          nb_patterns, filename, approx_factor ) ;

  buf = read_input_file( filename, &n_bytes ) ;
  if ( buf == NULL )
  {
      return 1 ;
  }

  /* Allocate the array of matches */
  n_matches = (int *)malloc( nb_patterns * sizeof( int ) ) ;
  if ( n_matches == NULL )
  {
      fprintf( stderr, "Error: unable to allocate memory for %ldB\n",
              nb_patterns * sizeof( int ) ) ;
      return 1 ;
  }


  /* Matching process takes place in GPU */

    int* d_n_matches;
    char * d_pattern;
    char* d_buf;
    int* offset = (int *)malloc( nb_patterns * sizeof( int ) ) ;
    int* lens = (int *)malloc( nb_patterns * sizeof( int ) ) ;
    int sum_lens;
    lens[0] = strlen(pattern[0]);
    offset[0] = 0;
    sum_lens = lens[0];
    for (i = 1; i < nb_patterns; i++) {
        offset[i] = offset[i-1] + lens[i-1];
        lens[i] = strlen(pattern[i]);
        sum_lens += lens[i];
    }
    char* concat_patterns = (char*) malloc( sum_lens * sizeof( char ) ) ;
    for (i = 0; i < nb_patterns; i++) {
        strcpy (concat_patterns + offset[i], pattern[i]);
    }
    
    gettimeofday(&t1, NULL);
    cudaMalloc((void **)&d_n_matches, nb_patterns*sizeof(int));
    cudaMalloc((void **)&d_pattern, sum_lens*sizeof(char));
    cudaMalloc((void **)&d_buf, n_bytes);
    cudaMemcpy(d_pattern, concat_patterns, sum_lens*sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_buf, buf, n_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_n_matches, n_matches, nb_patterns*sizeof(int), cudaMemcpyHostToDevice);

    gettimeofday(&t2, NULL);

    duration = (t2.tv_sec -t1.tv_sec)+((t2.tv_usec-t1.tv_usec)/1e6);

    printf( "Memory copy from host to device done in %lf s \n", duration) ;

    int Dg = 4;
    int Db = 256;
    for (i = 0; i < nb_patterns; i++) {
        kernelwrapper<<<Dg,Db>>>(d_n_matches, d_buf, d_pattern, i, lens[i], offset[i], n_bytes, approx_factor);
    }
    
    cudaMemcpy(n_matches, d_n_matches, nb_patterns*sizeof(int), cudaMemcpyDeviceToHost);

    gettimeofday(&t3, NULL);

    duration = (t3.tv_sec -t2.tv_sec)+((t3.tv_usec-t2.tv_usec)/1e6);

    printf( "Calculation on GPU done in %lf s \n", duration) ;


  for ( i = 0 ; i < nb_patterns ; i++ )
  {
      printf( "Number of matches for pattern <%s>: %d\n", 
              pattern[i], n_matches[i] ) ;
  }

  return 0 ;
}
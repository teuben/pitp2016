#include <stdio.h>
#include "omp.h"


void main()
{
  omp_set_num_threads(2);

// #pragma omp parallel num_threads(4) 
#pragma omp parallel
  //int nt =  omp_get_num_threads();
  //printf("Found %d num_threads\n",nt);
  {
    // int ID = 0;
    int ID = omp_get_thread_num();
    
    printf(" Hello(%d) ",ID);
    printf(" World(%d) \n",ID);
  }
}

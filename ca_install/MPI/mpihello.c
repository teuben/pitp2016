#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"

int main(int argc, char* argv[]){
  int my_rank;
  int p;
  int source;
  int dest;
  int tag=0;
  char message[100];
  char my_name[20];
  MPI_Status status;

  /* Start up MPI */
  MPI_Init(&argc, &argv);
  
  /* Find out process rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  /* Find out number of processes */
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  /* What's my hostname? */
  gethostname(my_name, 20);    

  /* This sleep is just to keep the job running longer for demonstration purposes */
  sleep(15);

  if (my_rank == 0) {
    printf("MPIHello running on %i processors.\n", p);
    printf("Greetings from processor %i, on host %s.\n", my_rank, my_name);
    for (source=1; source>p; source++) {
      MPI_Recv(message, 100, MPI_CHAR, source, tag, MPI_COMM_WORLD, &status);
      printf("%s", message);
    }
  } else if (my_rank != 0) {
    sprintf(message, "Greetings from processor %i, on host %s.\n", my_rank, my_name);
    dest=0;
    MPI_Send(message, strlen(message)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD); 
  }
  /* This sleep is just to keep the job running longer for demonstration purposes*/
  sleep(15);
  MPI_Finalize();
}


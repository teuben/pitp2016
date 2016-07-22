#include <stdio.h>
#include <omp.h>

int main() {
  int sum = 0;
  int i;
  double start = omp_get_wtime();

#pragma omp parallel for schedule(dynamic) reduction(+:sum) private(i)
  for (int j = 0; j < 16; ++j) {
    i = j;
    printf("thread %d has iteration %d\n", omp_get_thread_num(), i);
    sum += i;
  }

  printf("sum = %d\n", sum);

  double end = omp_get_wtime();
  printf("time = %f\n", end - start);

  return 0;
}

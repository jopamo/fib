#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <gmp.h>
#include <stdbool.h>

#define MAX_N 10000

void printUsage(const char *programName) {
  printf("Usage: %s [OPTIONS]\n", programName);
  printf("Example: %s -B -L 40 -H 55 -T 35\n", programName);
  printf("Example: %s -C -N 50\n\n", programName);

  printf("Options:\n");
  printf("  -B, --benchmark   Run benchmarking. If set, allows the following options:\n");
  printf("  -L, --low <value>   Set the low bound for binary search. (Requires -B)\n");
  printf("  -H, --high <value>  Set the high bound for binary search. (Requires -B)\n");
  printf("  -T, --time <value>  Set the target time (in seconds) for the search. (Requires -B)\n");
  printf("  -C, --compare   Run comparison of dynamic and iterative Fibonacci methods for a specified 'n'.\n");
  printf("  -S, --slow    Include recursive method in comparison. (Use with -C)\n");
  printf("  -N, --nvalue <value> Specify the 'n' value for Fibonacci calculation. (Required for -C)\n");
  printf("  -h, --help    Show this help message and exit.\n\n");

  printf("Description:\n");
  printf("  This program calculates the Fibonacci number for a given value using different methods.\n");
  printf("  The benchmark (-B) option runs a binary search to find an 'n' value for which the recursive\n");
  printf("  Fibonacci calculation takes approximately the specified time. The low, high, and time values\n");
  printf("  are used to control this search. The compare (-C) option compares the results and execution\n");
  printf("  times of dynamic and iterative methods for a specified 'n' value, and includes recursive\n");
  printf("  method if -S is specified. Without -B or -C, the program simply performs the iterative approach.\n");
}

void fibonacci_aux(mpz_t result, mpz_t a, mpz_t b, int n, mpz_t memo[], bool initialized[]) {
  if (n == 0) {
    mpz_set(result, a);
    return;
  }

  if (initialized[n]) {
    mpz_set(result, memo[n]);
    return;
  }

  mpz_t next_a, next_b;
  mpz_init_set(next_a, b);
  mpz_init(next_b);
  mpz_add(next_b, a, b);

  if (n < MAX_N) {
    mpz_set(memo[n], next_b);
    initialized[n] = true;
  }

  fibonacci_aux(result, next_a, next_b, n - 1, memo, initialized);

  mpz_clear(next_a);
  mpz_clear(next_b);
}

void fibonacci_recursive(mpz_t result, int n) {
  static mpz_t memo[MAX_N];
  static bool initialized[MAX_N] = {false};

  static bool is_first_call = true;
  if (is_first_call) {
    for (int i = 0; i < MAX_N; ++i) {
      initialized[i] = false;
    }
    is_first_call = false;
  }

  mpz_t a, b;
  mpz_init_set_ui(a, 0);
  mpz_init_set_ui(b, 1);

  fibonacci_aux(result, a, b, n, memo, initialized);

  mpz_clear(a);
  mpz_clear(b);
}

int calcMaxFib() {
  const double phi = (1 + sqrt(5)) / 2;
  return (int) floor(log(ULLONG_MAX * sqrt(5)) / log(phi));
}

typedef struct {
  mpz_t a, b, c, d;
} Matrix2x2;

void matrixInit(Matrix2x2 *m) {
  mpz_init(m->a);
  mpz_init(m->b);
  mpz_init(m->c);
  mpz_init(m->d);
}

void matrixSetIdentity(Matrix2x2 *m) {
  mpz_set_ui(m->a, 1);
  mpz_set_ui(m->b, 0);
  mpz_set_ui(m->c, 0);
  mpz_set_ui(m->d, 1);
}

void matrixCopy(Matrix2x2 *dest, Matrix2x2 *src) {
  mpz_set(dest->a, src->a);
  mpz_set(dest->b, src->b);
  mpz_set(dest->c, src->c);
  mpz_set(dest->d, src->d);
}

void matrixMultiply(Matrix2x2 *result, Matrix2x2 *m1, Matrix2x2 *m2) {
  Matrix2x2 temp;
  matrixInit(&temp);

  mpz_mul(temp.a, m1->a, m2->a);
  mpz_addmul(temp.a, m1->b, m2->c);

  mpz_mul(temp.b, m1->a, m2->b);
  mpz_addmul(temp.b, m1->b, m2->d);

  mpz_mul(temp.c, m1->c, m2->a);
  mpz_addmul(temp.c, m1->d, m2->c);

  mpz_mul(temp.d, m1->c, m2->b);
  mpz_addmul(temp.d, m1->d, m2->d);

  matrixCopy(result, &temp);

  mpz_clear(temp.a);
  mpz_clear(temp.b);
  mpz_clear(temp.c);
  mpz_clear(temp.d);
}

void matrixPower(Matrix2x2 *result, Matrix2x2 *m, int n) {
  Matrix2x2 temp;
  matrixInit(&temp);
  matrixSetIdentity(&temp);

  if (n == 0) {
    matrixCopy(result, &temp);
    return;
  } else if (n == 1) {
    matrixCopy(result, m);
    return;
  }

  matrixPower(result, m, n / 2);
  matrixMultiply(&temp, result, result);

  if (n % 2 != 0) {
    matrixMultiply(result, &temp, m);
  }
  else {
    matrixCopy(result, &temp);
  }

  mpz_clear(temp.a);
  mpz_clear(temp.b);
  mpz_clear(temp.c);
  mpz_clear(temp.d);
}

void fibonacci_matrix(mpz_t result, int n) {
  if (n == 0) {
    mpz_set_ui(result, 0);
    return;
  }

  Matrix2x2 base, res;
  matrixInit(&base);
  matrixInit(&res);

  mpz_set_ui(base.a, 0);
  mpz_set_ui(base.b, 1);
  mpz_set_ui(base.c, 1);
  mpz_set_ui(base.d, 1);

  matrixPower(&res, &base, n);

  mpz_set(result, res.b);

  mpz_clear(base.a);
  mpz_clear(base.b);
  mpz_clear(base.c);
  mpz_clear(base.d);
  mpz_clear(res.a);
  mpz_clear(res.b);
  mpz_clear(res.c);
  mpz_clear(res.d);
}

int parseIntArg(char *arg) {
  errno = 0;
  char *endptr;
  long val = strtol(arg, &endptr, 10);

  if (errno != 0 || *endptr != '\0' || val < 0 || val > INT_MAX) {
    fprintf(stderr, "Invalid integer value: %s\n", arg);
    exit(1);
  }

  return (int)val;
}

void processArgs(int argc, char *argv[], int *binary_low, int *binary_high, int *target_time_sec, int *benchmark_flag, int *compare_flag, int *slow_flag, long long unsigned int *n_value) {
  *benchmark_flag = 0;
  *compare_flag = 1;
  *slow_flag = 1;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-B") == 0 || strcmp(argv[i], "--benchmark") == 0) {
      if (*compare_flag) {
        fprintf(stderr, "Error: -B contains all functionality of -C\n");
        exit(1);
      }
      *benchmark_flag = 1;
    }
    else if (strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--compare") == 0) {
      if (*benchmark_flag) {
        fprintf(stderr, "Error: -B contains all functionality of -C\n");
        exit(1);
      }
      *compare_flag = 1;
    }
    else if (strcmp(argv[i], "-S") == 0 || strcmp(argv[i], "--slow") == 0) {
      *slow_flag = 1;
      if (!*compare_flag) {
        fprintf(stderr, "Error: -S (--slow) requires -C (--compare).\n");
        exit(1);
      }
    }
    else if ((strcmp(argv[i], "-L") == 0 || strcmp(argv[i], "--low") == 0) && i + 1 < argc) {
      if (!*benchmark_flag) {
        fprintf(stderr, "Error: -L (--low) requires -B (--benchmark).\n");
        exit(1);
      }
      *binary_low = parseIntArg(argv[++i]);
    }
    else if ((strcmp(argv[i], "-H") == 0 || strcmp(argv[i], "--high") == 0) && i + 1 < argc) {
      if (!*benchmark_flag) {
        fprintf(stderr, "Error: -H(--high) requires -B(--benchmark).\n");
        exit(1);
      }
      *binary_high = parseIntArg(argv[++i]);
    }
    else if ((strcmp(argv[i], "-T") == 0 || strcmp(argv[i], "--time") == 0) && i + 1 < argc) {
      if (!*benchmark_flag) {
        fprintf(stderr, "Error: -T(--time) requires -B(--benchmark).\n");
        exit(1);
      }
      *target_time_sec = parseIntArg(argv[++i]);
    }
    else if ((strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--nvalue") == 0) && i + 1 < argc) {
      *n_value = parseIntArg(argv[++i]);
    }
    else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      printUsage(argv[0]);
      exit(1);
    }
  }
}

double measureTime(void (*fibFunc)(mpz_t, int), int n, mpz_t result) {
  struct timespec start, end;

  clock_gettime(CLOCK_MONOTONIC, &start);
  fibFunc(result, n);
  clock_gettime(CLOCK_MONOTONIC, &end);

  double time_taken = end.tv_sec - start.tv_sec;
  time_taken += (end.tv_nsec - start.tv_nsec) / 1e9;

  return time_taken;
}

void printTime(double timeInSeconds) {
  if (timeInSeconds >= 1.0) {
    // seconds
    printf("%f seconds\n", timeInSeconds);
  }
  // 1e-3 represents 0.001 seconds (1 millisecond)
  else if (timeInSeconds >= 1e-3) {
    printf("%f milliseconds\n", timeInSeconds * 1e3);
  }
  // 1e-6 represents 0.000001 seconds (1 microsecond)
  else if (timeInSeconds >= 1e-6) {
    printf("%f microseconds\n", timeInSeconds * 1e6);
  }
  // Assume nanoseconds at this point (1 sec = 1,000,000,000 ns)
  else {
    printf("%f nanoseconds\n", timeInSeconds * 1e9);
  }
}

void fibonacci_dynamic(mpz_t result, int n) {
  // Handle the base cases
  if (n <= 1) {
    mpz_set_ui(result, n);
    return;
  }

  mpz_t fib_prev, fib_curr, fib_next;
  mpz_init_set_ui(fib_prev, 0); // fib[0] = 0
  mpz_init_set_ui(fib_curr, 1); // fib[1] = 1
  mpz_init(fib_next);

  // Calculate the Fibonacci numbers iteratively
  for (int i = 2; i <= n; i++) {
    mpz_add(fib_next, fib_prev, fib_curr); // fib_next = fib_prev + fib_curr
    mpz_set(fib_prev, fib_curr);       // fib_prev = fib_curr
    mpz_set(fib_curr, fib_next);       // fib_curr = fib_next
  }

  // Set the result
  mpz_set(result, fib_curr);

  // Clear the mpz_t variables
  mpz_clear(fib_prev);
  mpz_clear(fib_curr);
  mpz_clear(fib_next);
}

void fibonacci_iterative(mpz_t result, int n) {
  if (n <= 1) {
    mpz_set_ui(result, n);
    return;
  }

  mpz_t a, b, c;
  mpz_init_set_ui(a, 0);  // a = 0
  mpz_init_set_ui(b, 1);  // b = 1
  mpz_init(c);

  for (int i = 2; i <= n; i++) {
    mpz_add(c, a, b);  // c = a + b

    mpz_set(a, b);   // a = b
    mpz_set(b, c);   // b = c
  }

  mpz_set(result, b);  // result = b

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(c);
}

int find_n_for_target_time_binary(int low, int high, int target_time_sec, double *current_time_used, double *total_time_used, mpz_t result) {
  int mid;
  int bestN = low;
  double bestTime = 0.0;
  *current_time_used = 0.0;

  clock_t overall_start = clock();

  while (low <= high) {
    mid = low + (high - low) / 2;

    clock_t start = clock();
    fibonacci_recursive(result, mid);
    clock_t end = clock();

    double current_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    if (current_time_used <= target_time_sec) {
      if (current_time_used > bestTime) {
        bestN = mid;
        bestTime = current_time_used;
      }
      low = mid + 1;
    }
    else {
      high = mid - 1;
    }
  }

  clock_t overall_end = clock();
  *total_time_used = ((double) (overall_end - overall_start)) / CLOCKS_PER_SEC;
  *current_time_used = bestTime;

  return bestN;
}

void benchmark(int binary_low, int binary_high, int target_time_sec, int n_value, int benchmark_flag, int compare_flag, int slow_flag) {
  int n;

  double time_taken_recursive = 0.0, time_taken_dynamic = 0.0, time_taken_iterative = 0.0, time_taken_matrix = 0.0;
  double current_time_used = 0.0, time_taken_to_find_n = 0.0;

  mpz_t result_recursive, result_dynamic, result_iterative, result_matrix;
  mpz_init(result_recursive);
  mpz_init(result_dynamic);
  mpz_init(result_iterative);
  mpz_init(result_matrix);

  if (benchmark_flag) {
    n = find_n_for_target_time_binary(binary_low, binary_high, target_time_sec, &current_time_used, &time_taken_to_find_n, result_recursive);
    printf("Time taken to find n (range: low = %d, high = %d) using binary search: ", binary_low, binary_high);
    printTime(time_taken_to_find_n);

    gmp_printf("Recursive Approach:\t n = %d\nResult = %Zd\nTime taken: ", n, result_recursive);
    printTime(current_time_used);
  }
  else {
    n = n_value;
  }

  time_taken_matrix = measureTime(fibonacci_matrix, n, result_matrix);
  gmp_printf("Matrix Approach:\n\tn = %d\n\tResult = %Zd\n\tTime taken: ", n, result_matrix);
  printTime(time_taken_matrix);

  if (compare_flag || benchmark_flag) {
    time_taken_iterative = measureTime(fibonacci_iterative, n, result_iterative);
    gmp_printf("Iterative Approach:\n\tn = %d\n\tResult = %Zd\n\tTime taken: ", n, result_iterative);
    printTime(time_taken_iterative);

    if (slow_flag) {
      time_taken_recursive = measureTime(fibonacci_recursive, n, result_recursive);
      gmp_printf("Recursive Approach:\n\tn = %d\n\tResult = %Zd\n\tTime taken: ", n, result_recursive);
      printTime(time_taken_recursive);
    }

    time_taken_dynamic = measureTime(fibonacci_dynamic, n, result_dynamic);
    gmp_printf("Dynamic Approach:\n\tn = %d\n\tResult = %Zd\n\tTime taken: ", n, result_dynamic);
    printTime(time_taken_dynamic);

    if (mpz_cmp(result_matrix, result_iterative) != 0 || mpz_cmp(result_dynamic, result_iterative) != 0) {
      printf("Warning: Inconsistent results between methods.\n");
    }
    else {
      printf("The results match.\n");
    }
  }

  mpz_clear(result_recursive);
  mpz_clear(result_dynamic);
  mpz_clear(result_iterative);
  mpz_clear(result_matrix);
}

int main(int argc, char *argv[]) {
  int binary_low = 1000, binary_high = 5000, target_time_sec = 30;
  int benchmark_flag = 0, compare_flag = 0, slow_flag = 0;
  unsigned long long n_value = 10000;

  processArgs(argc, argv, &binary_low, &binary_high, &target_time_sec, &benchmark_flag, &compare_flag, &slow_flag, &n_value);


  benchmark(binary_low, binary_high, target_time_sec, n_value, benchmark_flag, compare_flag, slow_flag);

  int n_max = calcMaxFib();
  printf("\nThe largest problem size n_max for Fibonacci numbers using long int without overflow is: %d\n", n_max);

  return 0;
}

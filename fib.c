#include <errno.h>
#include <float.h>
#include <gmp.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#define MAX_N 1000000

#include <stdio.h>

void printUsage(const char *programName) {
  (void) fputs("Usage: ", stdout);
  (void) fputs(programName, stdout);
  (void) fputs(" [OPTIONS]\n", stdout);
  (void) fputs("Example: ", stdout);
  (void) fputs(programName, stdout);
  (void) fputs(" -B -L 40 -H 55 -T 35\n", stdout);
  (void) fputs("Example: ", stdout);
  (void) fputs(programName, stdout);
  (void) fputs(" -C -N 50\n", stdout);
  (void) fputs("Example: ", stdout);
  (void) fputs(programName, stdout);
  (void) fputs(" -R\n\n", stdout);

  (void) fputs("Options:\n", stdout);
  (void) fputs("  -B, --benchmark   Run benchmarking. If set, allows the following options:\n", stdout);
  (void) fputs("  -L, --low <value>   Set the low bound for binary search. (Requires -B)\n", stdout);
  (void) fputs("  -H, --high <value>  Set the high bound for binary search. (Requires -B)\n", stdout);
  (void) fputs("  -T, --time <value>  Set the target time (in seconds) for the search. (Requires -B)\n", stdout);
  (void) fputs("  -C, --compare   Run comparison of all the methods for a specified 'n'.\n", stdout);
  (void) fputs("  -S, --slow  Include recursive method in comparison. (Use with -C)\n", stdout);
  (void) fputs("  -N, --nvalue <value> Specify the 'n' value for Fibonacci calculation. (Required for -C)\n", stdout);
  (void) fputs("  -R, --results  Print the results of the Fibonacci calculations.\n", stdout);
  (void) fputs("  -h, --help  Show this help message and exit.\n\n", stdout);

  (void) fputs("Description:\n", stdout);
  (void) fputs("  This program calculates the Fibonacci number for a given value using different methods.\n", stdout);
  (void) fputs("  The benchmark (-B) option runs a binary search to find an 'n' value for which the recursive\n", stdout);
  (void) fputs("  Fibonacci calculation takes approximately the specified time. The low, high, and time values\n", stdout);
  (void) fputs("  are used to control this search. The compare (-C) option compares the results and execution\n", stdout);
  (void) fputs("  times of dynamic and iterative methods for a specified 'n' value, and includes recursive\n", stdout);
  (void) fputs("  method if -S is specified. The results option (-R) controls whether calculation results\n", stdout);
  (void) fputs("  are printed. Without -B, -C, or -R, the program simply performs the iterative approach.\n", stdout);
}

void printDigitCount(mpz_t number) {
  int digitCount = 0;

  if (mpz_sgn(number) == 0) { // Check if the number is zero
    digitCount = 1;
  } else {
    digitCount = mpz_sizeinbase(number, 10); // Count the digits in base 10
  }

  printf("The result is %d digits long. (Use -R to print to screen)\n\n", digitCount);
}

void clearFibMemo() {
  static mpz_t memo[MAX_N];
  static bool initialized[MAX_N];

  for (int i = 0; i < MAX_N; ++i) {
    if (initialized[i]) {
      mpz_clear(memo[i]);
    }
  }
}

void fibonacci_aux(mpz_t result, int n, mpz_t memo[], bool initialized[]) {
  for (int i = 0; i <= n; ++i) {
    if (!initialized[i]) {
      mpz_add(memo[i], memo[i - 1], memo[i - 2]);
      initialized[i] = true;
    }
  }
  mpz_set(result, memo[n]);
}

void fibonacci_recursive(mpz_t result, int n) {
  static mpz_t memo[MAX_N];
  static bool initialized[MAX_N] = {false};

  static bool is_first_call = true;
  if (is_first_call) {
    for (int i = 0; i < MAX_N; ++i) {
      mpz_init(memo[i]);
      initialized[i] = false;
    }
    is_first_call = false;
  }

  fibonacci_aux(result, n, memo, initialized);
}

unsigned long long int fibRecursiveBench(unsigned long long int n) {
  if (n <= 1)
    return n;

  return fibRecursiveBench(n-1) + fibRecursiveBench(n-2);
}

void fibIterativeBench(int n) {
  if (n < 0) {
    printf("Please enter a non-negative integer.\n");
    return;
  }

  long long a = 0, b = 1;
  for (int i = 2; i <= n; ++i) {
    if (a > LONG_MAX - b) {
      printf("Long Int Overflow occurred at Fibonacci position: %d\n\n", i);
      return;
    }

    long long next = a + b;
    a = b;
    b = next;
  }

  if (n <= 1) {
    printf("Fibonacci number at position %d: %lld\n", n, (long long)n);
  }
  else {
    printf("Fibonacci number at position %d: %lld\n", n, b);
  }
}

int calcMaxFib() {
    const double phi = (1 + sqrt(5)) / 2;
    int n = (int)floor(log((double)LONG_MAX * sqrt(5)) / log(phi));
    return n;
}

int calcMaxFibIterative() {
  long a = 0, b = 1;
  int n = 1;

  while (1) {
    long c = a + b;
    if (c < a || c < b) { // Overflow detect
      break;
    }
    a = b;
    b = c;
    n++;
  }

  return n;
}

void compareFib() {
  int maxFibIndex1 = calcMaxFib();
  int maxFibIndex2 = calcMaxFibIterative();

  printf("Long Int Max Fibonacci Index (Method 1): %d\n", maxFibIndex1);
  printf("Long Int Max Fibonacci Index (Method 2): %d\n", maxFibIndex2);

  if (maxFibIndex1 == maxFibIndex2) {
    printf("The results match.\n");
  }
  else {
    printf("Warning: Inconsistent results between methods.\n");
  }
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

// Power of a 2x2 matrix.
void matrixPower(Matrix2x2 *result, Matrix2x2 *m, int n) {
  Matrix2x2 temp; // Temp matrix
  matrixInit(&temp); // Initialize
  matrixSetIdentity(&temp); // Set to the identity matrix.

  // if n is 0, the result is the identity matrix.
  if (n == 0) {
    matrixCopy(result, &temp); // Copy the identity matrix to the result.
    return;
  }
  // if n is 1, the result is the matrix itself.
  else if (n == 1) {
    matrixCopy(result, m); // Copy input matrix to the result.
    return;
  }

  // Compute the matrix raised to n/2.
  matrixPower(result, m, n / 2);
  // Multiply result of recursion by itself
  matrixMultiply(&temp, result, result);

  // If n is odd, multiply the result by the matrix again.
  if (n % 2 != 0) {
    matrixMultiply(result, &temp, m);
  }
  // If n is even, just copy the squared matrix to the result.
  else {
    matrixCopy(result, &temp);
  }

  // Clear the mpz_t types
  mpz_clear(temp.a);
  mpz_clear(temp.b);
  mpz_clear(temp.c);
  mpz_clear(temp.d);
}


// nth Fibonacci number using matrix exponentiation.
void fibonacci_matrix(mpz_t result, int n) {
  // if n is 0, set the result to 0 and return.
  if (n == 0) {
    mpz_set_ui(result, 0); // Set the result to 0
    return;
  }

  Matrix2x2 base, res; // Declare/init base and the result.
  matrixInit(&base);
  matrixInit(&res);

  // Setup the base matrix
  mpz_set_ui(base.a, 0); // Set base[0][0] to 0.
  mpz_set_ui(base.b, 1); // Set base[0][1] to 1.
  mpz_set_ui(base.c, 1); // Set base[1][0] to 1.
  mpz_set_ui(base.d, 1); // Set base[1][1] to 1.

  // Power of the base matrix to the nth term.
  matrixPower(&res, &base, n);

  // n is located at res[0][1] (res.b).
  mpz_set(result, res.b); // Copy res.b to the result.

  // Clear all
  mpz_clear(base.a);
  mpz_clear(base.b);
  mpz_clear(base.c);
  mpz_clear(base.d);
  mpz_clear(res.a);
  mpz_clear(res.b);
  mpz_clear(res.c);
  mpz_clear(res.d);
}


unsigned long long int parseIntArg(char *arg) {
  errno = 0;
  char *endptr;
  long val = strtol(arg, &endptr, 10);

  if (errno != 0 || *endptr != '\0' || val < 0 || val > MAX_N) {
    fprintf(stderr, "Please enter a number between 0 and 1000000\n");
    exit(1);
  }

  return (int)val;
}

void processArgs(int argc, char *argv[], int *binary_low, int *binary_high,
         int *target_time_sec, int *benchmark_flag, int *compare_flag,
         int *slow_flag, int *results_flag, unsigned long long int *n_value) {
  *benchmark_flag = 0;
  *compare_flag = 0;
  *slow_flag = 0;
  *results_flag = 0;
  int high_flag = 0;
  int low_flag = 0;
  int time_flag = 0;
  int n_flag = 0;

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-B") == 0 || strcmp(argv[i], "--benchmark") == 0) {
      *benchmark_flag = 1;
    }
    else if (strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--compare") == 0) {
      *compare_flag = 1;
    }
    else if (strcmp(argv[i], "-S") == 0 || strcmp(argv[i], "--slow") == 0) {
      *slow_flag = 1;
    }
    else if ((strcmp(argv[i], "-L") == 0 || strcmp(argv[i], "--low") == 0) && i + 1 < argc) {
      *binary_low = parseIntArg(argv[++i]);
      low_flag = 1;
    }
    else if ((strcmp(argv[i], "-H") == 0 || strcmp(argv[i], "--high") == 0) && i + 1 < argc) {
      *binary_high = parseIntArg(argv[++i]);
      high_flag = 1;
    }
    else if ((strcmp(argv[i], "-T") == 0 || strcmp(argv[i], "--time") == 0) && i + 1 < argc) {
      *target_time_sec = parseIntArg(argv[++i]);
      time_flag = 1;
    }
    else if ((strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--nvalue") == 0) && i + 1 < argc) {
      *n_value = strtoull(argv[++i], NULL, 10);
      n_flag = 1;
    }
    else if (strcmp(argv[i], "-R") == 0 || strcmp(argv[i], "--results") == 0) {
      *results_flag = 1;
    }
    else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      printUsage(argv[0]);
      exit(0);
    }
    else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      printUsage(argv[0]);
      exit(1);
    }
  }

  // Rule enforcement
  if (*benchmark_flag && *compare_flag) {
    fprintf(stderr, "Error: -B (--benchmark) and -C (--compare) cannot be used together.\n");
    exit(1);
  }
  if ((low_flag || high_flag || time_flag) && !*benchmark_flag) {
    fprintf(stderr, "Error: -L (--low), -H (--high), and -T (--time) require -B (--benchmark).\n");
    exit(1);
  }
  if (*slow_flag && !*compare_flag) {
    fprintf(stderr, "Error: -S (--slow) requires -C (--compare).\n");
    exit(1);
  }
  if (n_flag && *benchmark_flag) {
    fprintf(stderr, "Error: -N (--nvalue) cannot be used with -B (--benchmark).\n");
    exit(1);
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
  if (n <= 1) {
    mpz_set_ui(result, n);
    return;
  }

  mpz_t *fib_cache = malloc((n + 1) * sizeof(mpz_t));
  mpz_init_set_ui(fib_cache[0], 0); // fib[0] = 0
  mpz_init_set_ui(fib_cache[1], 1); // fib[1] = 1

  for (int i = 2; i <= n; ++i) {
    mpz_init(fib_cache[i]);
    mpz_add(fib_cache[i], fib_cache[i - 1], fib_cache[i - 2]);
  }

  mpz_set(result, fib_cache[n]);

  for (int i = 0; i <= n; ++i) {
    mpz_clear(fib_cache[i]);
  }
  free(fib_cache);
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

int find_n_for_target_time_binary(int low, int high, int target_time_sec,
                                  double *cpu_time_used, double *total_time_used,
                                  unsigned long long int *benchResult) {
  int mid;
  int bestN = low;
  double bestTime = 0.0;
  *cpu_time_used = 0.0;

  clock_t overall_start = clock();

  printf("Finding 'n' value using binary search. Please wait...\n");

  while (low <= high) {
    mid = low + (high - low) / 2;

    clock_t start = clock();
    *benchResult = fibRecursiveBench(mid);
    clock_t end = clock();

    double current_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("\tLow: %d, High: %d, Mid: %d, Time: %f\n", low, high, mid, current_time_used);

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
  *cpu_time_used = bestTime;

  printf("'n' found. Using %d for further runs.\n\n", bestN);

  return bestN;
}

void benchmark(int binary_low, int binary_high, int target_time_sec, int n_value,
               int benchmark_flag, int compare_flag, int slow_flag, int results_flag) {
  int n;

  double time_taken_recursive = 0.0, time_taken_dynamic = 0.0, time_taken_iterative = 0.0, time_taken_matrix = 0.0;
  double current_time_used = 0.0, time_taken_to_find_n = 0.0;

  unsigned long long int benchResult = 0;

  mpz_t result_recursive, result_dynamic, result_iterative, result_matrix;
  mpz_init(result_recursive);
  mpz_init(result_dynamic);
  mpz_init(result_iterative);
  mpz_init(result_matrix);

  if (benchmark_flag) {
    n = find_n_for_target_time_binary(binary_low, binary_high, target_time_sec, &current_time_used, &time_taken_to_find_n, &benchResult);
    printf("Time taken to find n (range: low = %d, high = %d) using binary search: ", binary_low, binary_high);
    printTime(time_taken_to_find_n);
    printf("\n");

    if (results_flag) {
      printf("Recursive Approach (n = %d, Result = %llu): ", n, benchResult);
    }
    printf("(n = %d)\n", n);
    printf("Recursive Time:\t");
    printTime(current_time_used);
  }
  else {
    n = n_value;
    printf("(n = %d)\n", n);
  }

  // Matrix Approach
  time_taken_matrix = measureTime(fibonacci_matrix, n, result_matrix);
  if (results_flag) {
    gmp_printf("Matrix Approach (n = %d, Result = %Zd): ", n, result_matrix);
  }
  printf("Matrix Time:\t");
  printTime(time_taken_matrix);
  printDigitCount(result_matrix); // Call the function here

  // Iterative Approach
  if (compare_flag || benchmark_flag) {
    time_taken_iterative = measureTime(fibonacci_iterative, n, result_iterative);
    if (results_flag) {
      gmp_printf("Iterative Approach (n = %d, Result = %Zd): ", n, result_iterative);
    }
    printf("Iterative Time:\t");
    printTime(time_taken_iterative);

    // Recursive Approach
    if (slow_flag) {
      time_taken_recursive = measureTime(fibonacci_recursive, n, result_recursive);
      if (results_flag) {
        gmp_printf("Recursive Approach (n = %d, Result = %Zd): ", n, result_recursive);
      }
      printf("Recursive Time:\t");
      printTime(time_taken_recursive);
    }

    // Dynamic Approach
    time_taken_dynamic = measureTime(fibonacci_dynamic, n, result_dynamic);
    if (results_flag) {
      gmp_printf("Dynamic Approach (n = %d, Result = %Zd): ", n, result_dynamic);
    }
    printf("Dynamic Time:\t");
    printTime(time_taken_dynamic);

    // Comparing results
    if (mpz_cmp(result_matrix, result_iterative) != 0 || mpz_cmp(result_dynamic, result_iterative) != 0) {
      printf("Warning: Inconsistent results between methods.\n\n");
    }
    else {
      printf("The results match.\n\n");
    }
  }

  mpz_clear(result_recursive);
  mpz_clear(result_dynamic);
  mpz_clear(result_iterative);
  mpz_clear(result_matrix);
}

int main(int argc, char *argv[]) {
  int binary_low = 30, binary_high = 55, target_time_sec = 30;
  int benchmark_flag = 0, compare_flag = 0, slow_flag = 0, results_flag = 0;
  unsigned long long n_value = 1000000;

  processArgs(argc, argv, &binary_low, &binary_high, &target_time_sec, &benchmark_flag, &compare_flag, &slow_flag, &results_flag, &n_value);

  benchmark(binary_low, binary_high, target_time_sec, n_value, benchmark_flag, compare_flag, slow_flag, results_flag);

  fibIterativeBench(93);

  compareFib();

  clearFibMemo();

  return 0;
}

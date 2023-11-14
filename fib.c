#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>

void printUsage(const char *programName) {
  printf("Usage: %s [OPTIONS]\n", programName);
  printf("Example: %s -B -L 40 -H 55 -T 35\n", programName);
  printf("Example: %s -C -N 50\n\n", programName);

  printf("Options:\n");
  printf("  -B, --benchmark   Run benchmarking. If set, allows the following options:\n");
  printf("  -L, --low <value>   Set the low bound for binary search. (Requires -B)\n");
  printf("  -H, --high <value>  Set the high bound for binary search. (Requires -B)\n");
  printf("  -T, --time <value>  Set the target time (in seconds) for the search. (Requires -B)\n");
  printf("  -C, --compare     Run comparison of all three Fibonacci methods for a specified 'n'.\n");
  printf("  -N, --nvalue <value> Specify the 'n' value for Fibonacci calculation. (Required for -C)\n");
  printf("  -h, --help      Show this help message and exit.\n\n");

  printf("Description:\n");
  printf("  This program calculates the Fibonacci number for a given value using different methods.\n");
  printf("  The benchmark (-B) option runs a binary search to find an 'n' value for which the recursive\n");
  printf("  Fibonacci calculation takes approximately the specified time. The low, high, and time values\n");
  printf("  are used to control this search. The compare (-C) option compares the results and execution\n");
  printf("  times of recursive, dynamic, and iterative methods for a specified 'n' value. Without -B or -C,\n");
  printf("  the program simply performs the iterative approach.\n");
}

void processArgs(int argc, char *argv[], int *binary_low, int *binary_high, int *target_time_sec, int *benchmark_flag, int *compare_flag, int *n_value) {
  *benchmark_flag = 1;
  *compare_flag = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-B") == 0 || strcmp(argv[i], "--benchmark") == 0) {
      *benchmark_flag = 1;
    }
    else if (strcmp(argv[i], "-C") == 0 || strcmp(argv[i], "--compare") == 0) {
      *compare_flag = 1;
    }
    else if ((strcmp(argv[i], "-L") == 0 || strcmp(argv[i], "--low") == 0) && i + 1 < argc) {
      if (*benchmark_flag == 0) {
        printf("Error: -L/--low requires -B/--benchmark.\n");
        exit(1);
      }
      *binary_low = atoi(argv[++i]);
    }
    else if ((strcmp(argv[i], "-H") == 0 || strcmp(argv[i], "--high") == 0) && i + 1 < argc) {
      if (*benchmark_flag == 0) {
        printf("Error: -H/--high requires -B/--benchmark.\n");
        exit(1);
      }
      *binary_high = atoi(argv[++i]);
    }
    else if ((strcmp(argv[i], "-T") == 0 || strcmp(argv[i], "--time") == 0) && i + 1 < argc) {
      if (*benchmark_flag == 0) {
        printf("Error: -T/--time requires -B/--benchmark.\n");
        exit(1);
      }
      *target_time_sec = atoi(argv[++i]);
    }
    else if ((strcmp(argv[i], "-N") == 0 || strcmp(argv[i], "--nvalue") == 0) && i + 1 < argc) {
      *n_value = atoi(argv[++i]);
    }
    else {
      printf("Unknown option: %s\n", argv[i]);
      printUsage(argv[0]);
      exit(1);
    }
  }
}

double measureTime(long int (*func)(int, int, int), int arg1, int arg2, int arg3) {
  struct timespec start, end;

  // CLOCK_MONOTONIC provides a monotonic time since some arbitrary start point
  clock_gettime(CLOCK_MONOTONIC, &start);
  func(arg1, arg2, arg3);
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
    // milliseconds
    printf("%f milliseconds\n", timeInSeconds * 1e3);
  }
  // 1e-6 represents 0.000001 seconds (1 microsecond)
  else if (timeInSeconds >= 1e-6) {
    // microseconds
    printf("%f microseconds\n", timeInSeconds * 1e6);
  }
  // Assume nanoseconds at this point (1 sec = 1,000,000,000 ns)
  else {
    printf("%f nanoseconds\n", timeInSeconds * 1e9);
  }
}

long int fibonacci_recursive(int n) {
  // If n is 0 or 1, first two numbers in the Fibonacci sequence are 0 and 1.
  if (n <= 1)
    return n;

  // Calls itself twice. Once with (n-1) and once with (n-2), and then adds
  // the results of these two calls. This is based on the principle that the
  // nth Fibonacci number is the sum of the (n-1)th and (n-2)th Fibonacci numbers.
  return fibonacci_recursive(n-1) + fibonacci_recursive(n-2);
}

long int fibonacci_dynamic(int n) {
  // If n is 0 or 1, first two numbers in the Fibonacci sequence are 0 and 1.
  if (n <= 1) {
    return n;
  }

  // Create an array to store Fibonacci numbers up to n. (n+1) to accommodate
  // all values from 0 to n.
  long int fib[n+1];

  // fib[0] is the first Fibonacci number, which is 0.
  // fib[1] is the second Fibonacci number, which is 1.
  fib[0] = 0;
  fib[1] = 1;

  // Start at 2 because the first two Fibonacci numbers are already known.
  for (int i = 2; i <= n; i++) {
    // Calculate the ith Fibonacci number.
    // It is the sum of the two preceding numbers, fib[i-1] and fib[i-2].
    fib[i] = fib[i-1] + fib[i-2];
  }

  // nth Fibonacci number
  return fib[n];
}


long int fibonacci_iterative(int n) {
  // If n is 0 or 1, first two numbers in the Fibonacci sequence are 0 and 1.
  if (n <= 1) {
    return n;
  }

  // 'a' is initialized to 0 (the 0th Fibonacci number).
  // 'b' is initialized to 1 (the 1st Fibonacci number).
  long int a = 0, b = 1;

  // 'c' is tmp
  long int c;

  // Start at 2 because the first two Fibonacci numbers are already known.
  for (int i = 2; i <= n; i++) {
    // Next one
    c = a + b;

    // Shift values
    a = b;
    b = c;
  }

  // After the shift, 'b' holds the nth Fibonacci number.
  // Return 'b', which is the requested Fibonacci number.
  return b;
}


int find_n_for_target_time_binary(int low, int high, int target_time_sec, double *cpu_time_used, double *total_time_used, long int *result) {
  int mid;
  int bestN = low;
  double bestTime = 0.0;
  *cpu_time_used = 0.0;

  clock_t overall_start = clock();

  while (low <= high) {
    mid = low + (high - low) / 2;

    clock_t start = clock();
    *result = fibonacci_recursive(mid);
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
  *cpu_time_used = bestTime;

  return bestN;
}

long int fibonacci_recursive_wrapper(int n, int unused1, int unused2) {
  return fibonacci_recursive(n);
}

long int fibonacci_dynamic_wrapper(int n, int unused1, int unused2) {
  return fibonacci_dynamic(n);
}

long int fibonacci_iterative_wrapper(int n, int unused1, int unused2) {
  return fibonacci_iterative(n);
}

void benchmark(int binary_low, int binary_high, int target_time_sec, int n_value, int benchmark_flag, int compare_flag) {
  int n;
  long int result_recursive = 0, result_dynamic = 0, result_iterative = 0;
  double time_taken_recursive = 0.0, time_taken_dynamic = 0.0, time_taken_iterative = 0.0;

  if (benchmark_flag && !compare_flag) {
    double time_taken_to_find_n;
    n = find_n_for_target_time_binary(binary_low, binary_high, target_time_sec,
                      &time_taken_recursive, &time_taken_to_find_n, &result_recursive);
    printf("Time taken to find n (range: low = %d, high = %d) using binary search: ", binary_low, binary_high);
    printTime(time_taken_to_find_n);
  }
  else {
    n = n_value;
  }

  result_iterative = fibonacci_iterative(n);
  time_taken_iterative = measureTime(fibonacci_iterative_wrapper, n, 0, 0);
  printf("Iterative Approach for n = %d: Result = %ld, Time taken: ", n, result_iterative);
  printTime(time_taken_iterative);

  if (compare_flag && !benchmark_flag) {
    result_recursive = fibonacci_dynamic(n);
    time_taken_recursive = measureTime(fibonacci_recursive_wrapper, n, 0, 0);
  }

  if (benchmark_flag || compare_flag) {
    printf("Recursive Approach for n = %d: Result = %ld, Time taken: ", n, result_recursive);
    printTime(time_taken_recursive);

    result_dynamic = fibonacci_dynamic(n);
    time_taken_dynamic = measureTime(fibonacci_dynamic_wrapper, n, 0, 0);
    printf("Dynamic Approach for n = %d: Result = %ld, Time taken: ", n, result_dynamic);
    printTime(time_taken_dynamic);

    if (result_recursive != result_iterative || result_dynamic != result_iterative) {
      printf("Warning: Inconsistent results between methods.\n");
    }
    else {
      printf("The results match.\n");
    }
  }
}

int calcMaxFib() {
  const double phi = (1 + sqrt(5)) / 2; // golden ratio
  return (int) floor(log(LONG_MAX * sqrt(5)) / log(phi));
}

int main(int argc, char *argv[]) {
  int binary_low = 43, binary_high = 53, target_time_sec = 30, n_value = 51;
  int benchmark_flag = 0, compare_flag = 0;

  processArgs(argc, argv, &binary_low, &binary_high, &target_time_sec, &benchmark_flag, &compare_flag, &n_value);

  benchmark(binary_low, binary_high, target_time_sec, n_value, benchmark_flag, compare_flag);

  int n_max = calcMaxFib();
  printf("The largest problem size n_max for Fibonacci numbers using long int without overflow is: %d\n", n_max);

  return 0;
}

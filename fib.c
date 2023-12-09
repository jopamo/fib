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

#define MAX_N 10000000

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
  (void) fputs("\t  -L, --low <value>   Set the low bound for binary search. (Requires -B)\n", stdout);
  (void) fputs("\t  -H, --high <value>  Set the high bound for binary search. (Requires -B)\n", stdout);
  (void) fputs("\t  -T, --time <value>  Set the target time (in seconds) for the search. (Requires -B)\n", stdout);
  (void) fputs("  -C, --compare   Run comparison of all the methods for a specified 'n'.\n", stdout);
  (void) fputs("  -N, --nvalue <value> Specify the 'n' value for Fibonacci calculation. (Required for -C)\n", stdout);
  (void) fputs("  -S, --slow  Include recursive method in comparison. (Use with -C)\n", stdout);
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
  }
  else {
    digitCount = mpz_sizeinbase(number, 10); // Count the digits in base 10
  }

  printf("The result is %d digits long. (Use -R to print to screen)\n\n", digitCount);
}

/**
 * Clears the memoization array used in the Fibonacci calculation.
 * This function iterates over the static arrays used for memoization,
 * clearing the memory allocated for each mpz_t type element that was initialized.
 */
void clearFibMemo() {
    static mpz_t memo[MAX_N];  // Array for storing Fibonacci numbers
    static bool initialized[MAX_N];  // Array to track initialization status

    for (int i = 0; i < MAX_N; ++i) {
        if (initialized[i]) {
            mpz_clear(memo[i]);  // Clear memory if element was initialized
        }
    }
}


/**
 * Auxiliary function to compute Fibonacci numbers using memoization.
 * This function iteratively calculates Fibonacci numbers up to the nth index,
 * storing each computed number in a memoization array to avoid redundant calculations.
 *
 * @param result The mpz_t variable to store the nth Fibonacci number.
 * @param n The index of the Fibonacci number to compute.
 * @param memo Array of mpz_t types for storing computed Fibonacci numbers.
 * @param initialized Array indicating if a particular index in the memo array has been computed.
 */
void fibonacci_aux(mpz_t result, int n, mpz_t memo[], bool initialized[]) {
    for (int i = 0; i <= n; ++i) {
        if (!initialized[i]) {
            mpz_add(memo[i], memo[i - 1], memo[i - 2]);  // Calculate the Fibonacci number
            initialized[i] = true;  // Mark as computed
        }
    }

    mpz_set(result, memo[n]);  // Set the result to the nth Fibonacci number
}



/**
 * Calculates Fibonacci numbers using a recursive approach with memoization.
 * This function initializes static arrays for memoization on the first call and
 * uses an auxiliary function for the actual computation of the Fibonacci numbers.
 * The memoization technique avoids redundant calculations, enhancing efficiency.
 *
 * @param result The mpz_t variable to store the nth Fibonacci number.
 * @param n The index of the Fibonacci number to compute.
 */
void fibonacci_recursive(mpz_t result, int n) {
    static mpz_t memo[MAX_N];  // Array for storing Fibonacci numbers
    static bool initialized[MAX_N] = {false};  // Array to track initialization status
    static bool is_first_call = true;  // Flag to check the first function call

    // Initialize memoization array on the first call
    if (is_first_call) {
        for (int i = 0; i < MAX_N; ++i) {
            mpz_init(memo[i]);  // Initialize each mpz_t in the array
            initialized[i] = false;  // Mark as not initialized
        }
        is_first_call = false;  // Update the flag after first initialization
    }

    // Compute the Fibonacci number using the auxiliary function
    fibonacci_aux(result, n, memo, initialized);
}



/**
 * Calculates Fibonacci numbers using a recursive approach.
 * It follows the classic recursive definition of the Fibonacci
 * sequence but is inefficient for large n due to redundant calculations.
 * @param n The index of the Fibonacci number to calculate.
 * @return The Fibonacci number at the given index.
 */
unsigned long long int fibRecursiveBench(unsigned long long int n) {
    if (n <= 1) return n;  // Base case: return n for 0 or 1

    // Recursive case: sum of the two preceding numbers
    return fibRecursiveBench(n-1) + fibRecursiveBench(n-2);
}



/**
 * Calculates Fibonacci numbers using an iterative approach and tests for overflow.
 * It is designed to check the point at which computing a Fibonacci number
 * results in an overflow for a long long integer.
 *
 * @param n The index of the Fibonacci number to calculate.
 */
void fibIterativeBench(int n) {
    if (n < 0) {
        printf("Please enter a non-negative integer.\n");  // Validate input
        return;
    }

    long long a = 0, b = 1;  // Initialize first two Fibonacci numbers

    for (int i = 2; i <= n; ++i) {
        if (a > LONG_MAX - b) {  // Check for overflow
            printf("Long Int Overflow occurred at Fibonacci position: %d\n\n", i);
            return;  // Exit if overflow occurs
        }

        long long next = a + b;  // Compute next Fibonacci number
        a = b;  // Update previous number
        b = next;  // Update current number
    }

    // Output Fibonacci number for n
    printf("Fibonacci number at position %d: %lld\n", n, (n <= 1) ? (long long)n : b);
}



/**
 * Calculates the maximum Fibonacci index using the Golden Ratio.
 * It uses the properties of the Fibonacci sequence and the Golden Ratio to determine
 * the largest index of a Fibonacci number that can fit within a long integer range.
 *
 * @return The maximum index of a Fibonacci number within the long integer range.
 */
int calcMaxFib() {
    const double phi = (1 + sqrt(5)) / 2;  // Golden Ratio

    // Calculate and return the maximum Fibonacci index
    return (int)floor(log((double)LONG_MAX * sqrt(5)) / log(phi));
}


/**
 * Calculates the maximum Fibonacci index iteratively.
 * It iteratively computes Fibonacci numbers until an overflow occurs,
 * thereby finding the largest index that can be represented without overflow.
 *
 * @return The maximum index of a Fibonacci number that can be calculated without overflow.
 */
int calcMaxFibIterative() {
    long a = 0, b = 1;  // Initialize the first two Fibonacci numbers
    int n = 1;  // Start from index 1

    while (1) {
        long c = a + b;  // Compute the next Fibonacci number

        // Check for overflow
        if (c < a || c < b) break;  // Detect overflow and break the loop

        // Update Fibonacci numbers for next iteration
        a = b;
        b = c;
        n++;  // Increment index
    }

    return n;  // Return the maximum index found
}


/**
 * Compares the maximum Fibonacci index calculable using two different methods.
 * It calculates the maximum Fibonacci index using a formula-based approach and an iterative approach,
 * then compares and prints the results from both methods.
 */
void compareFib() {
    int maxFibIndex1 = calcMaxFib();  // Calculate using formula-based approach
    int maxFibIndex2 = calcMaxFibIterative();  // Calculate using iterative approach

    // Print results from both methods
    printf("Long Int Max Fibonacci Index (Method 1): %d\n", maxFibIndex1);
    printf("Long Int Max Fibonacci Index (Method 2): %d\n", maxFibIndex2);

    // Compare and print appropriate message
    if (maxFibIndex1 == maxFibIndex2) {
        printf("The results match.\n");  // Indicate matching results
    } else {
        printf("Warning: Inconsistent results between methods.\n");  // Warn about discrepancy
    }
}


// A structure to represent a 2x2 matrix using GMP's mpz_t type for each element
typedef struct {
  mpz_t a, b, c, d; // Elements of the matrix: a, b, c, d
} Matrix2x2;


// Initialize each element of the matrix with GMP's mpz_init function
// This is necessary for mpz_t types before they can be used
void matrixInit(Matrix2x2 *m) {
  mpz_init(m->a);
  mpz_init(m->b);
  mpz_init(m->c);
  mpz_init(m->d);
}


/**
 * Sets a 2x2 matrix to the identity matrix.
 * In an identity matrix, diagonal elements are 1 and off-diagonal elements are 0.
 *
 * @param m The matrix to be set to identity.
 */
void matrixSetIdentity(Matrix2x2 *m) {
    mpz_set_ui(m->a, 1);  // Set top-left element (diagonal) to 1
    mpz_set_ui(m->d, 1);  // Set bottom-right element (diagonal) to 1

    mpz_set_ui(m->b, 0);  // Set top-right element (off-diagonal) to 0
    mpz_set_ui(m->c, 0);  // Set bottom-left element (off-diagonal) to 0
}


/**
 * Copies the contents of one 2x2 matrix to another.
 *
 * @param dest The destination matrix.
 * @param src The source matrix.
 */
void matrixCopy(Matrix2x2 *dest, Matrix2x2 *src) {
    mpz_set(dest->a, src->a);  // Copy top-left element
    mpz_set(dest->b, src->b);  // Copy top-right element
    mpz_set(dest->c, src->c);  // Copy bottom-left element
    mpz_set(dest->d, src->d);  // Copy bottom-right element
}


// multiply two 2x2 matrices
void matrixMultiply(Matrix2x2 *result, Matrix2x2 *m1, Matrix2x2 *m2) {
  // Temporary matrix to hold the intermediate results of multiplication
  Matrix2x2 temp;
  matrixInit(&temp); // Initialize the temporary matrix

  // Perform multiplication of m1 and m2, and store in temp
  mpz_mul(temp.a, m1->a, m2->a);  // temp.a = m1->a * m2->a
  mpz_addmul(temp.a, m1->b, m2->c);  // temp.a += m1->b * m2->c

  mpz_mul(temp.b, m1->a, m2->b);  // temp.b = m1->a * m2->b
  mpz_addmul(temp.b, m1->b, m2->d);  // temp.b += m1->b * m2->d

  mpz_mul(temp.c, m1->c, m2->a);  // temp.c = m1->c * m2->a
  mpz_addmul(temp.c, m1->d, m2->c);  // temp.c += m1->d * m2->c

  mpz_mul(temp.d, m1->c, m2->b);  // temp.d = m1->c * m2->b
  mpz_addmul(temp.d, m1->d, m2->d);  // temp.d += m1->d * m2->d


  // Copy the computed result from temp to the result matrix
  matrixCopy(result, &temp);

  // Clear resources allocated for the temporary matrix
  mpz_clear(temp.a);
  mpz_clear(temp.b);
  mpz_clear(temp.c);
  mpz_clear(temp.d);
}


// Power of a 2x2 matrix
void matrixPower(Matrix2x2 *result, Matrix2x2 *m, int n) {
  Matrix2x2 temp; // Temp matrix
  matrixInit(&temp); // Initialize
  matrixSetIdentity(&temp); // Set to the identity matrix

  // if n is 0, the result is the identity matrix
  if (n == 0) {
    matrixCopy(result, &temp); // Copy the identity matrix to the result
    return;
  }
  // if n is 1, the result is the matrix itself
  else if (n == 1) {
    matrixCopy(result, m); // Copy input matrix to the result
    return;
  }

  // Compute the matrix raised to n/2
  matrixPower(result, m, n / 2);
  // Multiply result of recursion by itself
  matrixMultiply(&temp, result, result);

  // If n is odd, multiply the result by the matrix again
  if (n % 2 != 0) {
    matrixMultiply(result, &temp, m);
  }
  // If n is even, just copy the squared matrix to the result
  else {
    matrixCopy(result, &temp);
  }

  // Clear the mpz_t types
  mpz_clear(temp.a);
  mpz_clear(temp.b);
  mpz_clear(temp.c);
  mpz_clear(temp.d);
}


// nth Fibonacci number using matrix exponentiation
void fibonacci_matrix(mpz_t result, int n) {
  // if n is 0, set the result to 0 and return
  if (n == 0) {
    mpz_set_ui(result, 0); // Set the result to 0
    return;
  }

  Matrix2x2 base, res; // Declare/init base and the result
  matrixInit(&base);
  matrixInit(&res);

  // Setup the base matrix
  mpz_set_ui(base.a, 0); // Set base[0][0] to 0
  mpz_set_ui(base.b, 1); // Set base[0][1] to 1
  mpz_set_ui(base.c, 1); // Set base[1][0] to 1
  mpz_set_ui(base.d, 1); // Set base[1][1] to 1

  // Power of the base matrix to the nth term
  matrixPower(&res, &base, n);

  // n is located at res[0][1] (res.b)
  mpz_set(result, res.b); // Copy res.b to the result

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

/**
 * Prints the given time in a human-readable format.
 * This function displays the time in seconds, milliseconds, microseconds, or
 * nanoseconds based on the input value's magnitude. It provides a convenient
 * way to interpret and display time measurements.
 *
 * @param timeInSeconds The time duration to be printed, in seconds.
 */
void printTime(double timeInSeconds) {
    if (timeInSeconds >= 1.0) {
        printf("%f seconds\n", timeInSeconds);  // Print in seconds if >= 1 second
    } else if (timeInSeconds >= 1e-3) {
        printf("%f milliseconds\n", timeInSeconds * 1e3);  // Convert to milliseconds
    } else if (timeInSeconds >= 1e-6) {
        printf("%f microseconds\n", timeInSeconds * 1e6);  // Convert to microseconds
    } else {
        printf("%f nanoseconds\n", timeInSeconds * 1e9);  // Convert to nanoseconds
    }
}


/**
 * Computes the nth Fibonacci number using dynamic programming.
 * This function uses dynamic programming to efficiently calculate the Fibonacci number.
 * It avoids redundant calculations by storing previously computed Fibonacci numbers in an array.
 *
 * @param result mpz_t variable to store the nth Fibonacci number.
 * @param n The index of the Fibonacci number to be calculated.
 */
void fibonacci_dynamic(mpz_t result, int n) {
    if (n <= 1) {
        mpz_set_ui(result, n);  // Base case: directly set result for n = 0 or 1
        return;
    }

    mpz_t *fib_cache = malloc((n + 1) * sizeof(mpz_t));  // Allocate array for Fibonacci numbers

    mpz_init_set_ui(fib_cache[0], 0);  // Initialize fib[0]
    mpz_init_set_ui(fib_cache[1], 1);  // Initialize fib[1]

    // Calculate Fibonacci numbers iteratively
    for (int i = 2; i <= n; ++i) {
        mpz_init(fib_cache[i]);  // Initialize current element
        mpz_add(fib_cache[i], fib_cache[i - 1], fib_cache[i - 2]);  // fib[i] = fib[i-1] + fib[i-2]
    }

    mpz_set(result, fib_cache[n]);  // Set result to the nth Fibonacci number

    // Clear memory and free the array
    for (int i = 0; i <= n; ++i) {
        mpz_clear(fib_cache[i]);
    }
    free(fib_cache);
}



/**
 * Computes the nth Fibonacci number using an iterative approach.
 * This function iteratively calculates Fibonacci numbers, offering an efficient
 * solution without the need for recursion. It uses three variables to store
 * intermediate Fibonacci numbers.
 *
 * @param result mpz_t variable to store the nth Fibonacci number.
 * @param n The index of the Fibonacci number to be calculated.
 */
void fibonacci_iterative(mpz_t result, int n) {
    if (n <= 1) {
        mpz_set_ui(result, n);  // Directly set result for n = 0 or 1
        return;
    }

    mpz_t a, b, c;
    mpz_init_set_ui(a, 0);  // Initialize a as fib[0]
    mpz_init_set_ui(b, 1);  // Initialize b as fib[1]
    mpz_init(c);            // c will store the next Fibonacci number

    for (int i = 2; i <= n; i++) {
        mpz_add(c, a, b);  // Calculate next Fibonacci number
        mpz_set(a, b);     // Shift a to b's position
        mpz_set(b, c);     // Update b to the new Fibonacci number
    }

    mpz_set(result, b);  // Set the result as the nth Fibonacci number

    mpz_clear(a);  // Free memory allocated for a
    mpz_clear(b);  // Free memory allocated for b
    mpz_clear(c);  // Free memory allocated for c
}


/**
 * Finds an optimal 'n' value for Fibonacci calculation that approximates a target execution time.
 * It employs binary search to find the 'n' value where the recursive calculation of the Fibonacci
 * number is closest to a specified target time. It also measures the CPU time used and the total time taken.
 *
 * @param low The lower bound for the binary search.
 * @param high The upper bound for the binary search.
 * @param target_time_sec The target execution time in seconds.
 * @param cpu_time_used Pointer to store the CPU time used for the calculation.
 * @param total_time_used Pointer to store the total time taken for the process.
 * @param benchResult Pointer to store the Fibonacci number result of the benchmark.
 * @return The optimal 'n' value found.
 */
int find_n_for_target_time_binary(int low, int high, int target_time_sec,
                                  double *cpu_time_used, double *total_time_used,
                                  unsigned long long int *benchResult) {
    int mid;  // Middle index for binary search
    int bestN = low;  // Best 'n' value found so far
    double bestTime = 0.0;  // Best time closest to target time

    *cpu_time_used = 0.0;  // Initialize CPU time used
    clock_t overall_start = clock();  // Start time of overall process

    printf("Finding 'n' value using binary search. Please wait...\n");

    // Binary search to find optimal 'n' value
    while (low <= high) {
        mid = low + (high - low) / 2;
        clock_t start = clock();
        *benchResult = fibRecursiveBench(mid);  // Fibonacci calculation at mid index
        clock_t end = clock();
        double current_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

        printf("\tLow: %d, High: %d, Mid: %d, Time: %f\n", low, high, mid, current_time_used);

        // Adjust search bounds based on current time taken
        if (current_time_used <= target_time_sec + 5) {
            if (current_time_used > bestTime && current_time_used <= target_time_sec + 5) {
                bestN = mid;
                bestTime = current_time_used;
            }
            if (current_time_used < target_time_sec) low = mid + 1;
            else break;
        } else {
            high = mid - 1;
        }
    }

    clock_t overall_end = clock();
    *total_time_used = ((double) (overall_end - overall_start)) / CLOCKS_PER_SEC;
    *cpu_time_used = bestTime;

    printf("'n' found. Using %d for further runs.\n\n", bestN);

    return bestN;  // Return the best 'n' value found
}


/**
 * Prints the results of a Fibonacci number calculation.
 * It displays the method used, the 'n' value, the result, and the execution time.
 *
 * @param method The method used to calculate the Fibonacci number.
 * @param n The index of the Fibonacci number calculated.
 * @param result The calculated Fibonacci number.
 * @param time_taken The time taken to calculate the Fibonacci number.
 */
void printResults(const char* method, int n, mpz_t result, double time_taken) {
    printf("%s Approach (n = %d): ", method, n);
    gmp_printf("Result = %Zd, ", result);
    printTime(time_taken);
}


/**
 * Compares the results of Fibonacci number calculations from different methods.
 * It checks if the results obtained from matrix, iterative, and dynamic methods are consistent.
 *
 * @param matrixResult The result obtained from the matrix method.
 * @param iterativeResult The result obtained from the iterative method.
 * @param dynamicResult The result obtained from the dynamic method.
 */
void compareResults(mpz_t matrixResult, mpz_t iterativeResult, mpz_t dynamicResult) {
    if (mpz_cmp(matrixResult, iterativeResult) == 0 && mpz_cmp(matrixResult, dynamicResult) == 0) {
        printf("The results match.\n\n");
    } else {
        printf("Warning: Inconsistent results between methods.\n\n");
    }
}


/**
 * Benchmarks various methods of calculating Fibonacci numbers.
 * It allows comparison of performance across recursive, iterative, dynamic, and matrix-based approaches.
 * Optionally, it can find an 'n' value for which a recursive calculation takes a target time.
 *
 * @param binary_low Lower bound for binary search.
 * @param binary_high Upper bound for binary search.
 * @param target_time_sec Target time in seconds for the binary search.
 * @param n_value The 'n' value for Fibonacci calculation.
 * @param benchmark_flag Flag to run benchmarking.
 * @param compare_flag Flag to run comparison of methods.
 * @param slow_flag Flag to include slow (recursive) method in comparison.
 * @param results_flag Flag to print the results of calculations.
 */
void benchmark(int binary_low, int binary_high, int target_time_sec, int n_value,
               int benchmark_flag, int compare_flag, int slow_flag, int results_flag) {
    int n;
    double time_taken_recursive = 0.0, time_taken_dynamic = 0.0,
           time_taken_iterative = 0.0, time_taken_matrix = 0.0,
           time_taken_to_find_n = 0.0;
    unsigned long long int benchResult;
    mpz_t result_recursive, result_dynamic, result_iterative, result_matrix;

    // Initialize GMP variables
    mpz_init(result_recursive);
    mpz_init(result_dynamic);
    mpz_init(result_iterative);
    mpz_init(result_matrix);

    // Benchmarking mode
    if (benchmark_flag) {
        // Find 'n' using binary search
        n = find_n_for_target_time_binary(binary_low, binary_high, target_time_sec, &time_taken_recursive, &time_taken_to_find_n, &benchResult);
        printf("Time taken to find n using binary search: ");
        printTime(time_taken_to_find_n);
        printf("\nRecursive Approach (n = %d): ", n);
        printTime(time_taken_recursive);
    } else {
        // Use provided 'n' value
        n = n_value;
    }

    // Perform matrix method and print results
    time_taken_matrix = measureTime(fibonacci_matrix, n, result_matrix);
    printTime(time_taken_matrix);

    if (results_flag) printResults("Matrix", n, result_matrix, time_taken_matrix);

    // Perform iterative method and print results
    if (compare_flag || benchmark_flag) {
        time_taken_iterative = measureTime(fibonacci_iterative, n, result_iterative);
        if (results_flag) printResults("Iterative", n, result_iterative, time_taken_iterative);

        // Perform recursive method if slow_flag is set
        if (slow_flag) {
            time_taken_recursive = measureTime(fibonacci_recursive, n, result_recursive);
            if (results_flag) printResults("Recursive", n, result_recursive, time_taken_recursive);
        }

        // Perform dynamic method and print results
        time_taken_dynamic = measureTime(fibonacci_dynamic, n, result_dynamic);
        if (results_flag) printResults("Dynamic", n, result_dynamic, time_taken_dynamic);

        // Check if results match
        compareResults(result_matrix, result_iterative, result_dynamic);
    }

    // Clear GMP variables
    mpz_clear(result_recursive);
    mpz_clear(result_dynamic);
    mpz_clear(result_iterative);
    mpz_clear(result_matrix);
}

int main(int argc, char *argv[]) {
  int binary_low = 30, binary_high = 55, target_time_sec = 30;
  int benchmark_flag = 0, compare_flag = 0, slow_flag = 0, results_flag = 0;
  unsigned long long n_value = MAX_N;

  processArgs(argc, argv, &binary_low, &binary_high, &target_time_sec, &benchmark_flag, &compare_flag, &slow_flag, &results_flag, &n_value);

  benchmark(binary_low, binary_high, target_time_sec, n_value, benchmark_flag, compare_flag, slow_flag, results_flag);

  printf("This is part 2 of the assignment:\n");
  fibIterativeBench(92);
  fibIterativeBench(100);

  compareFib();

  clearFibMemo();

  return 0;
}

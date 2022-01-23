#include <cstdio>
#include <cstdint>
#include <ctime>
#include "flint/fmpz.h"
#include "flint/flint.h"
#include "flint/fq.h"
#include "omp.h"

#define USE_QUICKSILVER_PRIME
// #define USE_VIRGO_PRIME

#if (defined(USE_QUICKSILVER_PRIME) && defined(USE_VIRGO_PRIME)) || (!defined(USE_QUICKSILVER_PRIME) && !defined(USE_VIRGO_PRIME))
#error "Please select one and only one of the prime numbers."
#endif

void find_first_quadratic_nonresidue(fmpz_t res, fmpz_t p) {
	fmpz_set_ui(res, 2ull);

	while(fmpz_jacobi(res, p) != -1) {
		fmpz_add_ui(res, res, 1ull);
	}
}

void modsqrt(fmpz *out_sqrt, bool* out_is_qr, fmpz *in, int len, fmpz_t qnr, fmpz_t p) {
	int num_of_threads = omp_get_max_threads();
	int chunk_size = (len + num_of_threads - 1) / num_of_threads;

	fmpz* tmp = _fmpz_vec_init(num_of_threads);

	double outer_start_time = omp_get_wtime();
	#pragma omp parallel for default(shared)
	for(int i = 0; i < num_of_threads; i++) {
		int start = i * chunk_size;
		int end = (i + 1) * chunk_size;

		if(end > len) {
			end = len;
		}

		for(int j = start; j < end; j++) {
			fmpz_set(&tmp[i], &in[j]);

			// 1. Check if A is a quad residue.
			bool is_qr = fmpz_jacobi(&tmp[i], p) == 1;

			// 2. If not, A = A * QNR, multiply A with QNR that we found.
			if(!is_qr) {
				fmpz_mul(&tmp[i], &tmp[i], qnr);
			}

			// 3. Set the out_is_qr[j].
			out_is_qr[j] = is_qr;

			// 4. Compute the sqaure root of A.
			fmpz_sqrtmod(&out_sqrt[j], &tmp[i], p);
		}
	}
	double outer_end_time = omp_get_wtime();
	_fmpz_vec_clear(tmp, num_of_threads);

	double time_elapsed = outer_end_time - outer_start_time;
	printf("Total time taken: %.3f seconds for %d entries\n", time_elapsed, len);
}

int main() {
	flint_rand_t state;
	flint_randinit(state);

	// Set the prime number
	fmpz_t p;
#ifdef USE_QUICKSILVER_PRIME
	{
		uint64_t qs_prime;
		qs_prime = (1ULL << 62) - (1ULL << 16) + 1;
		fmpz_init_set_ui(p, qs_prime);
	}
#endif
#ifdef USE_VIRGO_PRIME
	{
		uint64_t virgo_prime;
		virgo_prime = (1ULL << 62) - (1ULL << 26) - (1ULL << 25) + 1;
		fmpz_init_set_ui(p, virgo_prime);
	}
#endif

	// Find the first QNR
	fmpz_t qnr;
	fmpz_init(qnr);

	printf("Number of available threads: %d\n", omp_get_max_threads());

	find_first_quadratic_nonresidue(qnr, p);

	printf("The first QNR found is: ");
	fmpz_print(qnr);
	printf("\n");

	// Generate a list of random numbers
	int len = 1ULL << 20; // 2^20 = 1048576
	printf("Benchmark %d square roots \n", len);
	fmpz* vec_in = _fmpz_vec_init(len);

	for(int i = 0; i < len; i++) {
		fmpz_randm(&vec_in[i], state, p);
	}

	// Allocate the space for the output
	fmpz* vec_out = _fmpz_vec_init(len);
	bool* vec_is_qr = static_cast<bool *>(malloc(sizeof(bool) * len));

	// Call the main program
	modsqrt(vec_out, vec_is_qr, vec_in, len, qnr, p);

	_fmpz_vec_clear(vec_in, len);
	_fmpz_vec_clear(vec_out, len);
	free(vec_is_qr);
	fmpz_clear(p);
	flint_randclear(state);

	return 0;
}
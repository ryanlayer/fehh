#include "fehh.h"

uint32_t pop_count(uint32_t *b, uint32_t n_words)
{
    uint32_t i, count=0;

    for (i = 0; i < n_words; ++i) {
        uint32_t x = b[i];
#ifndef __SSE4_2__
        for (; x; count++)
            x &= x-1;
#else
        count += _mm_popcnt_u32(x);
#endif
    }

    return count;
}

void and(uint32_t *r, uint32_t *a, uint32_t *b, uint32_t n_words)
{
    uint32_t i;

    for (i = 0; i < n_words; ++i) 
        r[i] = a[i] & b[i];
}

void set(uint32_t *b, uint32_t i)
{
    b[i/32] |= 1 << (i%32);
}

void set_bit_arrays(char *h, uint32_t *a0, uint32_t *a1, uint32_t num_samples)
{
    uint32_t i;
    for (i = 0; i < num_samples; ++i) {
        if (h[i*2] == '0') 
            set(a0,i);
        else
            set(a1,i);
    }
}



/**
 * @brief Compute the EHH considering the next marker.
 *
 * This function must be called in order (either up stream or down stream)
 * starting with the target locus to ensure that the disjoint sets in A are
 * correct.
 *
 * @param
 * @param A disjoint set 
 * @param len_A number of elements in the disjoint set
 * @param C1 mask defining the samples with the core haplotype at the targe
 * @param B list of bit arrays for variants that have been seen
 * @param B_i index within B of the marker to be considered here
 * @param S number of unique haplotypes up to this point
 * @param n_words number of words required for one node in the disjoint set
 * @param n_bytes number of bytes required for one node in the disjoint set
 * @param B0_c pre allocated space (n_bytes long)
 * @param B1_c pre allocated space  (n_bytes long)
 * @param R pre allocated space (n_bytes long)
 * @param t_A pre allocated space (num_samples*n_bytes long)
 */
long double ehh_step(uint32_t **A,
                     uint32_t *len_A,
                     uint32_t *C1,
                     uint32_t num_C1,
                     uint32_t *B,
                     uint32_t i,
                     uint32_t *S,
                     uint32_t n_words,
                     uint32_t n_bytes,
                     uint32_t *B0_c,
                     uint32_t *B1_c,
                     uint32_t *R,
                     uint32_t **t_A)
{
    uint32_t *B0, *B1;
    B0 = B + (i*2*n_words);
    B1 = B + ((i*2+1)*n_words);

    and(B0_c, B0, C1, n_words);
    and(B1_c, B1, C1, n_words);

    uint32_t pop, j, len_t_A = 0;
    for (j = 0; j < *len_A; ++j) {
        and(R, (*A) + (j*n_words), B0_c, n_words);
        pop = pop_count(R, n_words);

        if (pop == 1) 
            *S += 1;
        else if (pop > 0) {
            memcpy((*t_A)+len_t_A*n_words, R, n_bytes);
            len_t_A += 1;
        }

        and(R, (*A) + (j*n_words), B1_c, n_words);
        pop = pop_count(R, n_words);

        if (pop == 1) 
            *S += 1;
        else if (pop > 0) {
            memcpy((*t_A)+len_t_A*n_words, R, n_bytes);
            len_t_A += 1;
        }
    }

    *len_A = len_t_A;
    uint32_t *swap_A = *A;
    *A = *t_A;
    *t_A = swap_A;

    long double ehh = 0;
    for (j = 0; j < *len_A; ++j)
        ehh += n_choose_2[pop_count((*A)+(j*n_words), n_words)] /
               n_choose_2[num_C1];

    return ehh;
}

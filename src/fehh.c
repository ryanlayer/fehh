#include "fehh.h"

uint32_t pop_count(uint32_t *b, uint32_t num_words)
{
    uint32_t i, count=0;

    for (i = 0; i < num_words; ++i) {
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

void and(uint32_t *r, uint32_t *a, uint32_t *b, uint32_t num_words)
{
    uint32_t i;

    for (i = 0; i < num_words; ++i) 
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
 * @param num_words number of words required for one node in the disjoint set
 * @param num_bytes number of bytes required for one node in the disjoint set
 * @param B0_c pre allocated space (num_bytes long)
 * @param B1_c pre allocated space  (num_bytes long)
 * @param R pre allocated space (num_bytes long)
 * @param t_A pre allocated space (num_samples*num_bytes long)
 */
long double ehh_step(uint32_t **A,
                     uint32_t *len_A,
                     //uint32_t *C1,
                     uint32_t num_C1,
                     uint32_t *B,
                     uint32_t i,
                     //uint32_t *S,
                     uint32_t num_words,
                     uint32_t num_bytes,
                     //uint32_t *B0_c,
                     //uint32_t *B1_c,
                     uint32_t *R,
                     uint32_t **t_A)
{
    uint32_t *B0, *B1;
    B0 = B + (i*2*num_words);
    B1 = B + ((i*2+1)*num_words);

    //and(B0_c, B0, C1, num_words);
    //and(B1_c, B1, C1, num_words);

    uint32_t pop, j, len_t_A = 0;
    for (j = 0; j < *len_A; ++j) {
        //and(R, (*A) + (j*num_words), B0_c, num_words);
        and(R, (*A) + (j*num_words), B0, num_words);
        pop = pop_count(R, num_words);

        /*
        if (pop == 1) 
            *S += 1;
        else if (pop > 0) {
        */
        if (pop > 0) {
            memcpy((*t_A)+len_t_A*num_words, R, num_bytes);
            len_t_A += 1;
        }

        //and(R, (*A) + (j*num_words), B1_c, num_words);
        and(R, (*A) + (j*num_words), B1, num_words);
        pop = pop_count(R, num_words);

        /*
        if (pop == 1) 
            *S += 1;
        else if (pop > 0) {
        */
        if (pop > 0) {
            memcpy((*t_A)+len_t_A*num_words, R, num_bytes);
            len_t_A += 1;
        }
    }

    *len_A = len_t_A;
    uint32_t *swap_A = *A;
    *A = *t_A;
    *t_A = swap_A;

    long double ehh = 0;
    for (j = 0; j < *len_A; ++j)
        ehh += n_choose_2[pop_count((*A)+(j*num_words), num_words)] /
               n_choose_2[num_C1];

    return ehh;
}


ssize_t push_B(uint32_t *B,
               uint32_t **B0,
               uint32_t **B1,
               uint32_t *B_i,
               uint32_t num_words,
               uint32_t num_samples,
               char **locus_name,
               char **line, 
               size_t *linecap,
               FILE *fp)
{
    ssize_t linelen;
    char *p, *h;
    if (((linelen = getline(line, linecap, fp)) > 0)) {

        // skip some fields
        p = strtok(*line, " ");
        *locus_name = strtok(NULL, " ");
        // skip some fields
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");
        h = p + strlen(p) + 1;

        // parse the variants into bit arays 
        *B0 = B + ((*B_i)*2*num_words);
        *B1 = B + (((*B_i)*2+1)*num_words);
        set_bit_arrays(h, *B0, *B1, num_samples);
        *B_i += 1;
    }

    return linelen;
}

long double ihh_back(uint32_t locus_i,
                     uint32_t *B,
                     uint32_t **A0,
                     uint32_t **A1,
                     uint32_t num_words,
                     uint32_t num_bytes,
                     uint32_t *R,
                     double *genetic_pos,
                     uint32_t *physical_pos,
                     uint32_t **t_A) 
{
#if 1
    // work down stream with just C1
    double locus_genetic_pos = genetic_pos[locus_i];
    uint32_t locus_physical_pos = physical_pos[locus_i];

    uint32_t *B0 = B + (locus_i*2*num_words);
    uint32_t *B1 = B + ((locus_i*2+1)*num_words);
    uint32_t num_C0 = pop_count(B0, num_words);
    uint32_t num_C1 = pop_count(B1, num_words);

    memcpy(*A1, B1, num_bytes);
    uint32_t len_A1 = 1;
    memcpy(*A0, B0, num_bytes);
    uint32_t len_A0 = 1;

    uint32_t i = locus_i - 1;
    long double e0 = 1, e1 = 1;
    while ( (i > 0) && (e0 > 0.05) && (e1 > 0.05) ) {
        printf("-%u\t%f\t", locus_physical_pos-physical_pos[i],
                           genetic_pos[i]-locus_genetic_pos);
        e1 = ehh_step(A1,
                      &len_A1,
                      //C1,
                      num_C1,
                      B,
                      i,
                      //&S,
                      num_words,
                      num_bytes,
                      //B0_c,
                      //B1_c,
                      R,
                      t_A);

        printf("%Lf\t", e1);

        e0 = ehh_step(A0,
                      &len_A0,
                      //C0,
                      num_C0,
                      B,
                      i,
                      //&S,
                      num_words,
                      num_bytes,
                      //B0_c,
                      //B1_c,
                      R,
                      t_A);

        printf("%Lf\n", e0);

        i -= 1;
    }
    printf("--\n");
#endif
    return 0;
}

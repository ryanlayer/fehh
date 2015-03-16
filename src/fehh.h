#ifndef __FEHH_H__
#define __FEHH_H__

#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "n_choose_2.h"


uint32_t pop_count(uint32_t *b, uint32_t num_words);

void and(uint32_t *r, uint32_t *a, uint32_t *b, uint32_t num_words);

void set(uint32_t *b, uint32_t i);

void set_bit_arrays(char *h, uint32_t *a0, uint32_t *a1, uint32_t num_samples);

long double ehh_step(uint32_t **A,
                     uint32_t *len_A,
                     //uint32_t *C1,
                     uint32_t num_C1,
                     uint32_t *B,
                     uint32_t B_i,
                     //uint32_t *S,
                     uint32_t num_words,
                     uint32_t num_bytes,
                     //uint32_t *B0_c,
                     //uint32_t *B1_c,
                     uint32_t *R,
                     uint32_t **t_A);

ssize_t push_B(uint32_t *B,
               uint32_t *B_i,
               uint32_t num_words,
               uint32_t num_samples,
               char **locus_name,
               char **line, 
               size_t *linecap,
               FILE *fp);

long double ihh_back(uint32_t locus_i,
                     uint32_t *B,
                     uint32_t **A0,
                     uint32_t **A1,
                     uint32_t num_words,
                     uint32_t num_bytes,
                     uint32_t *R,
                     double *genetic_pos,
                     uint32_t *physical_pos,
                     uint32_t **t_A);

long double ihh_forward(uint32_t locus_i,
                        uint32_t *B,
                        uint32_t *B_i,
                        FILE *fp,
                        char **line,
                        uint32_t **A0,
                        uint32_t **A1,
                        uint32_t num_words,
                        uint32_t num_bytes,
                        uint32_t num_samples,
                        uint32_t *R,
                        double *genetic_pos,
                        uint32_t *physical_pos,
                        uint32_t **t_A);
#endif

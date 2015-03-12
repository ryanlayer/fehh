#include "fehh.h"

//{{{int usage(char *prog)
int usage(char *prog)
{
    fprintf(stderr, "usage:\t%s <options>\n"
                    "\t\t-f\tfile name\n"
                    "\t\t-v\tnumber of variants\n"
                    "\t\t-s\tnumber of samples\n"
                    "\t\t-l\tlocus ID\n"
                    "\t\t-h\tHelp\n",
                    prog);
    return 1;
}
//}}}

int main(int argc, char **argv)
{
    char *file_name, *locus_id;
    uint32_t num_samples, num_variants;

    //{{{ input params
    int s_is_set = 0,
        f_is_set = 0,
        l_is_set = 0,
        v_is_set = 0,
        h_is_set = 0;

    int c;

    while ((c = getopt (argc, argv, "hs:f:v:l:")) != -1) {
        switch (c) {
        case 'h':
            h_is_set = 1;
            break;
        case 'v':
            v_is_set = 1;
            num_variants = atoi(optarg);
            break;
        case 's':
            s_is_set = 1;
            num_samples = atoi(optarg);
            break;
        case 'f':
            f_is_set = 1;
            file_name = optarg;
            break;
        case 'l':
            l_is_set = 1;
            locus_id = optarg;
            break;
        case '?':
            if ( (optopt == 'f') || 
                 (optopt == 'v') || 
                 (optopt == 's') || 
                 (optopt == 'l') )
                fprintf (stderr, "Option -%c requires an argument.\n",
                         optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        default:
            return usage(argv[0]);
        }
    }

    if (h_is_set)
        return usage(argv[0]);

    if (!f_is_set) {
        fprintf(stderr, "Text file name is not set.\n");
        return usage(argv[0]);
    }

    if (!v_is_set) {
        fprintf(stderr, "Number of variants is not set.\n");
        return usage(argv[0]);
    }

    if (!l_is_set) {
        fprintf(stderr, "Locus ID is not set.\n");
        return usage(argv[0]);
    }

    if (!s_is_set) {
        fprintf(stderr, "Number of samples is not set.\n");
        return usage(argv[0]);
    }
    //}}}

    /*
    htsFile *fp = hts_open(file_name, "rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    uint32_t num_samples = bcf_hdr_nsamples(hdr);

    fprintf(stderr, "%u\n", num_samples);

    bcf1_t *rec = bcf_init();
    int32_t *gt = NULL;
    int ntmp = 0;
    while (bcf_read(fp, hdr, rec) >= 0) {
        bcf_unpack(rec, BCF_UN_ALL);
        int num_gts = bcf_get_genotypes(hdr, rec, &gt, &ntmp);
    }

    bcf_destroy(rec);
    */

    uint32_t n_words = num_samples / 32 + 1;
    uint32_t n_bytes =  n_words*sizeof(uint32_t);
    uint32_t *A0 = (uint32_t *)malloc(num_samples*n_bytes);
    uint32_t *A1 = (uint32_t *)malloc(num_samples*n_bytes);
    uint32_t *t_A = (uint32_t *)malloc(num_samples*n_bytes);
    uint32_t *swap_A;

    uint32_t *B0, *B1;

    // B will hold B0 and B1 for each variant
    uint32_t *B = (uint32_t *)calloc(num_variants*2, n_bytes);

    // C0 and C1 hold the masks that define the samples with and without
    // the core haplotype
    uint32_t *C0 = (uint32_t *)calloc(1, n_bytes);
    uint32_t *C1 = (uint32_t *)calloc(1, n_bytes);
    uint32_t num_C0, num_C1;

    uint32_t *B0_c = (uint32_t *)calloc(1, n_bytes);
    uint32_t *B1_c = (uint32_t *)calloc(1, n_bytes);

    uint32_t *R = (uint32_t *)malloc(n_bytes);

    FILE *fp = fopen(file_name, "r");
    
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", file_name);
        exit(1);
    } 

    char *p, *h, *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    uint32_t i, j, pop, len_A0 = 0, len_A1 = 0, len_t_A = 0, S = 0, B_i = 0;
    int locus_found = 0;


    // scan through the file until we either hit the locus of interest 
    // or reach the end
    while ( !(locus_found) && 
            ((linelen = getline(&line, &linecap, fp)) > 0)) {

        // skip some fields
        p = strtok(line, " ");
        p = strtok(NULL, " ");

        // if we hit the locus of interest, set C0 and C1 as the masks
        // for the ancestral and derived alleles
        if (strcmp(p, locus_id) == 0) {
            fprintf(stderr, "%s\n", p);
            set_bit_arrays(h, C0, C1, num_samples);
            num_C0 = pop_count(C0, n_words);
            num_C1 = pop_count(C1, n_words);
            locus_found = 1;
        }

        // skip some fields
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");
        h = p + strlen(p) + 1;

        // parse the variants into bit arays 
        B0 = B + (B_i*2*n_words);
        B1 = B + ((B_i*2+1)*n_words);
        set_bit_arrays(h, B0, B1, num_samples);

        B_i += 1;
    }

    // work down stream with just C1
    i = B_i - 1; //move back one to start at the targe locus

    B0 = B + (i*2*n_words);
    B1 = B + ((i*2+1)*n_words);

    and(R, B0, C1, n_words);
    memcpy(A1, R, n_bytes);
    and(R, B1, C1, n_words);
    memcpy(A1+n_words, R, n_bytes);
    len_A1 = 2;

    and(R, B0, C0, n_words);
    memcpy(A0, R, n_bytes);
    and(R, B1, C0, n_words);
    memcpy(A0+n_words, R, n_bytes);
    len_A0 = 2;


    for (i = B_i - 2; i > B_i - 12; i-=1) {
        fprintf(stderr, "i:%u\t", i);

#if 1
        long double e1 = ehh_step(&A1,
                                  &len_A1,
                                  C1,
                                  num_C1,
                                  B,
                                  i,
                                  &S,
                                  n_words,
                                  n_bytes,
                                  B0_c,
                                  B1_c,
                                  R,
                                  &t_A);
        printf("%Lf\t", e1);

        long double e0 = ehh_step(&A0,
                                  &len_A0,
                                  C0,
                                  num_C0,
                                  B,
                                  i,
                                  &S,
                                  n_words,
                                  n_bytes,
                                  B0_c,
                                  B1_c,
                                  R,
                                  &t_A);
        printf("%Lf\n", e0);
#endif

#if 0
        fprintf(stderr, "len_A:%u\t", len_A);
        B0 = B + (i*2*n_words);
        B1 = B + ((i*2+1)*n_words);

        and(B0_c, B0, C1, n_words);
        and(B1_c, B1, C1, n_words);

        len_t_A = 0;
        for (j = 0; j < len_A; ++j) {
            and(R, A + (j*n_words), B0_c, n_words);
            pop = pop_count(R, n_words);

            if (pop == 1) 
                S += 1;
            else if (pop > 0) {
                memcpy(t_A+len_t_A*n_words, R, n_bytes);
                len_t_A += 1;
            }

            and(R, A + (j*n_words), B1_c, n_words);
            pop = pop_count(R, n_words);

            if (pop == 1) 
                S += 1;
            else if (pop > 0) {
                memcpy(t_A+len_t_A*n_words, R, n_bytes);
                len_t_A += 1;
            }
        }

        len_A = len_t_A;
        swap_A = A;
        A = t_A;
        t_A = swap_A;

        //long double n_choose_2[100001];
        long double ehh_d = 0;
        for (j = 0; j < len_A; ++j) {
            printf("%u ", pop_count(A+(j*n_words), n_words));
            ehh_d += n_choose_2[pop_count(A+(j*n_words), n_words)] /
                     n_choose_2[num_C1];
        }
        for (j = 0; j < S; ++j) {
            printf("1 ");
        }
        printf(" %Lf\n", ehh_d);
#endif
    }
    

#if 0
        if (empty == 0) {
            empty = 1;
            memcpy(A, B0, n_bytes);
            memcpy(A+n_words, B1, n_bytes);
            len_A = 2;
        } else {
            len_t_A = 0;

            for (i = 0; i < len_A; ++i) {

                and(R, A + (i*n_words), B0, n_words);
                pop = pop_count(R, n_words);

                if (pop == 1) 
                    S += 1;
                else if (pop > 0) {
                    memcpy(t_A+len_t_A*n_words, R, n_bytes);
                    len_t_A += 1;
                }

                and(R, A + (i*n_words), B1, n_words);
                pop = pop_count(R, n_words);

                if (pop == 1) 
                    S += 1;
                else if (pop > 0) {
                    memcpy(t_A+len_t_A*n_words, R, n_bytes);
                    len_t_A += 1;
                }
            }

            len_A = len_t_A;
            swap_A = A;
            A = t_A;
            t_A = swap_A;
        }
    }
    
    uint32_t sum = 0;
    for (i = 0; i < len_A; ++i) {
        printf("%u\n", pop_count(A+(i*n_words), n_words));
        sum += pop_count(A+(i*n_words), n_words);
    }

    sum += S;
    printf("S:%u\n", S);
    printf("sum:%u\n", sum);

    free(line);
    fclose(fp);
#endif

    return 0;
}

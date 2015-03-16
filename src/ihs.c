#include "fehh.h"
#include <tgmath.h>

#define MIN(a,b) \
       ({ __typeof__ (a) _a = (a); \
               __typeof__ (b) _b = (b); \
             _a < _b ? _a : _b; })

//{{{int usage(char *prog)
int usage(char *prog)
{
    fprintf(stderr, "usage:\t%s <options>\n"
                    "\t\t-f\ttped file name\n"
                    "\t\t-m\tmap file name\n"
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
    char *map_file_name, *tped_file_name, *locus_id;
    uint32_t num_samples, num_variants;

    //{{{ input params
    int s_is_set = 0,
        m_is_set = 0,
        f_is_set = 0,
        v_is_set = 0,
        h_is_set = 0;

    int c;

    while ((c = getopt (argc, argv, "hs:f:v:m:")) != -1) {
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
            tped_file_name = optarg;
            break;
        case 'm':
            m_is_set = 1;
            map_file_name = optarg;
            break;
        case '?':
            if ( (optopt == 'f') || 
                 (optopt == 'v') || 
                 (optopt == 's') || 
                 (optopt == 'm') )
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

    if (!s_is_set) {
        fprintf(stderr, "Number of samples is not set.\n");
        return usage(argv[0]);
    }

    if (!m_is_set) {
        fprintf(stderr, "Map file name not set.\n");
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

    //{{{ get genetic_pos and physical_pos from the map file
    FILE *fp = fopen(map_file_name, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", map_file_name);
        exit(1);
    } 

    char *p, *line = NULL;
    size_t linecap = 0;
    ssize_t linelen;
    uint32_t i = 0;
    double *genetic_pos = (double *)malloc(num_variants * sizeof(double));
    uint32_t *physical_pos = (uint32_t *)
            malloc(num_variants * sizeof(uint32_t));
    while ( (linelen = getline(&line, &linecap, fp)) > 0 ) {
        p = strtok(line, " "); // chr
        p = strtok(NULL, " "); // name
        p = strtok(NULL, " "); // genetic pos
        genetic_pos[i] = atof(p);
        p = strtok(NULL, " "); // physical pos
        physical_pos[i] = atoi(p);
        i+=1;
    }
    fclose(fp);
    //}}}

    uint32_t maf_threshold = num_samples * 0.05;
    uint32_t num_words = num_samples / 32 + 1;
    uint32_t num_bytes =  num_words*sizeof(uint32_t);
    uint32_t *A0 = (uint32_t *)malloc(num_samples*num_bytes);
    uint32_t *A1 = (uint32_t *)malloc(num_samples*num_bytes);
    uint32_t *t_A = (uint32_t *)malloc(num_samples*num_bytes);
    uint32_t *swap_A;

    uint32_t *B0, *B1;

    // B will hold B0 and B1 for each variant
    uint32_t *B = (uint32_t *)calloc(num_variants*2, num_bytes);

    uint32_t num_C0, num_C1;

    uint32_t *R = (uint32_t *)malloc(num_bytes);

    fp = fopen(tped_file_name, "r");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", tped_file_name);
        exit(1);
    } 

    char *h, *locus_name;
    uint32_t af, j, pop, 
            len_A0 = 0, len_A1 = 0, len_t_A = 0, S = 0, 
            B_i = 0, L_i = 0;

    for ( L_i = 0; L_i < num_variants; ++L_i) {
        if (L_i == B_i)
            linelen = push_B(B,
                             &B_i,
                             num_words,
                             num_samples,
                             &locus_name,
                             &line, 
                             &linecap,
                             fp);
        if (linelen < 1)
            break;

        B0 = B + (L_i*2*num_words);
        B1 = B + ((L_i*2+1)*num_words);
 
        af = pop_count(B1, num_words);

        //af = MIN(af, num_samples - af);

        long double e0_d, e1_d, e0_u, e1_u;
        if (af >= maf_threshold) {
            // scan down stream (direction = -1)
            printf("%d\t",  physical_pos[L_i]);
            long double r0 = ihh_step(-1,
                                      &e0_d,
                                      &e1_d,
                                      L_i,
                                      B,
                                      &B_i,
                                      fp,
                                      &line,
                                      &A0,
                                      &A1,
                                      num_words,
                                      num_bytes,
                                      num_samples,
                                      R,
                                      genetic_pos,
                                      physical_pos,
                                      &t_A);

            // scan up stream (direction = 1)
            long double r1 = ihh_step(1,
                                      &e0_u,
                                      &e1_u,
                                      L_i,
                                      B,
                                      &B_i,
                                      fp,
                                      &line,
                                      &A0,
                                      &A1,
                                      num_words,
                                      num_bytes,
                                      num_samples,
                                      R,
                                      genetic_pos,
                                      physical_pos,
                                      &t_A);

            long double iHH_0 = e0_d + e0_u;
            long double iHH_1 = e1_d + e1_u;
            printf("%f\t%Lf\t%Lf\t%Lf\n", 
                    ((float)af)/((float)num_samples),
                    iHH_1,
                    iHH_0,
                    logl(iHH_1/iHH_0));
        }
    }
    return 0;
}

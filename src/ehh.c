#include "fehh.h"

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
        l_is_set = 0,
        v_is_set = 0,
        h_is_set = 0;

    int c;

    while ((c = getopt (argc, argv, "hs:f:v:l:m:")) != -1) {
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
        case 'l':
            l_is_set = 1;
            locus_id = optarg;
            break;
        case 'm':
            m_is_set = 1;
            map_file_name = optarg;
            break;
        case '?':
            if ( (optopt == 'f') || 
                 (optopt == 'v') || 
                 (optopt == 's') || 
                 (optopt == 'l') || 
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

    if (!l_is_set) {
        fprintf(stderr, "Locus ID is not set.\n");
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

    uint32_t num_words = num_samples / 32 + 1;
    uint32_t num_bytes =  num_words*sizeof(uint32_t);
    uint32_t *A0 = (uint32_t *)malloc(num_samples*num_bytes);
    uint32_t *A1 = (uint32_t *)malloc(num_samples*num_bytes);
    uint32_t *t_A = (uint32_t *)malloc(num_samples*num_bytes);
    uint32_t *swap_A;

    uint32_t *B0, *B1;

    // B will hold B0 and B1 for each variant
    uint32_t *B = (uint32_t *)calloc(num_variants*2, num_bytes);

    // C0 and C1 hold the masks that define the samples with and without
    // the core haplotype
    //uint32_t *C0 = (uint32_t *)calloc(1, num_bytes);
    //uint32_t *C1 = (uint32_t *)calloc(1, num_bytes);
    uint32_t num_C0, num_C1;

    //uint32_t *B0_c = (uint32_t *)calloc(1, num_bytes);
    //uint32_t *B1_c = (uint32_t *)calloc(1, num_bytes);

    uint32_t *R = (uint32_t *)malloc(num_bytes);

    //FILE *fp = fopen(tped_file_name, "r");
    fp = fopen(tped_file_name, "r");
    
    if (fp == NULL) {
        fprintf(stderr, "ERROR: Could not open %s\n", tped_file_name);
        exit(1);
    } 

    char *h;
    uint32_t j, pop, len_A0 = 0, len_A1 = 0, len_t_A = 0, S = 0, B_i = 0;
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
            //set_bit_arrays(h, C0, C1, num_samples);
            //num_C0 = pop_count(C0, num_words);
            //num_C1 = pop_count(C1, num_words);
            locus_found = 1;
        }

        // skip some fields
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");
        h = p + strlen(p) + 1;

        // parse the variants into bit arays 
        B0 = B + (B_i*2*num_words);
        B1 = B + ((B_i*2+1)*num_words);
        set_bit_arrays(h, B0, B1, num_samples);

        B_i += 1;
    }

    uint32_t locus_i = B_i - 1; //move back one to start at the targe locus
    double locus_genetic_pos = genetic_pos[locus_i];
    uint32_t locus_physical_pos = physical_pos[locus_i];
    long double e0 = 1, e1 = 1;
#if 0
    long double r = ihh_back(locus_i,
                             B,
                             &A0,
                             &A1,
                             num_words,
                             num_bytes,
                             R,
                             genetic_pos,
                             physical_pos,
                             &t_A);
#endif
#if 1
    // work down stream with just C1
    B0 = B + (locus_i*2*num_words);
    B1 = B + ((locus_i*2+1)*num_words);
    num_C0 = pop_count(B0, num_words);
    num_C1 = pop_count(B1, num_words);

    /*
    and(R, B0, C1, num_words);
    memcpy(A1, R, num_bytes);
    and(R, B1, C1, num_words);
    memcpy(A1+num_words, R, num_bytes);
    len_A1 = 2;

    and(R, B0, C0, num_words);
    memcpy(A0, R, num_bytes);
    and(R, B1, C0, num_words);
    memcpy(A0+num_words, R, num_bytes);
    len_A0 = 2;
    */
    memcpy(A1, B1, num_bytes);
    len_A1 = 1;
    memcpy(A0, B0, num_bytes);
    len_A0 = 1;

    i = locus_i - 1;
    //while ( (i > locus_i - 150) ) {
    while ( (i > 0) && (e0 > 0.05) && (e1 > 0.05) ) {
        printf("-%u\t%f\t", locus_physical_pos-physical_pos[i],
                           genetic_pos[i]-locus_genetic_pos);
        e1 = ehh_step(&A1,
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
                      &t_A);

        printf("%Lf\t", e1);

        e0 = ehh_step(&A0,
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
                      &t_A);

        printf("%Lf\n", e0);

        i -= 1;
    }
#endif
    printf("--\n");

    // work upstream
    B0 = B + (locus_i*2*num_words);
    B1 = B + ((locus_i*2+1)*num_words);
    num_C0 = pop_count(B0, num_words);
    num_C1 = pop_count(B1, num_words);

    /*
    and(R, B0, C1, num_words);
    memcpy(A1, R, num_bytes);
    and(R, B1, C1, num_words);
    memcpy(A1+num_words, R, num_bytes);
    len_A1 = 2;

    and(R, B0, C0, num_words);
    memcpy(A0, R, num_bytes);
    and(R, B1, C0, num_words);
    memcpy(A0+num_words, R, num_bytes);
    len_A0 = 2;
    */

    memcpy(A1, B1, num_bytes);
    len_A1 = 1;
    memcpy(A0, B0, num_bytes);
    len_A0 = 1;


    i = locus_i + 1;
    //while ( ((linelen = getline(&line, &linecap, fp)) > 0) &&
            //(i < locus_i + 150)) {
    //long double e0 = 1, e1 = 1;
    e0 = 1; 
    e1 = 1;
    while ( ((linelen = getline(&line, &linecap, fp)) > 0) &&
            (e0 > 0.05) && (e1 > 0.05) ) {
        // skip some fields
        p = strtok(line, " ");
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");
        h = p + strlen(p) + 1;

        // parse the variants into bit arays 
        B0 = B + (B_i*2*num_words);
        B1 = B + ((B_i*2+1)*num_words);
        set_bit_arrays(h, B0, B1, num_samples);

        B_i += 1;

        printf("%u\t%f\t", physical_pos[i]-locus_physical_pos,
                           genetic_pos[i]-locus_genetic_pos);
        e1 = ehh_step(&A1,
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
                      &t_A);

        printf("%Lf\t", e1);

        e0 = ehh_step(&A0,
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
                      &t_A);

        printf("%Lf\n", e0);

        i += 1;
    }

    return 0;
}

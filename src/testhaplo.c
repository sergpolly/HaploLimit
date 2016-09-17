#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "glpk.h"


#define BUFFER_SIZE 1024
#define LABEL_SIZE 32
#define ALL_HAPLO 0



// description of haplotypes: index and bounds ...
typedef struct {
    int index;
    double lower;
    double upper;
} HBound;



// description of the population ...
typedef struct {
    // alleles stats ...
    int loci;
    int haplotypes;
    int alleles_total;
    int independent_alleles; // added for clarity 
    int *allelicity;
    double *freqs;
    char **allele_labels;
    int **haplotype_description;
    // linkages ...
    int linkages;
    int *link_alleles_number; // used to be link_loci_number
    int *link_loci_number; // #_linked_loci <= #_linked_alleles (!!!)
    double *link_freq;
    int link_array_length;
    int *link_loci_indices;
    int *link_alleles_indices;
    //results ...
    HBound *hbounds;
} PopDes;



// Coordinate format of Sparse Matrix for the constraint matrix ...
typedef struct {
   int nonzero, rows, cols;
   int *irow;
   int *jcol;
   double *mval; // should become float afterwards ...
} SparseMat;


// reads input file with the name fname and tries to parse its contents filling in the population structure
// input file format assumed to be strict and any deviations from that format lead to the program termination.
// allele labels parsing attempted (optional)
void load_allele_frequency_file_safe(PopDes *population, const char *fname);

// // reads input file with the name fname and tries to parse its contents filling in the population structure
// // input file format assumed to be strict, deviations in the format lead to undefined behavior. clean code.
// // allele labels are not parsed - behavior undefined if labels are provided.
// void load_allele_frequency_file(PopDes *population, const char *fname);

// partial product of array elements starting from start(exclusive) till length-1 (C-style 0-numbering assumed);
// it'll return 1 if start == length-1 ...
int get_partial_product(int *array, int start, int stop);

// table of haplotypes is generated; loci are in columns, whereas rows are haplotypes
// alleles are 0-numbered (C-style).
void generate_haplotypes(PopDes *population);

// multiple alleles of the same locus can be involved in a linkage,
// implying that the combination of haplotypes involving those alleles sums up at a given frequency.
// this little function is needed to get the number of different loci involved in the linkage
// (which is <= number of involved alleles)
void generate_linked_loci_number(PopDes *population);

// generate constraint matrix, given population information 
void generate_constraint_matrix_linkage(PopDes *population, SparseMat *matrix);

// free memory allocated for the internals of the PopDes structure
// bunch of 2D and 1D arrays to deal with
void free_population(PopDes *population);

// free memory allocated for the internals of the SparseMat structure
// bunch of 1D arrays to deal with
void free_sparse_matrix(SparseMat *matrix);

// print all possible haplotypes using labels(optional) and display haplotype frequency limits(optional)
// void print_haplotypes_and_bounds(PopDes *population, int print_labels, int print_bounds, int num_haplo_print);
void print_haplotypes_and_bounds(PopDes *population, int print_labels, int print_bounds, int num_haplo_print, FILE *fp, const char *fname);

// simplex method for all haplos (both MIN and MAX) wrapped in a single function call 
void run_simplex_all_haplos(PopDes *population, glp_prob *local_lp, glp_smcp *local_parm);

// print matrix of constraints  
void print_constraints_matrix(SparseMat *matrix);

// comparison function for HBound that is used in qsort ...
int compare_hbounds(const void *, const void *);
// void print_union(HBound);

// function that calculates the product of haplo.-range ratios, either ALL of them or a number of TOP ones ...
double get_ratio_product(PopDes *population, HBound *hbounds_original, int fixed_haplo_index, int product_number);

// function sprintfs haplotype description into 'haplo_descr' for a given 'haplo_index' 
void get_haplotype_description(char *haplo_descr, PopDes *population, int haplo_index);



int main(int argc, char **argv) {


    PopDes *pop;
    SparseMat *mat;
    glp_prob *lp;
    HBound *hbounds_prescan; // we'd need that for scanning version of the program ...
    HBound *hbounds_prescan_sorted; // we'd need that for scanning version of the program ...
    // initialized later on ...


    pop = (PopDes *)malloc(sizeof(PopDes)); // allocate population description ...
    mat = (SparseMat *)malloc(sizeof(SparseMat)); // allocate constraint matrix structure ...


    load_allele_frequency_file_safe(pop, argv[1]);
    generate_haplotypes(pop);
    generate_constraint_matrix_linkage(pop, mat);

    // printf("HAPLOTYPES:\n\n");
    // print_haplotypes_and_bounds(pop, 1, 0);
    // printf("CONSTRAINTS MATRIX:\n\n");
    // print_constraints_matrix(mat);

    // launchnig and allocating the glpk ...
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_ERR;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "sample");

    // number of constraints:
    // allelicity_total + 1 (haplotype freqs normalization) + #_of_linkages ...
    int constraints = pop->independent_alleles + 1 + pop->linkages;
    // that was the same as mat->irow ...

    // submit number of constraints to glpk ...
    glp_add_rows(lp, constraints); // allelicity_total + 1 (haplotype freqs normalization) ...


    // initialize glpk ...
    // mind the 1-based numbering ...
    // just allele freqs ...
    for (int i = 1; i < pop->independent_alleles+1; i++){
        char buf[64];
        sprintf(buf,"F%d",i);
        glp_set_row_name(lp, i, buf);
        // glp_set_row_bnds(lp, i, GLP_FX, 0, 0);
    }
    // normalization constraint:
    glp_set_row_name(lp, pop->independent_alleles+1, "TOT");
    glp_set_row_bnds(lp, pop->independent_alleles+1, GLP_FX, 1.0, 1.0);
    // linkages ...
    int lnk_idx = 1;
    for (int i = pop->independent_alleles+2; i < constraints+1; i++){
        char buf[64];
        sprintf(buf,"LNK%d",lnk_idx);
        glp_set_row_name(lp, i, buf);
        lnk_idx++;
        // glp_set_row_bnds(lp, i, GLP_FX, 0, 0);
    }



    // add columns ...
    glp_add_cols(lp, pop->haplotypes); // number of haplotypes ...
    // initialize columns ...
    // mind 1-based numbering in glpk ...
    for (int h = 0; h < pop->haplotypes; h++){
        char buf[64];
        sprintf(buf,"H%d",h+1);
        glp_set_col_name(lp, h+1, buf);
        glp_set_col_bnds(lp, h+1, GLP_DB, 0.0, 1.0);
        // glp_set_obj_coef(lp, h+1, 0.0);
    }


    // now we'd need to make arrays in the SparseMatrix - 1-numbered and turn values into floating point numbers ...
    // and the arrays are already generated!
    glp_load_matrix(lp, mat->nonzero, mat->irow, mat->jcol, mat->mval);

    /////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    // THIS IS THE MAIN WORKING BLOC
    /////////////////////////////////////////////////////////////////////
    // THIS WAS INSIDE WHICH_HAP LOOP, BUT NOT ANYMORE - WORKS FINE SO FAR ...
    // set constraint values: frequency of alleles ...
    // ERROR WAS HERE, ALLELES WERE INSTEAD OF ALLELES+1 ! LAST FREQ WASN'T GETTING INTO THE SYSTEM ...
    for (int i = 1; i < pop->independent_alleles+1; i++){
        glp_set_row_bnds(lp, i, GLP_FX, pop->freqs[i-1], pop->freqs[i-1]);        
    }
    // normalization ...
    // all haplotypes sum up to 1.0 ... 
    glp_set_row_bnds(lp, pop->independent_alleles+1, GLP_FX, 1.0, 1.0);
    // linkages ...
    lnk_idx = 0; // zero based indexing for most of the internal arrays ...
    for (int i = pop->independent_alleles+2; i < constraints+1; i++){
        glp_set_row_bnds(lp, i, GLP_FX, pop->link_freq[lnk_idx], pop->link_freq[lnk_idx]);
        lnk_idx++;      
    }
    /////////////////////////////////////
    /////////////////////////////////////
    /////////////////////////////////////
    // allocate arrays for the haplotype frequency bounds ....
    pop->hbounds = (HBound *)malloc(pop->haplotypes*sizeof(HBound));
    hbounds_prescan = (HBound *)malloc(pop->haplotypes*sizeof(HBound));
    hbounds_prescan_sorted = (HBound *)malloc(pop->haplotypes*sizeof(HBound));

    // run simplex for all of the haplotypes, MIN and MAX ...
    run_simplex_all_haplos(pop, lp, &parm);

    // (UNSORTED) maybe one should be storing pop->hbounds array before proceeding ...
    memcpy(hbounds_prescan, pop->hbounds, pop->haplotypes*sizeof(HBound));

    // We'd need to determine most interesting haplotypes (haplotypes to scan) ourselves by sorting arrays jointly ...
    // then we'd need to update our linkage info, adding fixed constraint for haplotype under scanning, and perform simplex again and again ...
    // so, after the initial simplex, sort haplotypes by the upper limit and ...
    qsort(pop->hbounds, pop->haplotypes, sizeof(HBound), compare_hbounds);

    // (SORTED) maybe one should be storing pop->hbounds array before proceeding ...
    memcpy(hbounds_prescan_sorted, pop->hbounds, pop->haplotypes*sizeof(HBound));

    // AREA USED TO BE UNDER CONSTRUCTION ...
    int num_haplo_scan = 21;
    int num_freq_steps = 12;

    printf("HAPLOTYPES & BOUNDS:\n\n");
    print_haplotypes_and_bounds(pop, 1, 1, num_haplo_scan, stdout, NULL);
    print_haplotypes_and_bounds(pop, 1, 1, ALL_HAPLO, NULL, "original_out.dat");

    // We'd need to determine most interesting haplotypes (haplotypes to scan) ourselves by sorting arrays jointly ...
    // then we'd need to update our linkage info, adding fixed constraint for haplotype under scanning, and perform simplex again and again ...

    // PREPARING FOR THE SCANNING VERSION ...
    // 
    // DO THE OPERATION THAT HAS TO BE DONE JUST ONCE:

    // No additional constraints are needed, just apply boundary settings for the auxiliary variable corresponding to the haplotype of interest

    // FINAL OUTPUTTING
    FILE *fp;
    fp = fopen("violin.dat","w");
    if (fp == NULL) {
        fprintf(stderr,"Cannot open violin.dat for writing ...");
        exit(1);
    }
    fprintf(fp,"haplo_idx haplo fmin fmax scan_idx fscan ratio_prod\n");



    // next step would be to go through a ~dozen of "top" hbounds and fix corresponding haplotype frequencies, by adding extra linkage ...
    for (int i = 0; i < num_haplo_scan; i++)
    {
        // for each of the haplos do the scanning ...
        int haplo_idx = hbounds_prescan_sorted[i].index + 1;
        double freq_min = hbounds_prescan_sorted[i].lower;
        double freq_max = hbounds_prescan_sorted[i].upper;
        double freq_range = freq_max - freq_min;
        double freq_step = freq_range/num_freq_steps;
        // if (freq_range >= 0.2) { freq_step = 0.02; }
        // else if ((freq_range>=0.1)&&(freq_range<0.2)) { freq_step = 0.01; }
        // else if ((freq_range>=0.05)&&(freq_range<0.1)) { freq_step = 0.005; }
        // else { freq_step = 0.001; }

        char *haplo_descr;
        haplo_descr = (char *)malloc(20*pop->loci*sizeof(char));
        if (haplo_descr == NULL){fprintf(stderr,"Cannot allocate memory for haplo_descr\n"); exit(1);}
        ////////////////////////////////////////////////////////
        get_haplotype_description(haplo_descr, pop, haplo_idx-1);
        // print some info ...
        printf("\nPreparing to scan for haplo %d in a range from %.4lf to %.4lf, using %.4lf step ...\n", haplo_idx, freq_min, freq_max, freq_step);
        printf("haplotype %d description: %s\n", haplo_idx, haplo_descr);


        // the frequency we would use for fixation ...
        double freq_scan = freq_min;
        int dummy_counter = 1;
        ////////////////////////////////////////////// 
        while (freq_scan <= freq_max){
            // Simply fixing the auxilary variable works great!
            glp_set_col_bnds(lp, haplo_idx, GLP_FX, freq_scan, freq_scan);

            // the constraint matrix stays the same - no need to reassign and/or reload it!
            // RUN SIMPLEX FOR THE MODIFIED PROBLEM ...
            // run simplex for all of the haplotypes, MIN and MAX ...
            run_simplex_all_haplos(pop, lp, &parm);

            // // we need to calculate some metric on how much did this fixation affected the boundaries ..
            // //
            // // (1) variant: the product of all haplotype boundary ratios ...
            // // excluding the one we've fixed of course
            // double ratio_prod = get_ratio_product(pop, hbounds_prescan, haplo_idx, ALL_HAPLO);
            // //
            // // (2) variant: the product of only `num_haplo_scan` top haplo-ratios ...
            // // excluding the one we've fixed of course
            double ratio_prod = get_ratio_product(pop, hbounds_prescan_sorted, haplo_idx, num_haplo_scan);

            // // print current freq_scan ...
            // printf("Currently haplotype %d is fixed at %.4lf frequency\n", haplo_idx, freq_scan);
            printf("%d %.4lf %.4lf\n",dummy_counter,freq_scan,ratio_prod);
            fprintf(fp,"%d %s %.4lf %.4lf %d %.4lf %.4lf\n", haplo_idx, haplo_descr, freq_min, freq_max, dummy_counter, freq_scan, ratio_prod);


            // NOTE: in the master-script, those `ratio_products` were used as violin-plot widths ('Y'-coord)
            // along with the corresponding `freq_scan` as an 'X'-coord
            // That means, that we'd have to store this values somehow for each fixed `freq_scan` and each scanned haplo,
            // and then use this data to plot violin-plots using python or whatever ...

            // TO BE CONTINUED ...

            // // // DO WE NEED THAT???
            // // // so, after the initial simplex, sort haplotypes by the upper limit and ...
            // // qsort(pop->hbounds, pop->haplotypes, sizeof(HBound), compare_hbounds);
            // printf("HAPLOTYPES & BOUNDS:\n\n");
            // print_haplotypes_and_bounds(pop, 1, 1, 20, stdout, NULL);

            // incremet freq_scan until we reach freq_max ...
            freq_scan += freq_step;
            dummy_counter++;
            // that's a fix for 'freq_range' that is 0-width ...
            if (freq_step == 0.0) {
                break;
            }
        }

        // WE MUST RELEASE THAT FIXED HAPLOTYPE ...
        glp_set_col_bnds(lp, haplo_idx, GLP_DB, 0.0, 1.0);
        // NOW WE CAN MOVE ON TO SCANNING NEXT HAPLOTYPE ...
        printf("\n");
        // free it right away ...
        free(haplo_descr);


       
    }

    fclose(fp);


    glp_delete_prob(lp);

    // deallocate population structure description ....
    free_population(pop);
    free(pop);

    // deallocate sparse matrix description of the constraint matrix ...
    free_sparse_matrix(mat);
    free(mat);

    // success ...
    return 0;


}





// functions

int get_partial_product(int *array, int start, int stop){
    int par_product = 1;
    for (int i = start+1; i < stop; i++){
        par_product *= array[i];
    }
    return par_product;
}


void load_allele_frequency_file_safe(PopDes *population, const char *fname){
    // this function read input file and tries to interpret the contents.
    // input file format is intended to be strict to make the code cleaner.
    //
    // variables needed
    FILE *fp;
    int loci_number, haplotypes_number;
    int allelicity_total;
    int allelicity_independent;
    int *allelicity_numbers;
    double *allele_freqs;
    char buffer[BUFFER_SIZE];
    char **labels;
    //
    // opening file carefully ...
    fp = fopen(fname,"r");
    if (fp == NULL) { fprintf(stderr, "failed to open file: %s\n", fname); exit(1); }
    //
    printf("Reading input file: %s ...\n",fname);
    //
    printf("(1) Extracting # of loci and # of alleles at each loci. \n");
    // first line is LOCI #_of_loci
    if ( fgets(buffer, BUFFER_SIZE, fp) != NULL){
        if ( sscanf(buffer,"LOCI %d\n",&loci_number) != 1 ) {
            fprintf(stderr, "failed to parse: %s\n",buffer);
            exit(1);
        }
    } else {
        fprintf(stderr, "failed to read file: %s\n",fname);
        exit(1);        
    }
    //
    // iterate loci_number times to extract allele numbers ...
    allelicity_total = 0;
    allelicity_independent = 0;
    haplotypes_number = 1;
    allelicity_numbers = (int *)malloc(loci_number*sizeof(int));
    if (allelicity_numbers == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    for (int locus = 0; locus < loci_number; locus++){
        // reading number of alleles for each locus ...
        if ( fgets(buffer, BUFFER_SIZE, fp) != NULL ) {
            int index;
            if ( sscanf(buffer, "L%d %d\n", &index, &allelicity_numbers[locus]) != 2 ){
                fprintf(stderr, "failed to parse: %s\n",buffer);
                exit(1);
            }
            // check loci numbering if parsing is successful ...
            if ( index != locus+1) {
                fprintf(stderr, "loci numbering problem: %s\n", buffer);
                exit(1);
            }
        // if fgets failed ...
        } else {
            fprintf(stderr, "failed to read file: %s\n",fname);
            exit(1);        
        }
        //
        // calculate allelicity_total and haplotype_number ...
        allelicity_total += allelicity_numbers[locus]; // sum of all alleles ...
        allelicity_independent += (allelicity_numbers[locus]-1); // independent alleles only ...
        haplotypes_number *= allelicity_numbers[locus]; // prod of all alleles - # of haplotypes ...
    } // locus-loop is over
    // 
    // allocate memory for labels ...
    labels = (char **)malloc(allelicity_total*sizeof(char *));
    if (labels == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    for (int label_idx = 0; label_idx < allelicity_total; label_idx++) {
        labels[label_idx] = (char *)malloc(LABEL_SIZE*sizeof(char));
        if (labels[label_idx] == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    }
    // extract all the frequencies (only linearly independent frequencies!)
    allele_freqs = (double *)malloc(allelicity_independent*sizeof(double));
    if (allele_freqs == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    int label_through_index = 0;   
    int allele_through_index = 0;
    // locus by locus ...
    //
    printf("(2) Extracting alleles' frequencies and labels(optional) for %d loci\n", loci_number);
    //
    for (int locus = 0; locus < loci_number; locus++){
        // scan through all but the last allele, as it must be AUTO ...
        for (int a = 0; a < allelicity_numbers[locus]-1; a++){
            // extract actual frequencies ...
            //
            //
            if ( fgets(buffer, BUFFER_SIZE, fp) != NULL){
                int l_index; // locus index for verification
                int a_index; // allele index for verification
                // try extracting freq and label ...
                if (sscanf(buffer,"L%d A%d %lf %s\n", &l_index, &a_index, &allele_freqs[allele_through_index], labels[label_through_index])!=4) {
                    // try freq alone ...
                    if (sscanf(buffer,"L%d A%d %lf\n", &l_index, &a_index, &allele_freqs[allele_through_index])!=3){
                        fprintf(stderr, "failed to parse: %s\n",buffer);
                        exit(1);
                    } else {
                        // store L_iA_j in the labels[label_through_index] if none is provided ...
                        sprintf(labels[label_through_index],"L%dA%d",l_index,a_index);
                    }
                }
                // check locus and allele indexing if parsing succeeded 
                if ((l_index != locus+1)||(a_index != a+1)) {
                    fprintf(stderr, "loci/allele numbering problem: %s\n", buffer);
                    exit(1);
                }
            } else {
                fprintf(stderr, "failed to read file: %s\n",fname);
                exit(1);        
            }
            // incrementing through index
            allele_through_index++;
            label_through_index++;
        } // a-loop is over
        // skipping AUTO, which must be the last element ...
        if ( fgets(buffer, BUFFER_SIZE, fp) != NULL){
            int l_index; // locus index for verification
            int a_index; // allele index for verification
            char auto_check[32]; // array to read the word 'AUTO'
            // try AUTO and a label ...
            if (sscanf(buffer,"L%d A%d %s %s\n", &l_index, &a_index, auto_check, labels[label_through_index])!=4) {
                if (sscanf(buffer,"L%d A%d %s\n", &l_index, &a_index, auto_check)!=3){
                    fprintf(stderr, "failed to parse: %s\n",buffer);
                    exit(1);
                } else {
                    // store L_iA_j in the labels[label_through_index] if none is provided ...
                    sprintf(labels[label_through_index],"L%dA%d",l_index,a_index);
                }
            }
            // check AUTO and loci/allele indexing if parsing succeeded 
            if ((l_index != locus+1)||(a_index != allelicity_numbers[locus])||(strcmp(auto_check,"AUTO")!=0)) {
                fprintf(stderr, "loci/allele/AUTO problem: %s\n", buffer);
                exit(1);
            }
        } else {
            fprintf(stderr, "failed to read file: %s\n",fname);
            exit(1);        
        }
        //  increment label through index only as there is no freq associated with AUTO ...
        label_through_index++;
    } // locus-loop is over 
    // assertions for verification ..
    assert(label_through_index == allelicity_total);
    assert(allele_through_index == allelicity_independent);
    //
    /////////////////////////////////////////////////////////
    ////// LINKAGE STUFF ...
    /////////////////////////////////////////////////////////
    int linkage_number;
    int *link_alleles_number;
    double *link_freq;
    int link_array_length = 0;
    int *link_loci;
    int *link_alleles;
    //
    printf("(3) Extracting linkage information(optional)\n");
    //
    if ( fgets(buffer, BUFFER_SIZE, fp) != NULL){
        if ( sscanf(buffer,"LINKAGE %d\n",&linkage_number) != 1 ) {
            fprintf(stderr, "failed to parse: %s\n",buffer);
            exit(1);
        }
    } else {
        // it simply means that there is no linkage information ...
        linkage_number = 0;
    }
    // allocations ...
    link_alleles_number = (int *)malloc(linkage_number*sizeof(int));
    if (link_alleles_number == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    link_freq = (double *)malloc(linkage_number*sizeof(double));
    if (link_freq == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    // read info about links ...
    for (int link = 0; link < linkage_number; link++){
        //
        if ( fgets(buffer, BUFFER_SIZE, fp) != NULL){
            int lnk_index; // link index for verification ...
            if ( sscanf(buffer,"LINK%d %d %lf\n", &lnk_index, &link_alleles_number[link], &link_freq[link]) != 3 ) {
                fprintf(stderr, "failed to parse: %s\n",buffer);
                exit(1);
            }
            // check link numbering persistence ... 
            if (lnk_index != link+1) {
                fprintf(stderr, "link numbering problem: %s\n", buffer);
                exit(1);
            }
        } else {
            fprintf(stderr, "failed to read file: %s\n",fname);
            exit(1);        
        }
        // count total number of alleles(@ all loci) involved in all linkages ...
        link_array_length += link_alleles_number[link];
    } // link-loop is over ...
    // more allocations ...
    link_loci = (int *)malloc(link_array_length*sizeof(int));
    if (link_loci == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    link_alleles = (int *)malloc(link_array_length*sizeof(int));
    if (link_alleles == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    // collect involved loci and alleles ...
    // for each link ...
    int link_al_through = 0;
    for (int link = 0; link < linkage_number; link++) {
        // parse all involved loci ...
        for (int involved_locus = 0; involved_locus < link_alleles_number[link]; involved_locus++) {
            int locus_idx,allele_idx;
            if ( fgets(buffer, BUFFER_SIZE, fp) != NULL){
                int lnk_index; // link index for verification ...
                if ( sscanf(buffer,"LINK%d L%d A%d\n", &lnk_index, &locus_idx, &allele_idx) != 3 ) {
                    fprintf(stderr, "failed to parse: %s\n",buffer);
                    exit(1);
                }
                // check link numbering persistence ... 
                if (lnk_index != link+1) {
                    fprintf(stderr, "link numbering problem: %s\n", buffer);
                    exit(1);
                }
                //
            } else {
                fprintf(stderr, "failed to read file: %s\n",fname);
                exit(1);        
            }
            // storing the parsed loci and alleles indices:
            // we're using zero-numbering internally, thus ...
            link_loci[link_al_through] = locus_idx-1;
            link_alleles[link_al_through] = allele_idx-1;
            link_al_through++;
        }
    } // link-loop is over ...
    // verification 
    assert(link_array_length == link_al_through);
    //
    printf("(4) %d linkages processed\n", linkage_number);
    /////////////////////////////////////////////////////////
    ////// LINKAGE STUFF ...
    /////////////////////////////////////////////////////////
    //
    // now fill in the population structure ...
    // population...
    population->loci = loci_number;
    population->haplotypes = haplotypes_number;
    population->alleles_total = allelicity_total;
    population->independent_alleles = allelicity_independent;
    // allocate population arrays ...
    population->allelicity = (int *)malloc(population->loci*sizeof(int));
    if (population->allelicity == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    population->freqs = (double *)malloc(population->independent_alleles*sizeof(double)); // independent freqs only ...
    if (population->freqs == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    // labels allocation ...
    population->allele_labels = (char **)malloc(allelicity_total*sizeof(char *));
    if (population->allele_labels == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    for (int label_idx = 0; label_idx < allelicity_total; label_idx++) {
        population->allele_labels[label_idx] = (char *)malloc(LABEL_SIZE*sizeof(char));
        if (population->allele_labels[label_idx] == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    }
    //
    // copy local arrays to the population arrays ...
    memcpy(population->allelicity,allelicity_numbers,(size_t)population->loci*sizeof(int));
    memcpy(population->freqs,allele_freqs,(size_t)population->independent_alleles*sizeof(double)); // independent freqs only
    // copy labels ...
    for (int label_idx = 0; label_idx < allelicity_total; label_idx++) {
        strcpy(population->allele_labels[label_idx], labels[label_idx]);
    }
    //
    // linkages ...
    population->linkages = linkage_number;
    population->link_array_length = link_array_length;
    population->link_alleles_number = (int *)malloc(linkage_number*sizeof(int));
    if (population->link_alleles_number == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    population->link_freq = (double *)malloc(linkage_number*sizeof(double));
    if (population->link_freq == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    population->link_loci_indices = (int *)malloc(link_array_length*sizeof(int));
    if (population->link_loci_indices == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    population->link_alleles_indices = (int *)malloc(link_array_length*sizeof(int));
    if (population->link_alleles_indices == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    // linkages info copying ...
    memcpy(population->link_alleles_number,link_alleles_number,(size_t)linkage_number*sizeof(int));
    memcpy(population->link_freq,link_freq,(size_t)linkage_number*sizeof(double));
    memcpy(population->link_loci_indices,link_loci,(size_t)link_array_length*sizeof(int));
    memcpy(population->link_alleles_indices,link_alleles,(size_t)link_array_length*sizeof(int));
    //
    // freeing allocated memory ...
    free(allelicity_numbers);
    free(allele_freqs);
    for (int label_idx = 0; label_idx < allelicity_total; label_idx++) {
        free(labels[label_idx]);
    }
    free(labels);
    free(link_alleles_number);
    free(link_freq);
    free(link_loci);
    free(link_alleles);
    //
    // closing input file
    fclose(fp);
    //
    return;
}



void generate_haplotypes(PopDes *population){
    //
    // allocate a 2D array to enumerate all haplotypes ...
    population->haplotype_description = (int **)malloc(population->haplotypes*sizeof(int *));
    if (population->haplotype_description == NULL) {
        fprintf(stderr, "unable to allocate memory: %s\n", "generate_haplotypes");
        exit(1);
    }
    for (int i = 0; i < population->haplotypes; i++){
        population->haplotype_description[i] = (int *)malloc(population->loci*sizeof(int));
        if (population->haplotype_description[i] == NULL) {
            fprintf(stderr, "unable to allocate memory: %s\n", "generate_haplotypes");
            exit(1);
        }
    }
    // let's fill in the array ...
    for (int l = 0; l < population->loci; l++){
        // get the stride ...
        int stride = get_partial_product(population->allelicity, l, population->loci);
        for (int h = 0; h < population->haplotypes; h++){
            population->haplotype_description[h][l] = (h/stride)%population->allelicity[l];
        }
    }
    //
    return;
}


void generate_linked_loci_number(PopDes *population){
    //
    // actual number of loci involved per linkage to be allocated here ...
    population->link_loci_number = (int *)malloc(population->linkages*sizeof(int));
    if (population->link_loci_number == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    // temp array to store number of times each loci is involved in a particular linkage ...
    int *tmp_loci_times_involved = (int *)malloc(population->loci*sizeof(int));
    if (tmp_loci_times_involved == NULL) {fprintf(stderr,"unable to allocate memory\n"); exit(1); }
    //
    int link_item_through = 0;
    for (int link = 0; link < population->linkages; link++){
        //
        // for each link ...
        // purge loci_times_involved array:
        for (int locus = 0; locus < population->loci; locus++) { tmp_loci_times_involved[locus] = 0; }
        // track loci involved within those alleles involved in that particular linkage:
        for (int i = 0; i < population->link_alleles_number[link]; i++){
            int locus_involved = population->link_loci_indices[link_item_through + i];
            tmp_loci_times_involved[locus_involved]++;
        }
        //
        // now count loci involved (non-zero entries in tmp_loci_times_involved):
        int actual_loci_involved = 0;
        for (int locus = 0; locus < population->loci; locus++) {
            if (tmp_loci_times_involved[locus] >= 1) actual_loci_involved++;
        }
        // done counting ...
        // store that number in population->link_loci_number:
        population->link_loci_number[link] = actual_loci_involved;
        // after all shift the through-index of the involved loci&alleles ...
        link_item_through += population->link_alleles_number[link];
    }
    // release that temporary array 
    free(tmp_loci_times_involved);
}



void generate_constraint_matrix_linkage(PopDes *population, SparseMat *matrix){
    //
    // generate linked loci number for each link:
    // printf("before generate_linked_loci_number\n");
    generate_linked_loci_number(population);
    // printf("after generate_linked_loci_number\n");
    // this will populate population->link_loci_number array
    // which is uninitialized otherwise.
    // ... population->link_loci_number[link] ...
    //
    // calculate non-zero elements in the sparse matrix string the constraints ...
    int length = 0;
    // first, just the allele frequencies constraints:
    for (int l = 0; l < population->loci; l++){
        // number of independent equation for the locus BY ...
        // number of combination with a given allele at this locus:
        length += (population->allelicity[l]-1) * (population->haplotypes/(population->allelicity[l]));
    }
    // normalization condition, sum of all haplotype fractions is 1:
    length += population->haplotypes;
    //
    // linkage related items to be added ...
    // // to add linkage-related numbers we'd just need the number of alleles involved in each link ...
    // // through-index of the loci and alleles involved in the linkages ...
    // int link_item_through = 0;
    // for (int link = 0; link < population->linkages; link++){
    //     // haps involved is a product of allelicities of all of the loci except for those in the linkage ...
    //     int haps_involved = population->haplotypes;
    //     for (int link_item = 0; link_item < population->link_alleles_number[link]; link_item++){
    //         // extract number of alleles possible at that particular locus of interest (linked locus)...
    //         int involved_locus_allelicity = population->allelicity[population->link_loci_indices[link_item_through]];
    //         haps_involved /= involved_locus_allelicity;
    //         link_item_through++;
    //     }
    //     length += haps_involved;
    // }
    // //
    // simple enumeration of haplotypes involved in the linkage information ...
    // through-index of the loci and alleles involved in the linkages ...
    int link_item_through = 0;
    for (int link = 0; link < population->linkages; link++){
        // for each link ...
        // scan through haplotypes ...
        for (int h = 0; h < population->haplotypes; h++){
            //scan through haplotypes ...
            int matching_alleles = 0; // true ...
            for (int i = 0; i < population->link_alleles_number[link]; i++){
                int locus_involved = population->link_loci_indices[link_item_through + i];
                int allele_involved = population->link_alleles_indices[link_item_through + i];
                if (population->haplotype_description[h][locus_involved]==allele_involved){
                    matching_alleles ++;
                }
            }
            // check if that haplotype is involved ...
            // matching alleles must be equal to the number of involved loci (!) which could be smaller than the number of involved alleles... 
            if ( matching_alleles == population->link_loci_number[link] ){ length++; }
        }
        // after all haplotypes are scanned, shift the through-index of the involved loci&alleles ...
        link_item_through += population->link_alleles_number[link];
    }
    //
    // just a simple sanity check ...
    assert(population->link_array_length == link_item_through);
    //
    // allocating Sparse Matrix related array given the length:
    matrix->nonzero = length;
    matrix->rows = population->independent_alleles + 1 + population->linkages; // allele freqs (linear independence accounted) + normalization + linkages ...
    matrix->cols = population->haplotypes; // presumably all haplotypes are involved, it's the case 
    matrix->irow = (int *)malloc((length+1)*sizeof(int));
    if (matrix->irow == NULL) {
        fprintf(stderr, "unable to allocate memory: %s\n","generate_constraint_matrix_linkage");
        exit(1);
    }
    matrix->jcol = (int *)malloc((length+1)*sizeof(int));
    if (matrix->jcol == NULL) {
        fprintf(stderr, "unable to allocate memory: %s\n","generate_constraint_matrix_linkage");
        exit(1);
    }
    matrix->mval = (double *)malloc((length+1)*sizeof(double));
    if (matrix->mval == NULL) {
        fprintf(stderr, "unable to allocate memory: %s\n","generate_constraint_matrix_linkage");
        exit(1);
    }
    //
    // we'd also need two more running indicies here ///
    int row = 0; // row index of the constraint matrix ...
    int idx = 1; // index along the sparse matrix arrays ...
    for (int l = 0; l < population->loci; l++){
        // allele frequencies obey normalization, so use just allelicity-1 of them as a constraint ...
        for (int a = 0; a < population->allelicity[l]-1; a++){
            // scan the entire set of haplotypes ...
            for (int h = 0; h < population->haplotypes; h++){            
                // checking if haplotypes match ...
                if (population->haplotype_description[h][l] == a){
                    // indexes must updated for the glpk - Fortran-style indexing ...
                    matrix->irow[idx] = row+1; // that's where we need that row index of that actual constraint matrix ...
                    matrix->jcol[idx] = h+1; // column index is simply the haplotype number ...
                    matrix->mval[idx] = 1.0;
                    idx++; // count every nonzero element ...
                }
            }
            row++; // constraints matrix row index incremented every other accounted allele ... 
        }
    }
    //
    // now, should add the normalization condition here ...
    for (int h = 0; h < population->haplotypes; h++){
        // indexes must updated for the glpk - Fortran-style indexing ...
        matrix->irow[idx] = row+1; // that's where we need that row index of that actual constraint matrix ...
        matrix->jcol[idx] = h+1; // column index is simply the haplotype number ...
        matrix->mval[idx] = 1.0;
        idx++; // count every nonzero element ...
    }
    row++;
    //
    // through-index of the loci and alleles involved in the linkages ...
    link_item_through = 0;
    for (int link = 0; link < population->linkages; link++){
        // for each link ...
        // scan through haplotypes ...
        for (int h = 0; h < population->haplotypes; h++){
            //scan through haplotypes ...
            int matching_alleles = 0; // true ...
            for (int i = 0; i < population->link_alleles_number[link]; i++){
                int locus_involved = population->link_loci_indices[link_item_through + i];
                int allele_involved = population->link_alleles_indices[link_item_through + i];
                if (population->haplotype_description[h][locus_involved]==allele_involved){
                    matching_alleles ++;
                }
            }
            //
            //
            // check if that haplotype is involved ...
            // matching alleles must be equal to the number of involved loci (!) which could be smaller than the number of involved alleles... 
            if ( matching_alleles == population->link_loci_number[link] ){
            // the old version, where >=2 alleles at a given haplotype were prohibited.
            // if (matching_alleles == population->link_alleles_number[link]){
                // indexes must be updated for the glpk - Fortran-style indexing ...
                matrix->irow[idx] = row+1; // that's where we need that row index of that actual constraint matrix ...
                matrix->jcol[idx] = h+1; // column index is simply the haplotype number ...
                matrix->mval[idx] = 1.0;
                idx++; // count every nonzero element ...                
            }
        }
        // after all haplotypes are scanned, shift the through-index of the involved loci&alleles ...
        link_item_through += population->link_alleles_number[link];
        row++;
    }
    //
    // folowing assertions must be met:
    assert(matrix->rows == row);
    assert(matrix->nonzero+1 == idx);
    assert(population->link_array_length == link_item_through);
    //
    //
    return;
}


void free_population(PopDes *population){
    //
    free(population->allelicity);
    free(population->freqs);
    // free labels ...
    for (int i = 0; i < population->alleles_total; i++) {
        free(population->allele_labels[i]);
    }
    free(population->allele_labels);
    // haplotypes description freeing ...
    for (int i = 0; i < population->haplotypes; i++){
        free(population->haplotype_description[i]);
    }
    free(population->haplotype_description);
    // linkages ...
    free(population->link_alleles_number);
    // link_loci_number is not necessarily generated, so check it first:
    if (population->link_loci_number != NULL) free(population->link_loci_number);
    free(population->link_freq);
    free(population->link_loci_indices);
    free(population->link_alleles_indices);
    //results ...
    free(population->hbounds);
    // free(population->upper);
}


void free_sparse_matrix(SparseMat *matrix){
    // simple 1D arrays deallocation ...
    free(matrix->irow);
    free(matrix->jcol);
    free(matrix->mval);
}



void print_haplotypes_and_bounds(PopDes *population, int print_labels, int print_bounds, int num_haplo_print, FILE *fp, const char *fname){
    // where to output ...
    if (fp != NULL) {
        /* code */
        if (fileno(fp) != STDOUT_FILENO) {
            fprintf(stderr, "This is meant for terminal output: stdout\n");
            exit(1);
        }
        printf("terminal output ...\n");
    } else if (fname != NULL) {
        /* code */
        fp = fopen(fname,"w");
        if (fp == NULL) { fprintf(stderr, "failed to open file: %s\n", fname); exit(1); }
    } else {
        // no output destination provided ...
        fprintf(stderr, "No output destination provided for 'print_haplotypes_and_bounds'\n");
        exit(1);
    }
    // 
    // 
    // print every haplotype ...
    num_haplo_print = (num_haplo_print != ALL_HAPLO) ? num_haplo_print : population->haplotypes;
    // 
    for (int h = 0; h < num_haplo_print; h++){
        if (print_bounds) {
            fprintf(fp,"H%04d: %05d %.4lf %.4lf ",h+1,population->hbounds[h].index,population->hbounds[h].lower,population->hbounds[h].upper);
        } else {
            fprintf(fp,"H%04d: ",h+1);
        }
        // labels are indexed with a through index, so use the shift for different loci ...
        int haplo_hdx = population->hbounds[h].index;
        int allele_through_index_shift = 0;
        for (int l = 0; l < population->loci; l++){
            // get allele index at this locus:
            int allele_index_here = population->haplotype_description[haplo_hdx][l];
            if (print_labels) {
                // print label is specified ...
                fprintf(fp,"%s:",population->allele_labels[allele_through_index_shift + allele_index_here]);
            } else{
                // otherwise print index(1-based numbering)
                fprintf(fp,"%d ",population->haplotype_description[haplo_hdx][l]+1);
            }
            // increment the shift using allelicity:
            allele_through_index_shift += population->allelicity[l];
        }
        // next haplotype:
        fprintf(fp,"\n");
    }
    // 
    // 
    // closing FILE if it was a user created one ONLY ...
    if (fileno(fp) != STDOUT_FILENO) {
        fclose(fp);
    }

}





//////////////////////////////////////////////// wrap this into a separate standing function ...
void run_simplex_all_haplos(PopDes *population, glp_prob *local_lp, glp_smcp *local_parm){
    for (int which_hap = 0; which_hap < population->haplotypes; which_hap++){
        // for each haplotype calculate its bounds ...
        //
        //
        // set the objective function, i.e. which haplotype do you want to find limits for?
        for (int i = 1; i < population->haplotypes+1; i++){
            glp_set_obj_coef(local_lp, i, 0.0);
        }
        glp_set_obj_coef(local_lp, which_hap+1, 1.0); // haplotype number which_hp+1 !
        //
        // uncomment this line if you'd like to see how GLPK formulate optimization problem:
        // glp_write_mps(local_lp, GLP_MPS_FILE, NULL, "zzN.mps");
        //
        glp_set_obj_dir(local_lp, GLP_MAX); // maximize
        glp_simplex(local_lp, local_parm);
        double upper_bound = glp_get_obj_val(local_lp);
        //
        glp_set_obj_dir(local_lp, GLP_MIN); //minimize ...
        glp_simplex(local_lp, local_parm);
        double lower_bound = glp_get_obj_val(local_lp);
        // fill in the lower and upper arrays of the population structure ...
        population->hbounds[which_hap].index = which_hap;
        population->hbounds[which_hap].lower = lower_bound;
        population->hbounds[which_hap].upper = upper_bound;
    }

}




void print_constraints_matrix(SparseMat *matrix){
    //
    // allocate 2D constraint matrix, just for printing ...
    int **check_mat;
    check_mat = (int **)malloc(matrix->rows*sizeof(int *));
    if (check_mat == NULL) {
        fprintf(stderr, "unable to allocate memory: %s\n", "print_constraints_matrix");
        exit(1);
    }
    for (int i = 0; i < matrix->rows; i++){
        check_mat[i] = (int *)malloc(matrix->cols*sizeof(int));
        if (check_mat[i] == NULL) {
            fprintf(stderr, "unable to allocate memory: %s\n", "print_constraints_matrix");
            exit(1);
        }
    }
    // fill it in with zeroes ...
    for (int i = 0; i < matrix->rows; i++){
        for (int j = 0; j < matrix->cols; j++){
            check_mat[i][j] = 0;
        }
    }
    // fill in the nonzero constraint elements ...
    for (int i = 0; i < matrix->nonzero; i++){
        // mind the 1-based indexing in the matrix ...
        check_mat[matrix->irow[i+1]-1][matrix->jcol[i+1]-1] = matrix->mval[i+1];
    }
    printf("\nPrint the matrix of constraints:\n");
    // printing that matrix out ...
    for (int i = 0; i < matrix->rows; i++){
        for (int j = 0; j < matrix->cols; j++){
            printf("%d ",check_mat[i][j]);
        }
        printf("\n");
    }
    //
    // free allocated memory:
    for (int i = 0; i < matrix->rows; i++){
        free(check_mat[i]);
    }
    free(check_mat);
}





int compare_hbounds(const void *ha, const void *hb){
    // inspired by http://www.anyexample.com/programming/c/qsort__sorting_array_of_strings__integers_and_structs.xml
    const HBound *iha = (const HBound *)ha;
    const HBound *ihb = (const HBound *)hb;
    // descending, sorted by HBound.upper ...
    return (int) 100000.0f*(ihb->upper - iha->upper);
    // // ascending, sorted by HBound.upper ...
    // return (int) 100000.0f*(iha->upper - ihb->upper);
    // // dscending, sorted by HBound.lower ...
    // return (int) 100000.0f*(ihb->lower - iha->lower);
    // /////////////////////////////////////////////////
    // /* integer comparison: returns negative if b > a 
    // and positive if a > b  
};




double get_ratio_product(PopDes *population, HBound *hbounds_original, int fixed_haplo_index, int product_number){
    double ratio_product = 1.0;
    double updated_range, original_range, range_ratio;
    int haplo_jdx, num_haplo_scan;
    // there are two regimes, ALL or a given number of TOP haplos ...
    if (product_number == ALL_HAPLO) {
        // ALL HAPLOS ...
        num_haplo_scan = population->haplotypes;
        for (int j = 0; j < num_haplo_scan; j++) {
            // get j-top haplotype index in the nonsorted array or just ALL of them ...
            haplo_jdx = j;
            // excluding fixed haplotype ...
            if (haplo_jdx != fixed_haplo_index-1) {
                updated_range  = (population->hbounds[haplo_jdx].upper - population->hbounds[haplo_jdx].lower);
                // product for ALL (haplo_jdx=j) and hbounds_original is unsorted, 
                // otherwise (haplo_jdx=hbounds_original[j].index), where hbounds_original is sorted(!)
                original_range = (hbounds_original[haplo_jdx].upper - hbounds_original[haplo_jdx].lower);
                range_ratio = (original_range > 0.0) ? (updated_range/original_range) : 1.0;
                // multiply only if 'range_ratio' is defined ...  
                ratio_product *= range_ratio;
            }
        }
    } else {
        // FIXED NUMBER ...
        num_haplo_scan = product_number;
        for (int j = 0; j < num_haplo_scan; j++) {
            // get j-top haplotype index in the nonsorted array or just ALL of them ...
            haplo_jdx = hbounds_original[j].index;
            // excluding fixed haplotype ...
            if (haplo_jdx != fixed_haplo_index-1) {
                updated_range  = (population->hbounds[haplo_jdx].upper - population->hbounds[haplo_jdx].lower);
                // product for ALL (haplo_jdx=j) and hbounds_original is unsorted, 
                // otherwise (haplo_jdx=hbounds_original[j].index), where hbounds_original is sorted(!)
                original_range = (hbounds_original[j].upper - hbounds_original[j].lower);
                range_ratio = (original_range > 0.0) ? (updated_range/original_range) : 1.0;
                // multiply only if 'range_ratio' is defined ...  
                ratio_product *= range_ratio;
            }
        }
    }
    // returnring ...
    return ratio_product;
}




void get_haplotype_description(char *haplo_descr, PopDes *population, int haplo_index){
    // int string_prefactor = 20; // allowing `string_prefactor` chars to describe a single allele at a single loci, just in case ...
    // // EXAMPLE DESCRIPTION:    L1A1:L2A1:L3A1:L4A3:
    // haplo_descr = (char *)malloc(string_prefactor*population->loci*sizeof(char));
    // if (haplo_descr == NULL){fprintf(stderr,"Cannot allocate memory for 'haplo_descr' in 'get_haplotype_description'\n"); exit(1);}
    int allele_through_index_shift = 0;
    int printed_len = 0;
    for (int l = 0; l < population->loci; l++){
        // get allele index at this locus:
        int allele_index_here = population->haplotype_description[haplo_index][l];
        // print label is specified ...
        printed_len += sprintf(haplo_descr+printed_len,"%s:",population->allele_labels[allele_through_index_shift + allele_index_here]);
        // printf("%s:",population->allele_labels[allele_through_index_shift + allele_index_here]);
        // increment the shift using allelicity:
        allele_through_index_shift += population->allelicity[l];
    }
    // printf("\n");
    return;
}








// void load_allele_frequency_file(PopDes *population, const char *fname){
//     // file format description ...
//     // The self-explanatory allele frequency description file example is to follow:
//     // LOCI 3
//     // L1 4
//     // L2 2
//     // L3 2
//     // L4 4
//     // L1 A1 0.1 Chr2_2repeat
//     // L1 A2 0.9 Chr2_3repeat
//     // L1 A3 AUTO Chr2_4repeat
//     // L2 A1 0.8 H274H
//     // L2 A2 AUTO H274Y
//     // L3 A1 0.01
//     // L3 A2 AUTO
//     // L4 A1 0.8 Seg1:B59
//     // L4 A2 0.07 Seg1:B59:A500T
//     // L4 A3 0.05 Seg1:B10
//     // L4 A4 AUTO Seg1:B10:G250C
//     // LINKAGE 2
//     // LINK 2 0.08
//     // LINK 3 0.07
//     // LINK1 L1 A2
//     // LINK1 L4 A4
//     // LINK2 L2 A2
//     // LINK2 L3 A2
//     // LINK2 L4 A2
//     // END
//     //
//     // The format of that file is STRICT! so, that the C-code remain readable and clean
//     // Allele label is the only flexible thing here, though ...
//     //
//     // variables needed
//     FILE *fp;
//     int loci_number, haplotypes_number;
//     int allelicity_total;
//     int allelicity_independent;
//     int *allelicity_numbers;
//     double *allele_freqs;
//     //
//     // open file if it's possible ...
//     fp = fopen(fname,"r");
//     //
//     // next line must the LOCI L, where L is the number of Loci; read line & get loci_number ...
//     fscanf(fp,"LOCI %d\n", &loci_number);
//     //
//     // iterate loci_number times to extract allele numbers ...
//     allelicity_total = 0;
//     allelicity_independent = 0;
//     haplotypes_number = 1;
//     allelicity_numbers = (int *)malloc(loci_number*sizeof(int));
//     for (int locus = 1; locus < loci_number+1; locus++){
//         // reading allele variants number ...
//         fscanf(fp,"L%*d %d\n", &allelicity_numbers[locus-1]);
//         // calculate allelicity_total and haplotype_number ...
//         allelicity_total += allelicity_numbers[locus-1]; // sum of all alleles ...
//         allelicity_independent += (allelicity_numbers[locus-1]-1); // independent alleles only ...
//         haplotypes_number *= allelicity_numbers[locus-1]; // prod of all alleles - # of haplotypes ...
//     }
//     //
//     // extract all the frequencies ...
//     int allele_idx = 0;
//     allele_freqs = (double *)malloc(allelicity_independent*sizeof(double));
//     // locus by locus ...
//     for (int i = 0; i < loci_number; i ++){
//         // scan through all but the last allele, as it must be AUTO ...
//         for (int j = 0; j < allelicity_numbers[i]-1; j++){
//             // extract actual frequencies ...
//             fscanf(fp,"L%*d A%*d %lf\n", &allele_freqs[allele_idx]);
//             allele_idx++;
//         }
//         // skipping AUTO, which must be the last element ...
//         char buf[1024];
//         fgets(buf,1024,fp); 
//         // fscanf(fp,"L%*d A%*d AUTO\n"); // doesn't work this way, line isn't read.
//     }
//     // that's it!
//     //
//     /////////////////////////////////////////////////////////
//     ////// LINKAGE STUFF ...
//     /////////////////////////////////////////////////////////
//     // Linkage goes here ...
//     int linkage_number;
//     int *link_alleles_number;
//     double *link_freq;
//     int link_array_length = 0;
//     int *link_loci;
//     int *link_alleles;
//     fscanf(fp,"LINKAGE %d\n",&linkage_number);
//     // allocations ...
//     link_alleles_number = (int *)malloc(linkage_number*sizeof(int));
//     link_freq = (double *)malloc(linkage_number*sizeof(double));
//     // read info about links ...
//     for (int link = 0; link < linkage_number; link++){
//         fscanf(fp,"LINK%*d %d %lf\n",&link_alleles_number[link],&link_freq[link]);
//         // count total number of alleles(@ all loci) involved in all linkages ...
//         link_array_length += link_alleles_number[link];
//     }
//     // more allocations ...
//     link_loci = (int *)malloc(link_array_length*sizeof(int));
//     link_alleles = (int *)malloc(link_array_length*sizeof(int));
//     // using the known total array length do:
//     for (int link_al = 0; link_al < link_array_length; link_al++){
//         int locus_idx,allele_idx;
//         fscanf(fp,"LINK%*d L%d A%d\n",&locus_idx,&allele_idx);
//         // we're using zero-numbering internally, thus ...
//         link_loci[link_al] = locus_idx-1;
//         link_alleles[link_al] = allele_idx-1;
//     }
//     /////////////////////////////////////////////////////////
//     ////// LINKAGE STUFF ...
//     /////////////////////////////////////////////////////////
//     //
//     // now fill in the population structure ...
//     // population...
//     population->loci = loci_number;
//     population->haplotypes = haplotypes_number;
//     population->alleles_total = allelicity_total;
//     population->independent_alleles = allelicity_independent;
//     // allocate population arrays ...
//     population->allelicity = (int *)malloc(population->loci*sizeof(int));
//     population->freqs = (double *)malloc(population->independent_alleles*sizeof(double));
//     //
//     // linkages ...
//     population->linkages = linkage_number;
//     population->link_alleles_number = (int *)malloc(linkage_number*sizeof(int));
//     population->link_freq = (double *)malloc(linkage_number*sizeof(double));
//     population->link_array_length = link_array_length;
//     population->link_loci_indices = (int *)malloc(link_array_length*sizeof(int));
//     population->link_alleles_indices = (int *)malloc(link_array_length*sizeof(int));
//     //
//     // copy local arrays to the population arrays ...
//     memcpy(population->allelicity,allelicity_numbers,(size_t)population->loci*sizeof(int));
//     memcpy(population->freqs,allele_freqs,(size_t)population->independent_alleles*sizeof(double)); 
//     // linkages info copying ...
//     memcpy(population->link_alleles_number,link_alleles_number,(size_t)linkage_number*sizeof(int));
//     memcpy(population->link_freq,link_freq,(size_t)linkage_number*sizeof(double));
//     memcpy(population->link_loci_indices,link_loci,(size_t)link_array_length*sizeof(int));
//     memcpy(population->link_alleles_indices,link_alleles,(size_t)link_array_length*sizeof(int));
//     //
//     return;
// }





// // generate constraint matrix ...
// void generate_constraint_matrix(PopDes *population, SparseMat *matrix){
//     //
//     // calculate non-zero elements in the sparse matrix string the constraints ...
//     int length = 0;
//     // first, just the allele frequencies constraints:
//     for (int l = 0; l < population->loci; l++){
//         // number of independent equation for the locus BY ...
//         // number of combination with a given allele at this locus:
//         length += (population->allelicity[l]-1) * (population->haplotypes/(population->allelicity[l]));
//     }
//     // normalization condition, sum of all haplotype fractions is 1:
//     length += population->haplotypes;
//     // linkage related items to be added ...
//     //
//     // allocating Sparse Matrix related array given the length:
//     matrix->nonzero = length;
//     matrix->rows = population->alleles_total+1; // mind linear independence & normalization (linkage - ?)
//     matrix->cols = population->haplotypes; // presumably all haplotypes are involved, it's the case 
//     matrix->irow = (int *)malloc((length+1)*sizeof(int));
//     matrix->jcol = (int *)malloc((length+1)*sizeof(int));
//     matrix->mval = (double *)malloc((length+1)*sizeof(double));
//     //
//     // we'd also need two more running indicies here ///
//     int row = 0; // row index of the constraint matrix ...
//     int idx = 1; // index along the sparse matrix arrays ...
//     for (int l = 0; l < population->loci; l++){
//         // allele frequencies obey normalization, so use just allelicity-1 of them as a constraint ...
//         for (int a = 0; a < population->allelicity[l]-1; a++){
//             // scan the entire set of haplotypes ...
//             for (int h = 0; h < population->haplotypes; h++){            
//                 // checking if haplotypes match ...
//                 if (population->haplotype_description[h][l] == a){
//                     // indexes must updated for the glpk - Fortran-style indexing ...
//                     matrix->irow[idx] = row+1; // that's where we need that row index of that actual constraint matrix ...
//                     matrix->jcol[idx] = h+1; // column index is simply the haplotype number ...
//                     matrix->mval[idx] = 1.0;
//                     idx++; // count every nonzero element ...
//                 }
//             }
//             row++; // constraints matrix row index incremented every other accounted allele ... 
//         }
//     }
//     //
//     // now, should add the normalization condition here ...
//     for (int h = 0; h < population->haplotypes; h++){
//         // indexes must updated for the glpk - Fortran-style indexing ...
//         matrix->irow[idx] = row+1; // that's where we need that row index of that actual constraint matrix ...
//         matrix->jcol[idx] = h+1; // column index is simply the haplotype number ...
//         matrix->mval[idx] = 1.0;
//         idx++; // count every nonzero element ...
//     }
//     row++;
//     //
//     //
//     //
//     // folowing assertions must be met:
//     assert(matrix->rows == row);
//     assert(matrix->nonzero+1 == idx);
//     //
//     return;
// }





























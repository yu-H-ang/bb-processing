#include "main.h"

// Define global variables declared in header file
char ROOT_DIR[] = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization";
char RESULT_ROOT_DIR[] = "/home-4/yzhan175@jhu.edu/scratch/fluidized_bed/fluidization/processing";
double tStart = 0.;
double tEnd = 10000.;
int tt;

int main(void)
{
    double *volz;
    
    // Read and sort output directory for finding files within our time limits
    //init_part_files();
    init_flow_files();

    // Create output directory
    create_output_dir();

    // Get number of particles
    //nparts = cgns_read_nparts();

    // Initialize domain and flow arrays
    domain_init();

    // Initialize partstruct and flow vars
    //parts_init();
    //flow_init();
    uf = (double*) malloc(dom.Gcc.s3 * sizeof(double));
    vf = (double*) malloc(dom.Gcc.s3 * sizeof(double));
    wf = (double*) malloc(dom.Gcc.s3 * sizeof(double));
    phase = (int*) malloc(dom.Gcc.s3 * sizeof(int));
    for (int i = 0; i < dom.Gcc.s3; i++) {
        phase[i] = -2;
    }
    volz = (double*) malloc(dom.Gcc.kn * sizeof(double));
    
    // open file
    char filename[CHAR_BUF_SIZE] = "";
    sprintf(filename, "%s/output.txt", RESULT_ROOT_DIR);
    FILE *fdat = fopen(filename, "w");
    if (fdat == NULL) {
        printf("Error opening file %s!\n", filename);
        exit(EXIT_FAILURE);
    }
    
    // loop through all files
    for(tt = 0; tt < nFiles; tt++) {
        
        // Read in new data
        printf("%d\n", tt);
        cgns_fill_flow();
        
        // calculation
        for(int k = 0; k < dom.Gcc.kn; k++) {
            volz[k] = 0.;
            for(int j = 0; j < dom.Gcc.jn; j++) {
                for(int i = 0; i < dom.Gcc.in; i++) {
                    volz[k] += (phase[i+j*dom.Gcc.s1+k*dom.Gcc.s2] > -1);
                }
            }
            volz[k] /= (dom.Gcc.in*dom.Gcc.jn);
            
            // write to file
            fprintf(fdat, "%f ", volz[k]);
        }
    }
    
    // close file
    fclose(fdat);
    // write phaseAvereaged to file
    //write_part_data();
    
    // Free and exit
    free(volz);
    //free_vars();
    /*
    for (int i = 0; i < nFiles; i++) {
        free(partFiles[i]);
    }
    free(partFiles);
    free(partFileMap);
    free(partFileTime);
    */
    for (int i = 0; i < nFiles; i++) {
        free(flowFiles[i]);
    }
    free(flowFiles);
    free(flowFileMap);
    free(flowFileTime);
    free(parts);
    free(uf);
    free(vf);
    free(wf);
    free(phase);
    return EXIT_SUCCESS;
}

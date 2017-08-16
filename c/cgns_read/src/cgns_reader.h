#ifndef _READER_H
#define _READER_H

#include "main.h"
#include <cgnslib.h>
#include <dirent.h>

/**** STRUCTURES ****/
// part_struct
typedef struct part_struct {
  double r;
  double x;
  double y;
  double z;
  double u;
  double v;
  double w;
} part_struct;

// grid_info
typedef struct grid_info {
  int is;
  int js;
  int ks;
  int ie;
  int je;
  int ke;
  int in;
  int jn;
  int kn;
  int s1;
  int s2;
  int s3;
} grid_info;

// dom_struct
typedef struct dom_struct {
  grid_info Gcc;
  double xs;
  double ys;
  double zs;
  double xe;
  double ye; 
  double ze;
  double xl;
  double yl;
  double zl;
  int xn;
  int yn;
  int zn;
  double dx;
  double dy;
  double dz;
} dom_struct;

// boundary condition structure
typedef struct BC {
  int pW;
  int pE;
  int pS;
  int pN;
  int pB;
  int pT;
} BC;


/**** VARIABLES ****/
// File Variables
extern int nFiles;
extern char **partFiles;
extern double *partFileTime;
extern int *partFileMap;
extern char **flowFiles;
extern double *flowFileTime;
extern int *flowFileMap;

// Flow Variables
extern double *uf;
extern double *vf;
extern double *wf;
extern int *phase;

// Number of Particles
extern int nparts;

// host and dev part_struct parts;
extern part_struct *parts;

// host and dev dom_struct doms
extern dom_struct dom;

// host and dev bc struct bcs
extern BC bc;

/**** FUNCTIONS ****/
// read input file
void main_read_input(void);

// read and sort part/flow files
void init_part_files(void);
void init_flow_files(void);

// Merge sort
void merge_sort(double *A, int n, int *A2);
void merge(double *A, int n, int m, int *A2);

// read CUDA_VISIBLE_DEVICES
int read_devices(void);

// Create directory for output data
void create_output_dir(void);

// Read nparts
int cgns_read_nparts(void);

// initialize part_struct and flow vars
void parts_init(void);
void flow_init(void);

// initialize dom_struct binDom
void domain_init(void);

// Read part_struct data
void cgns_fill_parts(void);

// Read flow_struct data
void cgns_fill_flow(void);

// show binDom and bc structures
void show_domain(void);

// get sigfigs of last file 
void get_sigfigs(void);

// write each timestep
void write_part_data(void);

// Free parts
void free_vars(void);

// txt operations
void read_txt(char* fname, double *arr, int len);
void write_txt(char *fname, double *arr, int len);

#endif

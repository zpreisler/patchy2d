#ifndef INIT_H
#define INIT_H
#define CONF_STAMP "#Configurational_file"
#include "params.h"
extern int init_random(header *t);
extern int init_configuration(char *file,header *t);
extern int save_configuration(char *file,header *t);
extern int load_configuration(char *file,header *t);
extern void copy_configuration(header *t,__m128d box);
extern input_files *find_configurational_files(int argc,char *argv[]);
//extern input_files *mpi_find_configurational_files(mpi_world *mpi,int argc,char *argv[]);
extern void print_input_configurational_files(FILE *f,input_files *input);
extern int load_configuration_file(FILE *f,header *t);
extern int init_configuration_random(header *t);
extern int read_input(int argc,char *argv[],input_files *input,header *t);
extern void reset_particle(compound_particle *c,header *t);
#endif

#ifndef PARAMS_H
#define PARAMS_H
#define DIR_LENGTH 32
#define SEED 1346
#define MPI_SEED 913
#include "zargs.h"
/*typedef struct mpi_world{
	int numtasks;
	int rank;
	int length;
	int provided;
	char processor_name[NAME_LENGTH];
	unsigned int seed;
	int waiting_for_reply_from;
	int blocked;
}mpi_world;*/
typedef struct input_files{
	int n;
	FILE **f;
	char **file;
}input_files;
//extern mpi_world *init_mpi_world(void);
extern header *init_header(void);
#endif

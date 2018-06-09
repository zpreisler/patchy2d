#ifndef ZARGS_H
#define ZARGS_H
#define ARG_STRING 1
#define ARG_INT 2
#define ARG_LONG 3
#define ARG_DOUBLE 4
#define ARG_LFVEC 5
#define ARG_LFVEC2 6
#define ARG_SPECIES 7
#define ARG_UPDATE 8
#define NAME_LENGTH 256
#define COMMENT '#'
#include <immintrin.h>
typedef struct arg{
	char short_opt,*long_opt,*key;
	void *value;
	unsigned flag;
	size_t mask;
	size_t off;
	size_t n;
	char *help;
	void *sub;
}arg;
typedef struct arg_set{
	int n;
	arg *args;
}arg_set;
typedef struct patch{
	int id;
	double *angle; //angle with respect to particle orientation
	__m128d *q_angle; //
	__m128d *q;	
}patch;
typedef struct particle{
	//Energy
	unsigned int en_new,en_old;	
	unsigned nd;
	__m128d *q; //position
	__m128d *q_tmp; //old position
	__m128d *q_track; //track
	__m128d *q_well; //well
	__m128d *or; //orientation
	__m128d *or_well; //well orientation
	__m128d *qp_rij; //q-p orientational vector
	double qp_r2;
	//Parameters
	double sigma,sigma_well;
	//Patches
	patch *patch;
	int npatch;
	double patch_width;
	//ID
	unsigned n;
	void *specie;
	//HASH
	unsigned h;
	double *nd_d2;
	struct particle **nd_list;
	struct particle **new_list,**old_list;
	struct particle *next,**prev;
	//Flag
	int flag;
	int type;
	//Graphs
	int pass;
	int idx;
	int npcycles; //number of pcycles
	void *pcycles; //pcycle
}particle;
typedef struct hash_table{
	particle *p;
	particle ***list;
}hash_table;
typedef struct update{
	struct update *next;
	void *value;
}update;
typedef struct species{
	struct species *next;
	unsigned int N,Nalloc; //number of particles
	int flag;
	double sigma; //particle diameter
	double sigma_well; //particle diameter plus patch range
	int npatch;
	double patch_width; //patch width for kern-frenkel potentials
	char patch_type[16];
	double patch_angle;
	double *angles;
	int grand_canonical; //grand canonical switch
	double mu; //chemical potential 
	double *interaction_matrix;
	particle *p; //pointer to particle list
}species;
typedef struct header{
	char name[NAME_LENGTH];
	unsigned int N,Nalloc;
	unsigned int npatch,npatch_alloc;
	long long int step;
	species *specie;
	species *update;
	particle *p;
	patch *s;
	int nspecies;
	//HASH
	int ndir;
	__m128i *dir;
	__m128i h;
	__m128d h1;
	__m128d hash_cell;
	__m128d box;
	__m128d copy;
	hash_table *table;
	//Arguments
	arg_set argz;
	//Parameters
	double epsilon,end_epsilon;
	double pressure;
	double uy;
	//Printing and writing modulus
	int mod,pmod;
	int optimize;
	int snapshot;
	int verbose;
	//Max displacment
	__m128d max_displacement;
	double max_rotation;
	double max_vol;
	double max_xy;
	double max_dxdy;
	double max_uy;
	double max_dsigma;
	//Random seed
	int seed;
	//Einstein
	double lambda;
	double lambda_coupling;
	__m128d q_dcm;
	//Graphs
	int ncycles;//number of cycles detected
	particle **matrix_list;
	void *adj_matrix;
	void *graph_pcycle;
	time_t init_time;
	long long int steps_passed;
}header;
extern int command_line_args(int argc,char *argv[],arg_set argz);
extern void dump_args(FILE *f,arg_set argz);
extern void usage(FILE *f,arg_set argz,char *prog);
extern void set_args(header *t);
extern int file_args(FILE *f,arg_set argz);
extern int read_string(FILE *f,char *s,int n);
extern int read_stamp(FILE *f,char *s,int n);
#endif

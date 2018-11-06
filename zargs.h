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
	//linking
	/////////
	void *p;
	void *c;
	void *specie;
	//
	int id;
	double *angle; //angle with respect to particle orientation
	__m128d *q_angle; //
	__m128d *q;	
}patch;
//structure for a particle with patches
///////////////////////////////////////
typedef struct particle{
	patch *patch; //pointer to array of patches
	void *c; //pointer to the compound
	void *specie; //poiter to the beging of the array of the particles of the same specie
	//Energy
	////////
	unsigned int en_new,en_old;	
	unsigned nd; // number of cell neighbours
	double *nd_r2;
	__m128d *nd_rij;
	//
	__m128d *q; //position
	__m128d *q_tmp; //old position
	__m128d *q_track; //track -- for Einstein solid calculations
	__m128d *q_well; //well -- for Einstein solid calculations
	__m128d *or; //orientation
	__m128d *or_well; //well orientation -- for Einstein solid
	__m128d *qp_rij; //q-p orientational vector
	double qp_r2; //save q-p r^2 distance
	//for compounds
	///////////////
	__m128d *dq; //vector from the center of the compound
	//Parameters
	////////////
	double sigma,sigma_well; //particle diameter and interaction range
	//Patches
	/////////
	int npatch; // number of patches
	double patch_width; // patch witdh (if patch width kept fixed)
	//id
	////
	unsigned n; //ID
	//HASH
	//////
	unsigned h;
	//double *nd_d2;
	struct particle **nd_list;
	struct particle **new_list,**old_list;
	struct particle *next,**prev;
	//Flag
	//////
	int flag;
	int type;
	//Graphs
	////////
	int pass; // tags visited particles
	int idx; // particle index
	int id; //id
	int npcycles; //number of pcycles -- for cycle analysis
	void *pcycles; //pcycle -- pointer to the list of cycles
	//Colors
	////////
	float *color;
}particle;
//structure for a compound particle
///////////////////////////////////
typedef struct compound_particle{
	patch *s; //pointer to patches belonging to the compound
	particle *p; //pointer to particles belonging to the compound
	void *specie; //pointer to the beging of the array of the compounds of the same specie
	int id;
	int nparticle;
	int npatch;
	unsigned int flag,n; //flag, ID;
	__m128d *q; //compound position;
	__m128d *q_tmp; //compound position;
	__m128d *q_well;
	__m128d *or; //compound orientation;
	__m128d *or_well;
	void *cluster;
}compound_particle;
//hash table structure
typedef struct hash_table{
	particle *p;
	particle ***list;
}hash_table;
typedef struct update{
	struct update *next;
	void *value;
}update;
typedef struct species{
	//linking
	/////////
	patch *s; //pointer to patch array
	particle *p; //pointer to particle array
	compound_particle *c; //pointer to compound particle array
	//
	unsigned int compound; //0: simple particle; 1: it is a particle compound
	unsigned int nppc; //number of particles per compound
	double rod_length;
	unsigned int ncompound,ncompound_alloc; //number of compounds
	unsigned int nparticle,nparticle_alloc; //number of particles
	unsigned int npatches,npatches_alloc; //number of all patches
	//N -- read from input line or command line
	///////////////////////////////////////////
	unsigned int N;//,Nalloc; //number of particles, number of allocated particles
	struct species *next;
	int flag;
	//particle specs (if the individual particles are not modified)
	///////////////////////////////////////////////////////////////
	double sigma; //particle diameter
	double sigma_well; //particle diameter plus patch range
	int npatch; //number of patches per particle
	double patch_width; //patch width for kern-frenkel potentials
	char patch_type[256];
	double patch_angle;
	double *angles;
	//Simulation settings for the specie
	////////////////////////////////////
	int grand_canonical; //grand canonical switch
	double mu; //chemical potential 
	double *interaction_matrix;
}species;
typedef struct cluster{
	compound_particle **p; //pointer to a pointer in cluster_list
	int n; //number of compound particles in the cluster
}cluster;
typedef struct cluster_list{
	int n; //number of stored compounds
	int ncluster; //number of clusters
	double avg_size; //average cluster size
	double max_size; //max cluster size
	compound_particle **c; //store compound particles
	cluster *clusters;
}cluster_list;
typedef struct files{
	FILE *fepsilon; // inverse temperature
	FILE *fmu;
	FILE *fpressure;
	FILE *fen; // U/N
	FILE *frho; // N/V
	FILE *fvol; //volume
	FILE *fn; //number of particles 
	FILE *ftime; //passed time
	FILE *fstep; //passed MC steps
	FILE *fncluster; //number of clusters
	FILE *fcluster_avg_size; //average cluster size
	FILE *fcluster_max_size; //max cluster size
}files;
typedef struct header{
	char name[NAME_LENGTH];
	//unsigned int N,Nalloc;// number of particles, number of particles allocated
	unsigned int ncompound,ncompound_alloc; //number of compounds, number to compounds allocated
	unsigned int nparticle,nparticle_alloc; //number of particles, number of particles allocated
	unsigned int npatches,npatches_alloc;//number of all patches, number of patches allocated
	long long int step;
	species *specie;
	species *update;
	compound_particle *c;
	particle *p;
	patch *s;
	//Clusters
	//////////
	cluster_list *cluster;
	int nspecies;
	//HASH
	//////
	int ndir;
	__m128i *dir;
	__m128i h;
	__m128d h1;
	__m128d hash_cell;
	__m128d box;
	__m128d copy;
	hash_table *table;
	//Arguments
	///////////
	arg_set argz;
	//Energy
	////////
	int energy;
	//Parameters
	////////////
	double epsilon,end_epsilon;
	double pressure;
	double uy;
	//Time
	//////
	time_t t1,t2;
	//Printing and writing
	//////////////////////
	int mod,pmod;
	int optimize;
	int snapshot;
	int verbose;
	int display;
	int explore;
	__m128d screen_geometry;
	//Ensembles
	///////////
	int npt;
	int nptxy;
	int nptdxdy;
	int shape;
	//Max displacments
	//////////////////
	__m128d max_displacement;
	double max_rotation;
	double max_vol;
	double max_xy;
	double max_dxdy;
	double max_uy;
	double max_dsigma;
	//Random seed
	/////////////
	int seed;
	//Files
	///////
	files file;
	//Einstein 
	//////////
	double lambda;
	double lambda_coupling;
	__m128d q_dcm;
	//Graphs
	////////
	int ncycles;//number of cycles detected
	particle **matrix_list;
	void *adj_matrix;
	void *graph_pcycle;
	time_t init_time;
	long long int steps_passed;
	//Colors
	////////
	float *particle_colors;
}header;
extern int command_line_args(int argc,char *argv[],arg_set argz);
extern void dump_args(FILE *f,arg_set argz);
extern void usage(FILE *f,arg_set argz,char *prog);
extern void set_args(header *t);
extern int file_args(FILE *f,arg_set argz);
extern int read_string(FILE *f,char *s,int n);
extern int read_stamp(FILE *f,char *s,int n);
#endif

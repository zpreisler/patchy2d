#ifndef ALLOC_H
#define ALLOC_H
#define NLIST 48
#define MIN_ALLOC 4096
#define M_ALLOC 2
#define max(a,b) ((a)>(b)?(a):(b))
typedef struct patch_memory{
	double *angle;
	__m128d *q_angle;
	__m128d *q;
}patch_memory;
typedef struct particle_memory{
	__m128d *q;
	__m128d *q_tmp;
	__m128d *q_track;
	__m128d *q_well;
	__m128d *qp_rij;
	__m128d *or;
	__m128d *or_well;
}particle_memory;
typedef struct compound_memory{
	__m128d *q; //compound position;
	__m128d *q_tmp;
	__m128d *q_well;
	__m128d *or; //compound orientation;
	__m128d *or_well;
}compound_memory;
extern particle *alloc_particles(header *t);
#endif

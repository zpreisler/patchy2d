#ifndef GRAND_CANONICAL
#define GRAND_CANONICAL
//extern int grand_canonical(mpi_world *mpi __attribute__((unused)),header *t);
//extern int grand_canonical_n(mpi_world *mpi __attribute__((unused)),header *t);
extern int pre_insert_particle(header *t,particle *q,int *en);
extern int post_insert_particle(header *t,particle *q);
extern int delete_particle(header *t,particle *p,particle *q);
extern int delete_compound(header *t,compound_particle *c,compound_particle *d);
extern int pre_insert_compound(header *t,compound_particle *c,int *en);
extern int post_insert_compound(header *t,compound_particle *c);
extern int mc_gc(header *t,int *en);
#endif

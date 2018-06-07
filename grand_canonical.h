#ifndef GRAND_CANONICAL
#define GRAND_CANONICAL
//extern int grand_canonical(mpi_world *mpi __attribute__((unused)),header *t);
//extern int grand_canonical_n(mpi_world *mpi __attribute__((unused)),header *t);
extern int mc_gc(header *t,int *en);
#endif

#ifndef CANONICAL_H
#define CANONICAL_H
#define SWAP_ORDER(a,b) {int (c)=(a);(a)=(b);(b)=(c);}
extern int canonical(mpi_world *mpi,header *t);
extern int mpi_canonical(mpi_world *mpi,header *t);
extern int nvt_dsigma(mpi_world *mpi __attribute__((unused)),header *t);
extern int mc_rotate(particle *p,header *t,int *en);
extern int mc_move(particle *p,header *t,int *en);
extern unsigned int rnd_rank(mpi_world *mpi);
extern particle *rnd_particle(header *t);
extern particle *rnd_specie(species *s);
extern int mc_move(particle *p,header *t,int *en);
extern int mc_rotate(particle *p,header *t,int *en);
#endif

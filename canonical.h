#ifndef CANONICAL_H
#define CANONICAL_H
#define SWAP_ORDER(a,b) {int (c)=(a);(a)=(b);(b)=(c);}
//extern int canonical(mpi_world *mpi,header *t);
//extern int mpi_canonical(mpi_world *mpi,header *t);
//extern int nvt_dsigma(mpi_world *mpi __attribute__((unused)),header *t);
//extern int mc_rotate(particle *p,header *t,int *en);
//extern int mc_move(particle *p,header *t,int *en);
//extern unsigned int rnd_rank(mpi_world *mpi);
extern particle *rnd_particle(header *t);
extern particle *rnd_specie(species *s);
//extern int mc_move(particle *p,header *t,int *en);
//extern int mc_rotate(particle *p,header *t,int *en);
//
extern void c_list_swap(compound_particle *c);
extern void c_adjust_lists(compound_particle *c);

extern int c_en_old(compound_particle *c);
extern int c_en_new(compound_particle *c);

extern int mc_move(compound_particle *c,header *t,int *en);
extern int mc_rotate(compound_particle *c,header *t,int *en);
#endif

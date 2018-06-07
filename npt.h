#ifndef NPT_H
#define NPT_H
//extern int npt(mpi_world *mpi __attribute__((unused)),header *t);
//extern int mpi_npt(mpi_world *mpi __attribute__((unused)),header *t);
//extern int npt_xy(mpi_world *mpi __attribute__((unused)),header *t);
extern int mc_npt(header *t,int *en);
extern int mc_npt_xy(header *t,int *en);
extern int mc_npt_dxdy(header *t,int *en);
extern int mc_uy(header *t,int *en);
#endif

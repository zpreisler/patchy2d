#ifndef ENERGY_H
#define ENERGY_H
extern int overlap(particle *p,header *t);
extern int particle_energy_hash(particle *p,header *t);
extern int particle_energy_hash2(particle *p);
extern int all_particle_energy_hash(header *t,int *en);
extern int checksum(FILE *f,header *t,int energy);
#endif

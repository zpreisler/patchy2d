#ifndef ENERGY_H
#define ENERGY_H
extern int overlap(particle *p,header *t);
extern int compound_overlap(compound_particle *c,header *t);
extern int particle_energy_hash(particle *p,header *t);
extern int particle_energy_hash2(particle *p);

extern int particle_energy_hash2m(particle *p,header *t);

extern int all_particle_energy_hash(header *t,int *en);
extern int checksum(FILE *f,header *t,int energy);
#endif

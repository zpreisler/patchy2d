#ifndef RUN_H
#define RUN_H
extern int run(header *t,mySDL *s);
extern int mc_rotate(particle *p,header *t,int *en);
extern int mc_move(particle *p,header *t,int *en);
extern compound_particle *rnd_compound(header *t);
extern particle *rnd_particle(header *t);
extern particle *rnd_specie(species *s);
#endif

#ifndef RUN_H
#define RUN_H

#if defined SDL
extern int run(header *t,mySDL *s);
#else
extern int run(header *t);
#endif

extern int mc_rotate(particle *p,header *t,int *en);
extern int mc_move(particle *p,header *t,int *en);
#endif

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <utils.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include "dSFMT.h"
#include "zargs.h"
#include "params.h"
#include "mm_math.h"
#include "hash.h"
#include "energy.h"
#include "lists.h"
#include "init.h"
#include "patches.h"
#include "mpi.h"
extern dsfmt_t dsfmt;
particle *rnd_particle(header *t){
	species *s=t->specie;
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*t->nparticle);
	while(s&&s->nparticle<(rnd+1)){
		rnd-=s->nparticle;
		s=s->next;
	}
	return (particle*)s->p+rnd;
}
compound_particle *rnd_compound(header *t){
	species *s=t->specie;
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*t->ncompound);
	while(s&&s->ncompound<(rnd+1)){
		rnd-=s->ncompound;
		s=s->next;
	}
	return (compound_particle*)s->c+rnd;
}
particle *rnd_specie(species *s){
	unsigned int rnd=(unsigned)(dsfmt_genrand_open_open(&dsfmt)*s->nparticle);
	return (particle*)s->p+rnd;
}

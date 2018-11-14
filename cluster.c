#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>
#include <math.h>
#include <utils.h>
#include "dSFMT.h"
#include "zargs.h"
#include "hash.h"
#include "mm_math.h"
#include "cluster.h"
extern dsfmt_t dsfmt;
int new_cluster(compound_particle *c,header *t){
	int i=t->cluster->n++;
	int j=t->cluster->ncluster++;	
	cluster *cc;
	*(t->cluster->c+i)=c;
	cc=t->cluster->clusters+j;
	cc->n=1;
	cc->p=t->cluster->c+i;
	return 0;
}
int add2cluster(compound_particle *c,header *t){
	int i=t->cluster->n++;
	int j=t->cluster->ncluster;	
	cluster *cc;
	*(t->cluster->c+i)=c;
	cc=t->cluster->clusters+j-1;
	cc->n++;
	c->cluster=(cluster*)cc;
	return 0;
}
int find_cluster(compound_particle *c,header *t){
	int i;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		cluster_check_particle(p,t);
	}
	return 0;
}
int cluster_check_particle(particle *p,header *t){
	unsigned int h;
	int k;
	__m128d rij;
	double r2,d2;
	double dd,w=cos(15.0/180.0*M_PI);
	particle *q;
	compound_particle *c,*d;
	h=hash(*p->q,t->h1);
	p->nd=0;
	for(k=0;k<t->ndir;k++){
		for(q=*(t->table)[h].list[k];q;q=q->next){
			if(q->c!=p->c){ //Only particles belonging to different compounds
				if(q!=p){
					c=q->c;
					if(!c->cluster){
						//check if compound is already part of a cluster
						////////////////////////////////////////////////
						rij=_mm_dist_uy(*(p->q),*(q->q),t->box,t->uy);
						r2=length2(rij);
						d2=SQR((p->sigma_well+q->sigma_well)*0.5);
						//check if the particle is within a specified radius
						////////////////////////////////////////////////////
						if(r2<d2){
							d=p->c;
							//add cluster to the cluster list	
							dd=dot(*c->or,*d->or);
							//printf("%lf\n",dd);
							//if(fabs(dd)>w){
								//Checking the orientation
								add2cluster(c,t);
								find_cluster(c,t);
							//}
						}
					}
				}
			}
		}
	}
	return 0;
}
int clusters_reset(header *t){
	unsigned int i;
	t->cluster->ncluster=0;
	t->cluster->n=0;
	compound_particle *c;
	for(i=0;i<t->ncompound;i++){
		c=t->c+i;
		c->cluster=NULL;
	}
	return 0;
}
int find_new_cluster(compound_particle *c,header *t){
	int i;
	particle *p;
	new_cluster(c,t);
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		cluster_check_particle(p,t);
	}
	return 0;
}
int find_all_clusters(header *t){
	unsigned int i;
	clusters_reset(t);
	compound_particle *c;
	for(i=0;i<t->ncompound;i++){
		c=t->c+i;
		if(!c->cluster){
			find_new_cluster(c,t);
		}
	}
	avg_max_cluster_size(t);
	return 0;
}
int print_clusters(header *t){
	int i,j;
	printf("number of clusters [%d]\n",t->cluster->ncluster);
	cluster *cc;
	compound_particle *c;
	for(i=0;i<t->cluster->ncluster;i++){
		cc=t->cluster->clusters+i;
		printf("[%d][%d]\n",i,cc->n);
		for(j=0;j<cc->n;j++){
			c=*(cc->p+j);
			printf("%p\n",c);
		}
	}
	return 0;
}

void avg_max_cluster_size(header *t){
	int i;
	int all=0;
	t->cluster->max_size=1;
	cluster *cc;
	for(i=0;i<t->cluster->ncluster;i++){
		cc=t->cluster->clusters+i;
		all+=cc->n;
		if(cc->n>t->cluster->max_size)t->cluster->max_size=cc->n;
	}
	t->cluster->avg_size=(double)all/(double)t->cluster->ncluster;
}

#include <stdlib.h>
#include <stdio.h>
#include <immintrin.h>
#include <math.h>
#include <utils.h>
#include "zargs.h"
#include "hash.h"
#include "mm_math.h"
int overlap(particle *p,header *t){
	unsigned int h;
	int k;
	__m128d rij;
	double r2,d2;
	particle *q;
	h=hash(*p->q,t->h1);
	p->nd=0;
	for(k=0;k<t->ndir;k++){
		for(q=*(t->table)[h].list[k];q;q=q->next){
			if(q->c!=p->c){ //Only particles belonging to different compounds
				if(q!=p){
					rij=_mm_dist_uy(*(p->q),*(q->q),t->box,t->uy);
					r2=length2(rij);
					d2=SQR((p->sigma+q->sigma)*0.5);
					//check -*- hard sphere overlap -*-
					///////////////////////////////////
					if(r2<d2){
						return 1;
					}
					p->nd_r2[p->nd]=r2;
					p->nd_rij[p->nd]=rij;
					p->nd_list[p->nd++]=q;
					*q->qp_rij=rij;
					q->qp_r2=r2;
				}
			}
		}
	}
	return 0;
}
int compound_overlap(compound_particle *c,header *t){
	int i;
	particle *p;
	for(i=0;i<c->nparticle;i++){
		p=c->p+i;
		if(overlap(p,t))return 1;
	}
	return 0;
}
int specific_interaction(int i,int j,int n,double *m){
	if(*(m+i+n*j)!=0.0)return 1;
	return 0;
}
int patch_energy(particle *p,particle *q,double r2){
	int j,l;
	int b=0;
	double c=1.0/sqrt(r2);
	if(q->c==p->c)return 0; //Only particles belonging to different compounds
	__m128d u=*q->qp_rij*c;
	__m128d v=*q->qp_rij*-c;
	for(j=0;j<p->npatch;j++){
		if(dot(*(p->patch+j)->q,u)>p->patch_width){
			for(l=0;l<q->npatch;l++){
				if(dot(*(q->patch+l)->q,v)>q->patch_width){
					p->new_list[p->en_new++]=q;
					b++;
				}
			}
		}
	}
	return b;
}
//if(p->flag==q->flag)return 0;
//#if defined(TYPE)
//	if(p->type==q->type)return 0;
//#endif
/*#if defined(SPECIFIC)
				//if((j==0&&l==2)||(j==2&&l==0)||(j==1&&l==1)){
				if(specific_interaction(j,l,q->npatch,((species*)p->specie)->interaction_matrix)){
					if(dot(*(q->patch+l)->q,v)>q->patch_width){
						p->new_list[p->en_new++]=q;
						b++;
					}
				}
#else*/
//#endif
int patch_energy_all(particle *p,particle *q,double r2){
	int j,l;
	int b=0;
	double c=1.0/sqrt(r2);
	if(q->c==p->c)return 0; //Only particles belonging to different compounds
	__m128d u=*q->qp_rij*c;
	__m128d v=*q->qp_rij*-c;
	for(j=0;j<p->npatch;j++){
		if(dot(*(p->patch+j)->q,u)>p->patch_width){
			for(l=0;l<q->npatch;l++){
				if(dot(*(q->patch+l)->q,v)>q->patch_width){
					p->new_list[p->en_new++]=q;
					q->new_list[q->en_new++]=p;
					b++;
				}
			}
		}
	}
	return b;
}
//if(p->flag==q->flag)return 0;
//#if defined(TYPE)
//	if(p->type==q->type)return 0;
//#endif
/*#if defined(SPECIFIC)
				//if((j==0&&l==2)||(j==2&&l==0)||(j==1&&l==1)){
				if(specific_interaction(j,l,q->npatch,((species*)p->specie)->interaction_matrix)){
					if(dot(*(q->patch+l)->q,v)>q->patch_width){
						p->new_list[p->en_new++]=q;
						q->new_list[q->en_new++]=p;
						b++;
					}
				}
#else*/
//#endif
int particle_energy_hash(particle *p,header *t){
	unsigned int h;
	int k;
	int b;
	double r2,d2;
	__m128d rij;
	particle *q;
	p->en_new=0;
	h=hash(*p->q,t->h1);
	unsigned int count=0;
	for(k=0,b=0;k<t->ndir;k++){
		for(q=*(t->table)[h].list[k];q;q=q->next){
			if(q->c!=p->c){ //only particle of different compounds interact
				count++;
				if(p!=q){
					rij=_mm_dist_uy(*(p->q),*(q->q),t->box,t->uy);
					r2=length2(rij);
					*q->qp_rij=rij;
					d2=SQR((p->sigma_well+q->sigma_well)*0.5);
					if(r2<d2){
						b+=patch_energy(p,q,r2);
					}
				}
			}
		}
	}
	if(p->nd!=count){
		printf("p->nd %d  count %d\n",p->nd,count);
	}
	return b;
}
/*
int particle_energy_hash2(particle *p){
	unsigned int i;
	int b;
	double r2,d2;
	particle *q;
	p->en_new=0;
	for(i=0,b=0;i<p->nd;i++){
		q=p->nd_list[i];
		//if(q->c!=p->c){ //only particle of different compounds interact
			r2=q->qp_r2;
			d2=SQR((p->sigma_well+q->sigma_well)*0.5);
			if(r2<d2){
				b+=patch_energy(p,q,r2);
			}
		//}
	}
	return b;
}
*/
int particle_energy_hash2(particle *p){
	unsigned int i;
	int b;
	double r2,d2;
	particle *q;
	//__m128d rij;
	p->en_new=0;
	for(i=0,b=0;i<p->nd;i++){
		q=p->nd_list[i];
		r2=p->nd_r2[i];
		*q->qp_rij=p->nd_rij[i];
		d2=SQR((p->sigma_well+q->sigma_well)*0.5);
		if(r2<d2){
			b+=patch_energy(p,q,r2);
		}
	}
	return b;
}
		//if(q->c!=p->c){ //only particle of different compounds interact
			//rij=_mm_dist_uy(*(p->q),*(q->q),t->box,t->uy);
			//r2=length2(rij);
			//*q->qp_rij=rij;
			//r2=q->qp_r2;
		//}
int all_particle_energy_hash(header *t,int *en){
	unsigned int i,h;
	int k;
	__m128d rij;
	double r2,d2;
	particle *p,*q;
	species *s=t->specie;
	int b=0;
	*en=0;
	while(s){
		p=s->p;
		for(i=0;i<s->nparticle;i++){
			(p+i)->en_new=0;
		}
		s=s->next;
	}
	s=t->specie;
	while(s){
		for(i=0;i<s->nparticle;i++){
			p=(particle*)s->p+i;
			h=hash(*p->q,t->h1);
			for(k=0;k<t->ndir;k++){
				for(q=*(t->table)[h].list[k];q;q=q->next){
					if(q->c!=p->c){ //Only particles belonging to different compounds
						if(q->n>p->n){
							rij=_mm_dist_uy(*(p->q),*(q->q),t->box,t->uy);
							r2=length2(rij);
							*q->qp_rij=rij;
							q->qp_r2=r2;
							if(r2<SQR((p->sigma+q->sigma)*0.5)){ //overlap
								return -1;
							}
							d2=SQR((p->sigma_well+q->sigma_well)*0.5);
							if(r2<d2){
								b+=patch_energy_all(p,q,r2);
							}
						}
					}
				}
			}
		}
		s=s->next;
	}
	(*en)=b;
	return b;
}
int checksum(FILE *f,header *t,int energy){
	int energy_check=0;
	int overlap;
	overlap=all_particle_energy_hash(t,&energy_check);
	if(overlap==-1)fprintf(f,"detected "RED"overlap\n"RESET);
	fprintf(f,RESET"energy="RED"%d"RESET" energy_check="RED"%d\n"RESET,energy,energy_check);
	if(energy==energy_check){
		fprintf(f,RESET"Checksum: "GREEN"pass\n"RESET);
		return 0;
	}
	else{
		fprintf(f,RESET"Checksum: "RED"fails\n"RESET);
		return 1;
	}
}

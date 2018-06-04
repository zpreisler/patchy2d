#include <stdlib.h>
#include <stdio.h>
#include <utils.h>
#include "params.h"
#include "zargs.h"
#include "lists.h"
void list_swap(particle *p){
	LIST_SWAP(p->new_list,p->old_list);
	SWAP(p->en_old,p->en_new);
}
void adjust_lists(particle *p){
	unsigned int i,j,in;
	particle *q;
	for(i=0;i<p->en_new;i++){
		for(j=0,in=0;!in&&j<p->en_old;j++){
			if(p->new_list[i]==p->old_list[j]){
				p->old_list[j]=p->old_list[--p->en_old];
				in=1;
			}
		}
		if(!in){
			q=p->new_list[i];
			q->new_list[q->en_new++]=p;
		}
	}
	for(i=0;i<p->en_old;i++){
		q=p->old_list[i];
		for(j=0;j<q->en_new;j++)
			if(q->new_list[j]==p){
				q->new_list[j]=q->new_list[--q->en_new];
				break;
			}
	}
	return;
}
void delete_lists(particle *p){
	unsigned int i,j;
	particle *q;
	for(i=0;i<p->en_new;i++){
		q=p->new_list[i];
		for(j=0;j<q->en_new;j++)
			if(q->new_list[j]==p){
				q->new_list[j]=q->new_list[--q->en_new];
				break;
			}
	}
	return;
}
void swap_lists(particle *p,particle *w){
	unsigned int i,j;
	particle *q;
	for(i=0;i<p->en_new;i++){
		q=p->new_list[i];
		for(j=0;j<q->en_new;j++)
			if(q->new_list[j]==p){
				q->new_list[j]=w;
				break;
			}
	}
	return;
}

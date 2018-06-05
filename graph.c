#include <stdlib.h>
#include <stdio.h>
#include <utils.h>
#include "zargs.h"
#include "graph.h"
//Calculates connectivity graphs
int init_path(path *p,int size){
	int i;
	p->pleft=(path**)alloc(sizeof(path*)*size);
	p->pright=(path**)alloc(sizeof(path*)*size);
	p->e=NULL;
	p->length=0;
	p->nbranch=0;
	for(i=0;i<size;i++){
		p->pleft[i]=NULL;
		p->pright[i]=NULL;
	}
	return 0;
}
int clear_path(path *p,int size){
	int i;
	p->e=NULL;
	p->length=0;
	p->nbranch=0;
	for(i=0;i<size;i++){
		p->pleft[i]=NULL;
		p->pright[i]=NULL;
	}
	return 0;
}
int count_edges(matrix *m){
	int i;
	for(i=0,m->nedges=0;i<SQR(m->size);i++){
		if(*(*m->graph+i))m->nedges++;
	}
	return m->nedges;
}
int assign_edges(matrix *m){
	int i,j,k;
	unsigned x;
	edge *e;
	path *p,*q;
	k=count_edges(m);
	k=0,x=0;
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			p=&m->first[i][j];
			if(m->graph[i][j]){
				e=m->edges+k;
				e->i=i;
				e->j=j;
				e->num=k++;
				e->inverse=NULL;
				e->xor=0;

				p->e=e;
				p->nbranch=0;
				p->length=0;
			}
			else p->e=NULL;
		}
	}
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			p=&m->first[i][j];
			if(p->e){
				q=&m->first[j][i];
				p->e->inverse=q->e;
				if(j>i){
					p->e->xor=x;
					p->e->inverse->xor=x;
					x++;
				}
			}
		}
	}
	m->nxor=x;
	return 0;
}
int init_edges(matrix *m){
	int i,j,k;
	unsigned x;
	edge *e;
	path *p,*q;
	k=count_edges(m);
	m->edges=(edge*)alloc(sizeof(edge)*k);
	k=0,x=0;
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			p=&m->first[i][j];
			if(m->graph[i][j]){
				e=m->edges+k;
				e->i=i;
				e->j=j;
				e->num=k++;
				e->inverse=NULL;
				e->xor=0;

				p->e=e;
				p->nbranch=0;
				p->length=0;
			}
			else p->e=NULL;
		}
	}
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			p=&m->first[i][j];
			if(p->e){
				q=&m->first[j][i];
				p->e->inverse=q->e;
				if(j>i){
					p->e->xor=x;
					p->e->inverse->xor=x;
					x++;
				}
			}
		}
	}
	m->nxor=x;
	return 0;
}
branch *init_branch(int size){
	int i;
	branch *b=(branch*)alloc(sizeof(branch));
	b->taken=(int*)alloc(sizeof(int)*size);
	b->all=(int*)alloc(sizeof(int)*size);
	b->length=0;
	b->size=size;
	for(i=0;i<size;i++){
		*(b->taken+i)=0;
		*(b->all+i)=0;
	}
	return b;
}
stack_edge *init_stack_edge(int size){
	int i;
	stack_edge *s=(stack_edge*)alloc(sizeof(stack_edge)*size);
	for(i=0;i<size;i++){
		(s+i)->length=0;
		(s+i)->trace=(unsigned int*)alloc(sizeof(unsigned int)*TRACE);
		(s+i)->e=(edge**)alloc(sizeof(edge*)*size);
	}
	return s;
}
cycle *init_cycle(int size){
	int i;
	//int size2=size*size;
	int size2=CSET_ALLOC;
	cycle *c=(cycle*)alloc(sizeof(cycle)*size2);
	for(i=0;i<size2;i++){
		(c+i)->trace=(unsigned int*)alloc(sizeof(unsigned int)*TRACE);
		(c+i)->e=(edge**)alloc(sizeof(edge*)*size);
	}
	return c;
}
cycle_path_set *init_cycle_path_set(int size){
	int i,size2=MAX(size*size,CSET_ALLOC);
	cycle_path_set *cset=(cycle_path_set*)alloc(sizeof(cycle_path_set));
	cset->c=(cycle_path**)alloc(sizeof(cycle_path*)*size2);
	*cset->c=(cycle_path*)alloc(sizeof(cycle_path)*size2);
	cset->size_allocated=size2;
	cset->size=0;
	cset->head=NULL;
	for(i=0;i<size2;i++){
		*(cset->c+i)=*cset->c+i;
		(*(cset->c+i))->length=0;
	}
	return cset;
}
int clear_cycle_path_set(cycle_path_set *cset){
	int i;
	for(i=0;i<cset->size;i++){
		(*(cset->c+i))->length=0;
		(*(cset->c+i))->left=NULL;
		(*(cset->c+i))->right=NULL;
		(*(cset->c+i))->first=NULL;
		(*(cset->c+i))->second=NULL;
	}
	cset->size=0;
	cset->head=NULL;
	return 0;
}
matrix *init_matrix(int size){
	int i;
	matrix *m=(matrix*)alloc(sizeof(matrix));
	m->size=size;
	m->size_alloc=size;
	m->p=(node*)alloc(sizeof(node)*size);
	m->graph=(int**)alloc(sizeof(int*)*size);
	m->dist=(int**)alloc(sizeof(int*)*size);
	m->order=(int**)alloc(sizeof(int*)*size);
	m->first=(path**)alloc(sizeof(path*)*size);
	m->second=(path**)alloc(sizeof(path*)*size);
	//single alloc
	*m->graph=(int*)alloc(sizeof(int)*size*size);
	*m->dist=(int*)alloc(sizeof(int)*size*size);
	*m->order=(int*)alloc(sizeof(int)*size*size);

	*m->first=(path*)alloc(sizeof(path)*size*size);
	*m->second=(path*)alloc(sizeof(path)*size*size);
	for(i=0;i<size;i++){
		*(m->graph+i)=*m->graph+i*size;
		*(m->dist+i)=*m->dist+i*size;
		*(m->order+i)=*m->order+i*size;

		*(m->first+i)=*m->first+i*size;
		*(m->second+i)=*m->second+i*size;
	}
	for(i=0;i<size*size;i++){
		*(*m->graph+i)=0;
		*(*m->dist+i)=0;
		*(*m->order+i)=0;
		init_path(*m->first+i,PATH_SIZE);
		init_path(*m->second+i,PATH_SIZE);
	}
	m->c=init_cycle(size);
	m->b=init_branch(size);
	m->s=init_stack_edge(size);
	m->cset=init_cycle_path_set(size);
	return m;
}
void alloc_edges(matrix *m,int size){
	m->edges=(edge*)alloc(sizeof(edge)*size);
}
int graph2dist(matrix *m){
	unsigned int i,j;
	unsigned int n=m->size;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			m->dist[i][j]=m->graph[i][j];
			(&m->first[i][j])->length=m->graph[i][j];
			(&m->second[i][j])->length=0;
			(&m->second[i][j])->nbranch=0;
			if(i!=j){
				if(!m->dist[i][j]){
					m->dist[i][j]=DIST_MAX;
					(&m->first[i][j])->length=DIST_MAX;
				}
			}
			else{
				m->dist[i][j]=0;
				(&m->first[i][j])->length=0;
			}
		}
	}
	return 0;
}
//printg matrices
int print_graph(FILE *f,matrix *m){
	int i,j;
	fprintf(f,"%d\n",m->size);
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			fprintf(f,"%d",m->graph[i][j]);
		}
		fputc('\n',f);
	}
	return 0;
}
int print_dist(FILE *f,matrix *m){
	int i,j;
	fprintf(f,"%d\n",m->size);
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			if(m->dist[i][j]==DIST_MAX)fprintf(f,"*");
			else fprintf(f,"%d",m->dist[i][j]);
		}
		fputc('\n',f);
	}
	return 0;
}
//managing paths
int path_change(path *p,path *ik,path *kj){
	p->length=ik->length+kj->length;
	p->nbranch=1;
	*(p->pleft)=ik;
	*(p->pright)=kj;
	return 0;
}
int path_append(path *p,path *ik,path *kj){
	p->length=ik->length+kj->length;
	p->pleft[p->nbranch]=ik;
	p->pright[p->nbranch]=kj;
	p->nbranch++;
	return 0;
}
int path_copy(path *q,path *p){
	int i;
	for(i=0;i<p->nbranch;i++){
		q->pleft[i]=p->pleft[i];
		q->pright[i]=p->pright[i];
	}
	q->length=p->length;
	q->nbranch=p->nbranch;
	return 0;
}
int path_delete(path *p){
	p->length=0;
	p->nbranch=0;
	return 0;
}
void path_delete_matrix(matrix *m){
	int i;
	for(i=0;i<m->size*m->size;i++){
		path_delete(*m->first+i);
		path_delete(*m->second+i);
	}
}
void clear_all_path(matrix *m){
	int i;
	for(i=0;i<m->size_alloc*m->size_alloc;i++){
		*(*m->graph+i)=0;
		path_delete(*m->first+i);
		path_delete(*m->second+i);
	}
}
int dump_path(path *p){
	int i;
	path *s,*t;
	printf("p->nbranch %d p->length %d\n",p->nbranch,p->length);
	if(p->length){
		for(i=0;i<p->nbranch;i++){
			s=p->pleft[i];
			t=p->pright[i];
			printf("[%2d,%2d][%2d,%2d] {%2d}{%2d}(%d)\n",s->i,s->j,t->i,t->j,s->length,t->length,p->length);
		}
	}
	return 0;
}
int floydwarshall(matrix *m){
	int i,j,k;
	int n=m->size;
	int dikj;
	for(k=0;k<n;k++){
		for(i=0;i<n;i++){
			for(j=0;j<n;j++){
				if(i!=k&&k!=j){ //Excluding loops
					dikj=m->dist[i][k]+m->dist[k][j];
					if(m->dist[i][j]>dikj){
						if(m->dist[i][j]==(dikj+1)){
							path_copy(&m->second[i][j],&m->first[i][j]);
						}
						else{
							path_delete(&m->second[i][j]);
						}
						path_change(&m->first[i][j],&m->first[i][k],&m->first[k][j]);
						m->dist[i][j]=dikj;
					}
					else if(m->dist[i][j]==(dikj-1)){
						path_append(&m->second[i][j],&m->first[i][k],&m->first[k][j]);
					}
					else if(m->dist[i][j]==dikj){
						path_append(&m->first[i][j],&m->first[i][k],&m->first[k][j]);
					}
				}
			}
		}
	}
	return 0;
}
//branch management
int branch_clear(branch *b){
	int k;
	for(k=0,b->length=0;k<b->size;k++){
		b->taken[k]=0;
		b->all[k]=0;
	}
	return 0;
}
int branch_adjust(branch *b){
	int k=b->length,l;
	for(k=b->length-1;k>=0;k--){
		if(b->taken[k]<b->all[k]-1){
			b->taken[k]++;
			for(l=k+1;l<b->length;l++){
				b->taken[l]=0;
			}
			b->length=0;
			return k+1;
		}
	}
	b->length=0;
	return 0;
}
int store_branch(path *p,stack_edge *s,branch *b){
	int k;
	if(p){
		if(!p->e){
			k=b->taken[b->length];
			b->all[b->length]=p->nbranch;
			b->length++;
			store_branch(p->pleft[k],s,b);
			store_branch(p->pright[k],s,b);
		}
		else *(s->e+s->length++)=p->e;
		return 0;
	}
	else return 1;
}
//traces
void trace_null(unsigned int *t){
	unsigned int i;
	for(i=0;i<TRACE;i++){
		*(t+i)=0;
	}
}
void set_trace(unsigned int *t,unsigned int n){
	unsigned int word=n>>0x5;
	unsigned int s=n&0x1f;
	*(t+word)|=1<<s;		
}
int set_trace_check(unsigned int *t,unsigned int n){
	unsigned int word=n>>0x5;
	unsigned int s=n&0x1f;
	if(*(t+word)&1<<s)return 1;
	*(t+word)|=1<<s;		
	return 0;
}
void print_trace(unsigned int *t){
	int i,j;
	unsigned int s;
	for(j=TRACE-1;j>=0;j--){
		s=*(t+j);
		for(i=TRACE*8-1;i>=0;i--){
			putchar('0'+((s>>i)&0x1));
		}
	}
	putchar('\n');
}
int trace_stack(stack_edge *s){
	int i;
	trace_null(s->trace);
	for(i=0;i<s->length;i++){
		set_trace(s->trace,s->e[i]->xor);
	}
	return 0;
}
int trace_cycle(cycle *c){
	int i;
	trace_null(c->trace);
	for(i=0;i<c->length;i++){
		if(set_trace_check(c->trace,c->e[i]->xor))return 1;
	}
	return 0;
}
int cmp_trace(unsigned int *t,unsigned int *s){
	unsigned int i;
	for(i=0;i<TRACE;i++){
		if(*(t+i)^*(s+i))return 1; // they are not equal, they must not be equal
	}
	return 0;
}
int cmp_tstacks(stack_edge *s,int k){
	int i;
	trace_stack(s+k);
	for(i=0;i<k;i++){
		if(!cmp_trace((s+i)->trace,(s+k)->trace))return 1;
	}
	return 0;
}
int cmp_cycles(cycle *c,int k){
	int i;
	for(i=0;i<k;i++){
		if(!cmp_trace((c+i)->trace,(c+k)->trace))return 1;
	}
	return 0;
}
//edge stack manipulations
int cmp_stacks(stack_edge *s,stack_edge *t){
	int i;
	for(i=0;i<s->length;i++){
		if(s->e[i]!=t->e[i]){
			return 0;
		}
	}
	return 1;
}
int cmp_kstacks(stack_edge *s,int k){
	int i;
	stack_edge *t=s+k;
	for(i=0;i<k;i++){
		if(cmp_stacks(s+i,t))return 1;
	}
	return 0;
}
int store_stack(path *p,stack_edge *s,branch *b){
	int k=0;
	branch_clear(b);
	s->length=0;
	if(p->length){
		do{
			store_branch(p,s+k,b);
			if(p->length==(s+k)->length){
				if(!cmp_tstacks(s,k)){
					k++;
				}
			}
			(s+k)->length=0;
		}while(branch_adjust(b));
	}
	return k;
}
int print_stack_edge(FILE *f,stack_edge *s,int k){
	int i,j;
	stack_edge *t;
	for(i=0;i<k;i++){
		t=s+i;
		fprintf(f,"%d",(*t->e)->i);		
		for(j=0;j<t->length;j++){
			fprintf(f,"-%d",t->e[j]->j);		
		}
		fputc('\n',f);
	}
	return 0;
}
//cycle path set manipulation
void btree_insert(cycle_path **node,cycle_path *c){
	if(!*node)*node=c;
	else{
		if((*node)->length>c->length){
			btree_insert(&(*node)->left,c);
		}
		else btree_insert(&(*node)->right,c);
	}
}
int btree_inorder(cycle_path_set *cset,cycle_path *c){
	if(!c)return 1;
	else{
		btree_inorder(cset,c->left);
		*(cset->c+cset->size++)=c;
		btree_inorder(cset,c->right);
		return 0;
	}
}
int btree_sort(cycle_path_set *cset){
	cset->size=0;
	btree_inorder(cset,cset->head);
	return 0;
}
int cycle_path_set_append(cycle_path_set *cset,path *first,path *second){
	cycle_path *c;
	c=*(cset->c+cset->size);
	if(first->nbranch<2&&!second->length)return 1;
	else{
		if(!second->length)c->length=2*first->length;
		else c->length=2*first->length+1;
		c->first=first;
		c->second=second;
		if(!cset->size)cset->head=c;
		else btree_insert(&cset->head,c);
		cset->size++;
		return 0;
	}
}
int cycle_path_store(matrix *m){
	int i,j;
	for(i=0;i<m->size;i++){
		for(j=0;j<m->size;j++){
			if(m->dist[i][j]>0){
				cycle_path_set_append(m->cset,&m->first[i][j],&m->second[i][j]);
			}
		}
	}
	btree_sort(m->cset);
	return 0;
}
//finding cycles
int stack2cycle(stack_edge *s,cycle *c){
	int i;
	for(i=0;i<s->length;i++){
		c->e[c->length++]=s->e[i];
	}
	return 0;
}
int stack2cycle_inverse(stack_edge *s,cycle *c){
	int i;
	for(i=0;i<s->length;i++){
		c->e[c->length++]=s->e[s->length-i-1]->inverse;
	}
	return 0;
}
int print_cycle(FILE *f,cycle *c,int k){
	int i,j;
	cycle *t;
	for(i=0;i<k;i++){
		t=c+i;
		fprintf(f,"[%d] %d",t->length,(*t->e)->i);		
		for(j=0;j<t->length;j++){
			fprintf(f,"-%d",t->e[j]->j);		
		}
		fputc('\n',f);
		for(j=0;j<t->length;j++){ //modified
			fprintf(f,"-%d",t->e[j]->xor);		
		}
		fputc('\n',f);
		//print_trace((c+i)->trace);
	}
	return 0;
}
pcycle *alloc_pcycle(int size,int length){
	int i;
	pcycle *pc=(pcycle*)alloc(size*sizeof(pcycle));
	pcycle *pp;
	for(i=0;i<size;i++){
		pp=pc+i;
		pp->length=0;
		pp->p=(particle**)alloc(sizeof(particle*)*length);
	}
	return pc;
}
int print_pcycle(pcycle *pc,int k){
	int i,j;
	particle *q;
	for(i=0;i<k;i++){
		printf("k %d length %d  ",i,(pc+i)->length);
		for(j=0;j<(pc+i)->length;j++){
			q=*((pc+i)->p+j);
			printf("-%d-",q->n);
		}
		printf("\n");
	}
	return 0;
}
int print_pcycle_length_histogram(FILE *f,pcycle *pc,int k){
	int i;
	int h[MAX_CYCLE_LENGTH];
	for(i=0;i<MAX_CYCLE_LENGTH;i++){
		h[i]=0;
	}
	for(i=0;i<k;i++){
		h[(pc+i)->length]++;
	}
	for(i=0;i<MAX_CYCLE_LENGTH;i++){
		fprintf(f,"%d ",h[i]);
	}
	fprintf(f,"\n");
	return 0;
}
int print_plist_cycle(FILE *f,particle **p,cycle *c,int k){
	int i,j;
	int l;
	cycle *t;
	particle *q;
	for(i=0;i<k;i++){
		t=c+i;
		l=(*t->e)->i;
		q=*(p+l);
		fprintf(f,"[%d] %d",t->length,q->n);		
		for(j=0;j<t->length;j++){
			l=t->e[j]->j;
			q=*(p+l);
			fprintf(f,"-%d",q->n);		
		}
		fputc('\n',f);
	}
	return 0;
}
int cycle_odd(cycle *c,int *n,path *first,stack_edge *s,branch *b){
	int i,j,k;
	k=store_stack(first,s,b);
	for(i=0;i<k;i++){
		for(j=i+1;j<k;j++){
			(c+*n)->length=0;
			stack2cycle(s+i,c+*n);
			stack2cycle_inverse(s+j,c+*n);
			if(!trace_cycle(c+*n)){
				if(!cmp_cycles(c,*n)){
					(*n)++;
				}
			}
		}
	}
	return *n;
}
int cycle_even(cycle *c,int *n,path *first,path *second,stack_edge *s,branch *b){
	int i,j,k,l;
	k=store_stack(first,s,b);
	l=store_stack(second,s+k,b);
	for(i=0;i<k;i++){
		for(j=k;j<l+k;j++){
			(c+*n)->length=0;
			stack2cycle(s+i,c+*n);
			stack2cycle_inverse(s+j,c+*n);
			if(!trace_cycle(c+*n)){
				if(!cmp_cycles(c,*n)){
					(*n)++;
				}
			}
		}
	}
	return 0;
}
int cycle_from_set(matrix *m){
	int i;
	//int t=m->nxor-m->size-1;
	int t=m->nxor-m->size+2;
	m->ncycles=0;
	//printf("cset->size=%d/%d\n",m->cset->size,m->cset->size_allocated);
	for(i=0;(i<m->cset->size)&&(m->ncycles!=t);i++){
		if(m->cset->c[i]->length<MAX_CYCLE_LENGTH){
			if(m->cset->c[i]->length%2){
				cycle_even(m->c,&m->ncycles,m->cset->c[i]->first,m->cset->c[i]->second,m->s,m->b);
			}
			else{
				cycle_odd(m->c,&m->ncycles,m->cset->c[i]->first,m->s,m->b);
			}
		}
	}
	return 0;
}
void treewalk(particle *p,particle **t,int *n,int depth,int ndepth){
	unsigned i;
	particle *q;
	if(p->pass<ndepth-1){
		if(depth<ndepth){
			if(!p->pass){
				p->idx=(*n);
				*(t+*n)=p;
				(*n)++;
			}
			p->pass++;
			for(i=0;i<p->en_new;i++){
				q=p->new_list[i];
				treewalk(q,t,n,depth+1,ndepth);
			}
		}
	}
}
void clear_pass(particle **p,int n){
	int i;
	particle *q;
	for(i=0;i<n;i++){
		q=*(p+i);
		q->pass=0;
		q->idx=-1;
	}
}
void plist2matrix(matrix *m,particle **p,int n){
	int i,j;
	unsigned k;
	m->size=n;
	particle *q;
	for(i=0;i<n;i++){
		q=*(p+i);
		for(k=0;k<q->en_new;k++){
			j=q->new_list[k]->idx;
			if(j!=-1){
				m->graph[i][j]=1;
			}
		}
	}
}
int cycle2pcycle(particle **p,cycle *c,pcycle *pc){
	int i,l;
	particle *q;
	pc->length=c->length;
	for(i=0;i<c->length;i++){
		l=c->e[i]->j;
		q=*(p+l);
		*(pc->p+i)=q;
	}
	return 0;
}
int cmp_pcycle(pcycle *a,pcycle *b){
	int i,j;
	int status=0;
	if(a->length!=b->length)return 0;
	for(i=0;i<a->length;i++){
		status=0;
		for(j=0;j<b->length;j++){
			if((*(a->p+i))->n==(*(b->p+j))->n)status=1;
		}
		if(!status)return 0;
	}
	return 1;
}
int add_pcycle(pcycle *pcycles){
	int i;
	particle *q;
	pcycle **qcycles;
	for(i=0;i<pcycles->length;i++){
		q=*(pcycles->p+i);
		qcycles=q->pcycles;
		qcycles[q->npcycles]=pcycles;
		q->npcycles++;
	}
	return 0;
}
int c2pcycle(particle **p,cycle *c,pcycle *pc,int nc,int k){
	int i,j;
	int l=0;
	int status=0;
	for(i=0;i<k;i++){
		cycle2pcycle(p,c+i,pc+nc+l);
		status=0;
		for(j=0;j<nc;j++){
			if(cmp_pcycle(pc+j,pc+nc+l)){
				status=1;
				break;
			}
		}
		if(!status){
			add_pcycle(pc+nc+l);
			l++;
		}
	}
	return l;
}
void print_treewalk(particle *p,header *t){
	int i;
	int n=0;
	particle **l=t->matrix_list;
	particle *q;
	treewalk(p,l,&n,0,TREE_DEPTH);
	for(i=0;i<n;i++){
		q=*(l+i);
		printf("[%d] idx %d n %d (%p)\n",i,q->idx,q->n,p);
	}
	clear_pass(l,n);
}
int find_cycles(particle *p,particle **l,pcycle *pcycles,matrix *m,int ncycles){
	int n=0,k=0;
	treewalk(p,l,&n,0,TREE_DEPTH);
	while(n>MATRIX_ALLOC-1){
		clear_pass(l,n);
		n=0;
		k++;
		treewalk(p,l,&n,0,TREE_DEPTH-k);
	}
	m->size=n;
	plist2matrix(m,l,n);
	clear_pass(l,n);
	if(n>2){
		assign_edges(m);
		graph2dist(m);

		floydwarshall(m);
		cycle_path_store(m);
		cycle_from_set(m);
		
		k=c2pcycle(l,m->c,pcycles,ncycles,m->ncycles);
		
		clear_cycle_path_set(m->cset);
		clear_all_path(m);
	}
	return k;
}
int find_all_cycles(header *t){
	unsigned int i;
	t->ncycles=0;
	particle *p;
	species *s=t->specie;
	while(s){
		for(i=0;i<s->Nalloc;i++){
			p=(particle*)s->p+i;
			p->npcycles=0;
		}
		s=s->next;
	}
	s=t->specie;
	while(s){
		for(i=0;i<s->N;i++){
			p=(particle*)s->p+i;
			if(p->npcycles<p->npatch){ // to avoid covered particles already
				t->ncycles+=find_cycles(p,t->matrix_list,t->graph_pcycle,t->adj_matrix,t->ncycles);
			}
		}
		s=s->next;
	}
	//printf("[Number of cycles found %d]\n",t->ncycles);
	//print_pcycle(t->graph_pcycle,t->ncycles);
	//print_pcycle_length_histogram(stdout,t->graph_pcycle,t->ncycles);
	return t->ncycles;
}
void alloc_graph(header *t){
	t->ncycles=0;
	t->matrix_list=(particle**)alloc(sizeof(particle*)*t->Nalloc);
	t->adj_matrix=init_matrix(MATRIX_ALLOC);
	t->graph_pcycle=alloc_pcycle(t->Nalloc,PATH_SIZE);
	alloc_edges(t->adj_matrix,t->npatch_alloc);
}

#ifndef GRAPH_H
#define GRAPH_H
#define SQR(x) ((x)*(x))
#define DIST_MAX 99
#define MATRIX_ALLOC 64
#define MAX(a,b) ((a)>(b)?(a):(b))
#define CSET_ALLOC 65536*2
#if defined (TYPE)
	#define PATH_SIZE 32
	#define TREE_DEPTH 24
	#define TRACE 20
	#define MAX_CYCLE_LENGTH 24
#else
	#define PATH_SIZE 16
	#define TREE_DEPTH 12
	#define TRACE 16
	#define MAX_CYCLE_LENGTH 10
#endif
typedef struct branch{
	int length; //branching length
	int size;
	int *taken;
	int *all;
}branch;
typedef struct node{
	int flag;
	int pass;
	int orientation;
	int nbond; //number of bonds
	double *angle;
	void **c; //pointer to a cycle
	struct node **next;
}node;
typedef struct edge{
	int num;
	int flag; //xor for a particular edge
	int i,j;	
	unsigned int xor; //xor for a particular edge
	node *pi,*pj; //nodes pi and pj
	struct edge *inverse; //if an edge in opposite direction exists
}edge;
typedef struct stack_edge{
	int length;
	unsigned int *trace;
	edge **e;
}stack_edge;
typedef struct path{
	int length; //length of the path
	int nbranch; //number of branching
	int i,j;	
	edge *e;
	struct path **pleft,**pright;
}path;
typedef struct cycle_path{
	int length;
	path *first;
	path *second;
	struct cycle_path *left,*right;
}cycle_path;
typedef struct cycle_path_set{
	int size;
	unsigned int size_allocated;
	cycle_path **c;
	cycle_path *head;
}cycle_path_set;
typedef struct cycle{
	int flag;
	int length;
	int *na; //angle along the cycle
	unsigned int *trace;
	node **p;
	edge **e;
}cycle;
typedef struct matrix{
	int size,size_alloc;
	node *p;
	branch *b;
	stack_edge *s;
	cycle_path_set *cset;
	int ncycles;
	cycle *c;
	int nedges,nxor;
	edge *edges;
	int **graph; //graph, adjency matrix
	int **dist; //distance matrix
	int **order;
	path **first;
	path **second;
}matrix;
typedef struct pcycle{
	int length;
	particle **p;
}pcycle;
extern void alloc_graph(header *t);
extern int find_all_cycles(header *t);
extern int print_pcycle_length_histogram(FILE *f,pcycle *pc,int k);
#endif

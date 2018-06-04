#ifndef LISTS_H
#define LISTS_H
#define SWAP(a,b) {unsigned int (c)=(a);(a)=(b);(b)=(c);}
#define LIST_SWAP(a,b) {particle **c=(a);(a)=(b);(b)=c;}
extern void list_swap(particle *p);
extern void adjust_lists(particle *p);
extern void delete_lists(particle *p);
extern void swap_lists(particle *p,particle *w);
#endif

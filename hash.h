#ifndef HASH_H
#define HASH_H
#define M 1048576//1024x1024 hash table size (<<10)
#define MAX(a,b) ((a)>(b)?(a):(b))
extern void hash1(header *t);
extern void hash_dirs(header *t);
extern void hash_alloc(header *t);
extern void hash_lists(header *t);
extern void set_hash(header *t);

extern void hash_insert(particle *p,__m128d h1,hash_table *table);
extern void hash_reinsert(particle *p,__m128d h1,hash_table *table);
extern void hash_delete(particle *p);

extern unsigned hash(__m128d x,__m128d h1);
#endif

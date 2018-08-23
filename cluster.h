#ifndef CLUSTER_H
#define CLUSTER_H
extern int new_cluster(compound_particle *c,header *t);
extern int add2cluster(compound_particle *c,header *t);
extern int print_clusters(header *t);
extern void set_all_particle_color(header *t);
extern void color_cluster(cluster *cc,float *color);

extern int cluster_check_particle(particle *p,header *t);
extern int find_new_cluster(compound_particle *c,header *t);
extern int clusters_reset(header *t);
extern int find_all_clusters(header *t);

extern void color_all_clusters(header *t);
extern void avg_max_cluster_size(header *t);

#endif

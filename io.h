#ifndef IO_H
#define IO_H
extern void print_log(FILE *f,header *t,long long int i,double time,double frac[4]);
extern void open_files(header *t);
extern void close_files(header *t);
extern void flush_files(header *t);
extern void write_files(header *t);
#endif

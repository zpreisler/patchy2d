#include <stdlib.h>
#include <stdio.h>
#define OPT_MIN 0.0000001
void gfrac(int *acc[2],double *frac,int n){
	int i,tot;
	for(i=0;i<n;i++){
		tot=(acc[i][0]+acc[i][1]);
		if(!tot)tot=1;
		*(frac+i)=(double)acc[i][0]/tot;
		acc[i][0]=acc[i][1]=0;
	}
	return;
}
void optimize(double **max,double *frac,int n){
	int i;
	for(i=0;i<n;i++){
		if(*(max+i)){
			*max[i]*=*(frac+i)*2.0;
			if(*max[i]>1.0)*max[i]=1.0;
			else if(*max[i]<OPT_MIN)*max[i]=OPT_MIN;
		}
	}	
	return;
}

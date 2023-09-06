//modello macroscopico come Bagnoli_deg1 ma con degenerazione approssimata

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void degenerazione(int *g, int L, int *D);
double normalize(double *p, int L);
void compute_fitness(double *H, double *H0, double *p, int L);

double b, r, J, R, a;

int main () {
	int L;  // dimensione del fenotipo
	int *g;   // degenerazione
	double *c;  //"coefficienti di degenerazione"
	int D;  //degenerazione massima
	double *p, *p1;  // popolazione
	double *H0;      // fitness statica
	double *H;       // fitness 
	int i, j;
	int t, T;
	double mu;
  
	L=32;
	b=0.04;
	R=3.2;
	a=2;
	r=0.96;
	J=0.6;
	T=1000;
	mu=0.003;

	p = calloc(L, sizeof(double));
	p1 = calloc(L, sizeof(double));
	g = calloc(L, sizeof(int));
	c = calloc(2*L, sizeof(double));
	H0 = calloc(L, sizeof(double));
	H = calloc(L, sizeof(double));
  
	degenerazione(g, L, &D);
	
	c[0]=0;  //calcolo coefficienti
	c[1]=1;
	for(i=2; i<L+1; i+=2){
		j=(int)(0.5*(double)i);
		c[i]=0.5*(double)g[j]/(double)D;
		c[i+1]=1-c[i];
	}
	for(i=L+2; i<2*L-2; i+=2){
		j=(int)(0.5*(double)i);
		c[i+1]=0.5*(double)g[j]/(double)D;
		c[i]=1-c[i+1];
	}
	c[2*L-2]=1;
	c[2*L-1]=0;
	
	for (i=0; i<L; i++) {
		p[i] = 1;
		H0[i] = b *(1-fabs(i-L/2.)/r - 1 / (1 + fabs(i-L/2.)/r));
	} 
	normalize(p, L);
	
	for (t=0; t<T; t++ ) {
		compute_fitness(H, H0, p, L);
		// survival
		for (i=0; i<L; i++) {
			p[i] *= exp(H[i]);
		}
		//mutations
		p1[0]=(1-c[1]*mu)*p[0]+c[2]*mu*p[1];  //mu=probabilità totale di avere mutazione
		for (i=2; i<2*L-2; i+=2){
			j=(int)(0.5*(double)i);
			p1[j]=(1-(c[i]+c[i+1])*mu)*p[j]+c[i+2]*mu*p[j+1]+c[i-1]*mu*p[j-1];
		}
		p1[L-1]=(1-c[2*L-2]*mu)*p[L-1]+c[2*L-3]*mu*p[L-2];
		
		for (i=0; i<L; i++) {
			p[i] = p1[i];
		}
		
		normalize(p, L);
    
	}
  
	FILE *dati;
	dati=fopen("B.dat", "w");
	fprintf(dati, "#fenotipo distribuzione:\n");
	for(i=0; i<L; i++){
		fprintf(dati, "%d %f\n", i, *(p+i));
	}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness_statica valore:\n");
	//for(i=0; i<L; i++){
		//fprintf(dati, "%d %f\n", i, *(H0+i));
	//}
	fprintf(dati, "\n");
	fprintf(dati, "#fitness_tout_court valore:\n");
	for(i=0; i<L; i++){
		fprintf(dati, "%d %f\n", i, *(H+i));
	}
	fprintf(dati, "\n");
	//fprintf(dati, "#degenerazione valore:\n");
	//for(i=0; i<L; i++){
		//fprintf(dati, "%d %d\n", i, *(g+i));
	//}
	fprintf(dati, "\n");
	fprintf(dati, "#coeff_degenerazione valore:\n");
	for(i=0; i<2*L-2; i++){
		fprintf(dati, "%d %f\n", i, *(c+i));
	}
	fclose(dati);
  
	return 0;
}

double normalize(double *p, int L) {
	int i; 
	double norm = 0;
  
	for (i=0; i<L; i++) {
		norm += p[i];
	}
	for (i=0; i<L; i++) {
		p[i]/=norm;
	}
	return(norm);
}

void compute_fitness(double *H, double *H0, double *p, int L){
	int i, j;
	double tmp;
  
	for (i=0; i<L; i++) {
		// static fitness
		tmp = H0[i];
		for (j=0; j<L; j++) {
			tmp -= J * exp(-pow(fabs((j-i)/R),a)/a)* p[j];
		}
		H[i] = tmp;
	}
}

void degenerazione(int *g, int L, int *D){
	
	*D=601080390;
	double K=6;
	int h=2;
	int i;
	int app;
	
	for(i=0; i<L; i++){
		*(g+i)=*D*exp(-pow(fabs((double)(i-L/2.))/K, h)/h);
		app=(int)(*(g+i));
		if(*(g+i)-(double)app>=0.5){
			*(g+i)=app+1;
		}else{
			*(g+i)=app;
		}
		if(*(g+i)==0){
			*(g+i)=1;
		}
	}
}






void degenerazioneALT(int *g, int L){  //per L pari
	
	int N=10000;  //scelta della dimensione del campione su cui stimare la degenerazione
	int D=3432;  //degenerazione massima =L!/((L/2)!)^2
	int i, j;
	int *v, *n;
	int nMAX=0;
	int provv, app;
	
	v=calloc(N, sizeof(int));
	n=calloc(L+1, sizeof(int));
	
	for(i=0; i<N; i++){
		*(v+i)=0;
		for(j=0; j<L; j++){
			(*(v+i))+=(int)(drand48()*2.0);
		}
	}
	
	for(i=0; i<L+1; i++){
		for(j=0; j<N; j++){
			if(*(v+j)==i){
				(*(n+i))++;
			}
		}
	}
	
	for(i=0; i<L+1; i++){
		if(*(n+i)>nMAX){
			nMAX=*(n+i);
		}
	}
	
	for(i=0; i<L+1; i++){
		provv=(int)((double)D/nMAX**(n+i));
		app=(int)provv;
		if((provv-(double)app>=0.5)||(app==0)){
			*(g+i)=app+1;
		}else{
			*(g+i)=app;
		}
	}
}


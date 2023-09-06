//modello macroscopico per riproduzione asessuata senza degenerazione

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double normalize(double *p, int L);
void compute_fitness(double *H, double *H0, double *p, int L);

double b, r, J, R, a;

int main () {
	int L; // dimensione del genoma e quindi del fenotipo
	int *g;      // degenerazione
	double *p, *p1;  // popolazione
	double *H0;      // fitness statica
	double *H;       // fitness 
	int i;
	int t, T;                                //b=0.04, J=0.6, R=10, r=3 (per L=100)
	double mu;                               //b=0.04, J=0.6, R=5, r=3, mu=0.0001 per 5 quasispecie (not really)
  
	L=100;
	b=0;
	R=25;
	a=5;
	r=3;
	J=0.6;
	T=1000;
	mu=0.001;

	p=calloc(L, sizeof(double));
	p1=calloc(L, sizeof(double));
	g=calloc(L, sizeof(double));
	H0=calloc(L, sizeof(double));
	H=calloc(L, sizeof(double));
  
  // inizializzazione
	for(i=0; i<L; i++){
		g[i]=1; // senza degenerazione
		p[i]=1;
		//p[i]=0;
		H0[i]=b*(1-fabs(i-L/2.)/r-1/(1+fabs(i-L/2.)/r));
		//H0[i]=exp(-pow((i-L/2.)/r,2)/2)/sqrt(2*M_PI*r*r);
	}
	//p[0]=1;
	//p[L/2]=1;
	normalize(p, L);
	
	for(t=0; t<T; t++){
		compute_fitness(H, H0, p, L);
		// survival
		for(i=0; i<L; i++){
			p[i] *= exp(H[i]);
		}
		// mutations
		p1[0]=(1-mu)*p[0]+mu*p[1];
		for(i=1; i<L-1; i++){
			p1[i]=(1-mu)*p[i]+mu*p[i+1]/2+mu*p[i-1]/2;
		}
		p1[L-1]=(1-mu)*p[L-1]+mu*p[L-2];
		
		for(i=0; i<L; i++){
			p[i]=p1[i];
		}
		normalize(p, L);
    
	}
  
	FILE *dati;
	dati=fopen("B.dat", "w");
	fprintf(dati, "#fenotipo distribuzione fitness_statica fitness_tout_court:\n");
	for(i=0; i<L; i++){
		fprintf(dati, "%d %f  %f  %f\n", i, *(p+i), *(H0+i), *(H+i));
	}
	fclose(dati);
  
	return 0;
}

double normalize(double *p, int L){
	int i; 
	double norm=0;
  
	for(i=0; i<L; i++){
		norm+=p[i];
	}
	for (i=0; i<L; i++) {
		p[i]/=norm;
	}
	return(norm);
}

void compute_fitness(double *H, double *H0, double *p, int L){
	int i, j;
	double tmp;
  
	for(i=0; i<L; i++){
		// static fitness
		tmp=H0[i];
		for(j=0; j<L; j++){
			tmp-=J*exp(-pow(fabs((j-i)/R),a)/a)*p[j];
		}
		H[i]=tmp;
	}
}


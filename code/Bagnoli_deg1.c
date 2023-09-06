#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double normalize(double *p, int L);
void compute_fitness(double *H, double *H0, double *p, int L);

double b, r, J, R, a;

int main(){
	int L;  // dimensione del genoma
	double *p, *p1;  // popolazione
	double *H0;      // fitness statica
	double *H;       // fitness 
	int i;
	int t, T;
	double mu;
  
	L=100;
	b=0.04;
	R=20;
	a=2;
	r=6;
	J=0.6;
	T=1000;
	mu=0.01;

	p=calloc(L+1, sizeof(double));
	p1=calloc(L+1, sizeof(double));
	H0=calloc(L+1, sizeof(double));
	H=calloc(L+1, sizeof(double));
	
	// inizializzazione
	for (i=0; i<L+1; i++) {
		p[i]=1;
		//p[i]=0;
		//H0[i]=b*(1-fabs(i-L/2.)/r-1/(1+fabs(i-L/2.)/r));
		H0[i] = exp(-pow((i-L/2.)/r,2)/2)/sqrt(2*M_PI*r*r);
	} 
	//p[0]=1;
	//p[L/2]=1;
	normalize(p, L);
	
	for (t=0; t<T; t++) {
		compute_fitness(H, H0, p, L);
		// survival
		for (i=0; i<L+1; i++) {
			p[i]*=exp(H[i]);
		}
		// mutations
		p1[0]=(1-mu)*p[0]+mu*p[1]/L;
		for (i=0; i<L+1; i++) {
			p1[i]=(1-mu)*p[i]+(i+1)*mu*p[i+1]/L+(L-i+1)*mu*p[i-1]/L;
		}
		p1[L]=(1-mu)*p[L]+mu*p[L-1]/L;
		
		for (i=0; i<L+1; i++) {
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
  
	for(i=0; i<L+1; i++){
		norm+=p[i];
	}
	for(i=0; i<L+1; i++){
		p[i]/=norm;
	}
	return(norm);
}

void compute_fitness(double *H, double *H0, double *p, int L){
	int i, j;
	double tmp;
  
	for(i=0; i<L+1; i++){
		// static fitness
		tmp=H0[i];
		for(j=0; j<L+1; j++){
			tmp-=J*exp(-pow(fabs((j-i)/R),a)/a)*p[j];
		}
		H[i]=tmp;
	}
}

int count_maxima(double *p, int L){
	int i, n;
  n=0;
  
  if(p[0]>p[1]) n++;
  for(i=1; i<L-2; i++){
    if(p[i]>p[i-1]&&p[i]>p[i+1]) n++;
  }
  if(p[L-1]>p[L-2]) n++;
  
  return(n);
}

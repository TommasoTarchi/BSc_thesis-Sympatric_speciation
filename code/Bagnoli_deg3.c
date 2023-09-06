//codice corretto da Bagnoli

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double normalize(double *p, int L);

void compute_fitness(double *H, double *H0, double *p, int L);

int count_maxima(double *p, int L);

double J, R, r, a;
 

int main (int argc, char ** argv) {
	int L = 100; // dimensione del genoma e quindi del fenotipo
	double *p, *p1;  // popolazione
	double *H0;      // fitness statica
	double *H;       // fitness
	int i;
	int t, T;
	double mu;
	double b=0.04;
  
	L = 100;
	T = 2000;
	mu = 0.01;
	a=2.; // kernel gaussiano
  
	if (argc < 4) {
		printf("usage: %s <R> <J> <r>\n", argv[0]);
		exit(1);
	}
	r = atof(argv[1]);
	J = atof(argv[2]);
	R = atof(argv[3]);
 
	p = calloc(L, sizeof(double));
	p1 = calloc(L, sizeof(double));
	H0 = calloc(L, sizeof(double));
	H = calloc(L, sizeof(double));
  
	// inizializzazione
	for (i=0; i<L; i++) {
		p[i] = 1;
		//H0[i] = b *(1-fabs(i-L/2.)/r - 1 / (1 + fabs(i-L/2.)/r));
		//H0[i] = exp(-pow((i-L/2.)/r,2)/2)/sqrt(2*M_PI*r*r);
	} 
	normalize(p, L);
  
	for (t=0; t<T; t++ ) {
		compute_fitness(H, H0, p, L);
		// survival
		for (i=0; i<L; i++) {
			p[i] *= exp(H[i]);
		}
    // mutations
    p1[0] = (1-mu) * p[0] + mu * p[1];
    for (i=1; i<L-1; i++) {
      // evoluzione senza degenerazione
      //p1[i] = (1-2*mu) * p[i] + mu * p[i+1] + mu * p[i-1];
      // evoluzione con degenerazione
      p1[i] = (1-mu) * p[i] + mu/L * (p[i+1]*(L-i) + p[i-1]*i); // **** degeneracy
    }
    p1[L-1] = (1-mu) * p[L-1] + mu * p[L-2];
    
    //p1[0]=(1-L*mu)*p[0]+mu*p[1];  //mutazione come in Bagnoli_deg1
		//for (i=0; i<L+1; i++) {
			//p1[i]=(1-L*mu)*p[i]+(i+1)*mu*p[i+1]+(L-i+1)*mu*p[i-1];
		//}
	//p1[L]=(1-L*mu)*p[L]+mu*p[L-1];

    
    for (i=0; i<L; i++) {
		p[i] = p1[i];
		}
		normalize(p, L);
	}
    
	// visualization
	FILE *dati;                            //STAMPA LA ROBA SBAGLIATA
	dati=fopen("B.dat", "w");
	fprintf(dati, "#fenotipo distribuzione fitness_statica fitness_tout_court:\n");
	for(i=0; i<L; i++){
		fprintf(dati, "%d %f  %f  %f\n", i, *(p+i), *(H0+i), *(H+i));
	}
	fclose(dati);
	
	return 0;
}

int count_maxima(double *p, int L) {
	int i, n;
	n=0;
  
	if (p[0]>p[1]) n++;
	for (i=1; i<L-2; i++) {
		if (p[i] > p[i-1] && p[i] > p[i+1]) n++;
	}
	if (p[L-1] > p[L-2]) n++;
	
	return(n);
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


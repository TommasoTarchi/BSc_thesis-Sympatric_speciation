//modello microscopico per riproduzione asessuata

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void genotipi_random(int **g);
void fenotipi(int **g, int *f);
void distribuzione(double *p, int *f);
void fitness(double *H0, double *H, double *A, double *Am, double *p, int **g, int *f, int t);
void riproduzione_mutazione(int *gG, int *gF);
void sopravvivenza(double A, double Am, int *k);

int L, N, D;
double T, mu, b, r, J, a, R;

int main(){
	
	L=20;  //lunghezza genoma
	N=1000;  //numero individui
	T=50000;  //generazioni
	b=0.04;
	r=0.6;
	J=0.6;
	a=2;
	R=2;
	mu=0.001;
	
	int i, j, k;
	int t;
	int **g;
	int *f;
	int **gF;
	int fF;
	double *p;
	double *H0;
	double *H;
	double *A;
	double Am=0;
	int *mark;  //marker per la riproduzione
	int markTOT;
	int v;
	
	srand48(time(NULL));
	
	g=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		*(g+i)=calloc(L, sizeof(int));
	}
	f=calloc(N, sizeof(int));
	gF=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		*(gF+i)=calloc(L, sizeof(int));
	}
	p=calloc(L+1, sizeof(double));
	H0=calloc(L+1, sizeof(double));
	H=calloc(L+1, sizeof(double));
	A=calloc(L+1, sizeof(double));
	mark=calloc(N, sizeof(int));
	
	genotipi_random(g);  //inizializza casualmente i genotipi
	fenotipi(g, f);  //calcola i fenotipi a partire dai genotipi
	
	for(i=0; i<L+1; i++){  //calcolo fitness statica
		H0[i]=b*(1-fabs(i-L/2.)/r-1/(1+fabs(i-L/2.)/r));
		//H0[i] = exp(-pow((i-L/2.)/r,2)/2)/sqrt(2*M_PI*r*r);
	}
	
	for(t=0; t<T; t++){
		
		for(i=0; i<N; i++){
			mark[i]=0;
		}
		distribuzione(p, f);  //calcola la distribuzione della popolazione nello spazio fenotipico
		fitness(H0, H, A, &Am, p, g, f, t);  //calcola fitness
		
		k=0;
		while(k<N){  //riproduzione-mutazione-sopravvivenza
			
			markTOT=0;  //reset dei marker
			for(i=0; i<N; i++){
				markTOT+=mark[i];
			}
			if(markTOT==N){
				for(i=0; i<N; i++){
					mark[i]=0;
				}
			}
			
			v=(int)(drand48()*N);  //scelta del genitore
			while(mark[v]==1){
				v=(int)(drand48()*N);
			}
			mark[v]=1;
			
			riproduzione_mutazione(*(g+v), *(gF+k));
			fF=0;
			for(j=0; j<L; j++){
				fF+=*(*(gF+k)+j);
			}
			sopravvivenza(*(A+fF), Am, &k);
		}
		
		for(i=0; i<N; i++){
			for(j=0; j<L; j++){
				*(*(g+i)+j)=*(*(gF+i)+j);
			}
		}
		fenotipi(g, f);
	}
	
	distribuzione(p, f);
	fitness(H0, H, A, &Am, p, g, f, t);
	
	FILE *dati;
	dati=fopen("ra.dat", "w");
	fprintf(dati, "#fenotipo distribuzione fitness_statica fitness_tout_court:\n");
	for(i=0; i<L+1; i++){
		fprintf(dati, "%d  %f  %f  %f\n", i, *(p+i), *(H0+i), *(H+i));
	}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness valore:\n");
	//for(i=0; i<L+1; i++){
		//fprintf(dati, "%d %f\n", i, *(A+i));
	//}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness_media:\n");
	//fprintf(dati, "%f\n\n", Am);
	//fprintf(dati, "#tasso_mutazione:\n");
	//fprintf(dati, "%f\n", mu);
	fclose(dati);
	
	return 0;
}

void genotipi_random(int **g){
	
	int i, j;
	
	for(i=0; i<N; i++){
		for(j=0; j<L; j++){
			*(*(g+i)+j)=(int)(drand48()*2.0);
		}
	}
}

void fenotipi(int **g, int *f){
	
	int i, j;
	
	for(i=0; i<N; i++){
		*(f+i)=0;
		for(j=0; j<L; j++){
			*(f+i)+=*(*(g+i)+j);
		}
	}
}

void distribuzione(double *p, int *f){
	
	int i, j;
	
	for(i=0; i<L+1; i++){
		*(p+i)=0;
		for(j=0; j<N; j++){
			if(*(f+j)==i){
				*(p+i)+=1.0/(double)N;
			}
		}
	}
}

void fitness(double *H0, double *H, double *A, double *Am, double *p, int **g, int *f, int t){
	
	int i, j;
	
	for(i=0; i<L+1; i++){
		H[i]=H0[i];
		for(j=0; j<L+1; j++){
			*(H+i)-=J*exp(-pow(fabs((j-i)/R),a)/a)*p[j];
		}
		*(A+i)=exp(*(H+i));
	}
	*Am=0;
	for(i=0; i<L+1; i++){
		*Am+=*(A+i)**(p+i);
	}
}

void riproduzione_mutazione(int *gG, int *gF){
	
	int i;
	int w=(int)(drand48()*L);
	double m=drand48();
	
	for(i=0; i<L; i++){
		*(gF+i)=*(gG+i);
	}
	*(gF+w)=(m<mu)*(1-*(gF+w))+(m>mu)**(gF+w);
}

void sopravvivenza(double A, double Am, int *k){
	
	if(drand48()<A/Am){
		(*k)++;
	}
}

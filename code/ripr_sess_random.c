#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void genotipi_random(int **g);
void fenotipi(int **g, int *f);
void distribuzione(double *p, int *f);
void fitness(double *H0, double *H, double *A, double *Am, double *p, int **g, int *f, int t);
void riproduzione_mutazione(int *g1, int *g2, int *gF);
void sopravvivenza(double A, double Am, int *k);

int L, N, D;
double T, mu, b, GAMMA, J, a, R, mu0, muINF, d, tau;

int main(){
	
	L=14;  //lunghezza genoma
	N=100;  //numero individui
	T=100000;  //generazioni
	D=14;  //assortatività
	b=1;
	GAMMA=10;
	J=16;
	a=2;
	R=7;
	mu0=0.1;
	muINF=0.000001;
	d=1000;
	tau=20000;
	
	int i, j;
	int t;
	int **g;
	int *f;
	int **gF;
	double *p;
	double *H0;
	double *H;
	double *A;
	double Am=0;
	
	int k;  //per la riproduzione
	int w, v;
	int fv, fw, fF;
	
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
	
	genotipi_random(g);  //inizializza casualmente i genotipi
	fenotipi(g, f);  //calcola i fenotipi a partire dai genotipi
	
	for(i=0; i<L+1; i++){
		*(H0+i)=exp(-1.0/b*pow(i/GAMMA, b));
	}
	
	for(t=0; t<T; t++){
		
		distribuzione(p, f);  //calcola la distribuzione della popolazione nello spazio fenotipico
		fitness(H0, H, A, &Am, p, g, f, t);  //calcola fitness
		mu=(mu0-muINF)/2*(1-tanh((t-tau)/d))+muINF;

		k=0;
		while(k<N){  //riproduzione-mutazione-selezione
			
			v=0;  //scelta dei genitori
			w=0;
			fv=0;
			fw=D+1;
			while((v==w)||(fabs((double)fv-(double)fw)>D)){
				w=(int)(drand48()*N);
				v=(int)(drand48()*N);
				fv=0;
				for(i=0; i<L; i++){
					fv+=*(*(g+v)+i);
				}
				fw=0;
				for(i=0; i<L; i++){
					fw+=*(*(g+w)+i);
				}
			}
			
			riproduzione_mutazione(*(g+v), *(g+w), *(gF+k));
			fF=0;
			for(i=0; i<L; i++){
				fF+=*(*(gF+k)+i);
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
	mu=(mu0-muINF)/2*(1-tanh((t-tau)/d))+muINF;
	
	FILE *dati;
	dati=fopen("rs.dat", "w");
	fprintf(dati, "#fenotipo distribuzione:\n");
	for(i=0; i<L+1; i++){
		fprintf(dati, "%d %f\n", i, *(p+i));
	}
	fprintf(dati, "\n");
	fprintf(dati, "#fitness_statica valore:\n");
	for(i=0; i<L+1; i++){
		fprintf(dati, "%d %f\n", i, *(H0+i));
	}
	fprintf(dati, "\n");
	//fprintf(dati, "#fitness_tout_court valore:\n");
	//for(i=0; i<L+1; i++){
		//fprintf(dati, "%d %f\n", i, *(H+i));
	//}
	//fprintf(dati, "\n");
	fprintf(dati, "#fitness valore:\n");
	for(i=0; i<L+1; i++){
		fprintf(dati, "%d %f\n", i, *(A+i));
	}
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
		*(H+i)=*(H0+i);
		for(j=0; j<L+1; j++){
			*(H+i)-=J*exp(-1.0/a*fabs(pow((i-j)/R, a)))**(p+j);
		}
		*(A+i)=exp(*(H+i));
	}
	*Am=0;
	for(i=0; i<L+1; i++){
		*Am+=*(A+i)**(p+i);
	}
}

void riproduzione_mutazione(int *g1, int *g2, int *gF){
	
	int i;
	double k=drand48();
	int w=(int)(drand48()*L);
	double m=drand48();
	
	for(i=0; i<L; i++){
		*(gF+i)=(k<0.5)**(g1+i)+(k>0.5)**(g2+i);
	}
	*(gF+w)=(m<mu)*(1-*(gF+w))+(m>mu)**(gF+w);
}

void sopravvivenza(double A, double Am, int *k){
	
	if(drand48()<A/Am){
		(*k)++;
	}
}

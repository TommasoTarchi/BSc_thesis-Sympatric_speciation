//Dieckman e Doebeli per organimsi diploidi con solo crattere ecologico e random mating in presenza di mutazione
//riproduzione sessuata ma ermafrodita
//meiosi senza crossing over
//monogamia

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void fenotipi(int **g, int *f);
void distribuzione(double *p, int *f);
void fitness(double *H0, double *H, double *A, double *Am, double *p, int **g, int *f, int t);
void meiosi(int *g, int *gam);
void riproduzione_mutazione(int *gam1, int *gam2, int *gF);
void sopravvivenza(double A, double Am, int *k);

int N, L, T;
double b, r, J, R, a, mu;

int main(){
	
	N=100;  //numero individui
	L=100; //lunghezza di ogni cromosoma
	T=10000;  //generazioni
	b=50;
	r=10;
	J=70;
	a=2;
	R=0.1;
	mu=0.001;
	
	int i, j, k, v, w;
	int t;
	int **g;
	int *f;
	int *gam1, *gam2;
	int **gF;
	int fF;
	double *p;
	double *H0;
	double *H;
	double *A;
	double Am;
	int *mark;  //marker per l'accoppiamento
	int markTOT;
		
	srand48(time(NULL));
	
	g=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		g[i]=calloc(2*L, sizeof(int));
	}
	f=calloc(N, sizeof(int));
	gam1=calloc(L, sizeof(int));
	gam2=calloc(L, sizeof(int));
	gF=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gF[i]=calloc(2*L, sizeof(int));
	}
	p=calloc(L+1, sizeof(double));
	H0=calloc(L+1, sizeof(double));
	H=calloc(L+1, sizeof(double));
	A=calloc(L+1, sizeof(double));
	mark=calloc(N, sizeof(int));
	
	for(i=0; i<L+1; i++){  //calcolo fitness
		//H0[i]=b*(1-fabs(i-L/2.)/r-1/(1+fabs(i-L/2.)/r));
		H0[i]=b*exp(-pow(fabs((i-L/2)/r),a)/a);
	}
	fenotipi(g, f);
	
	for(t=0; t<T; t++){
		
		for(i=0; i<N; i++){
			mark[i]=0;
		}
		distribuzione(p, f);
		fitness(H0, H, A, &Am, p, g, f, t);
		
		k=0;
		while(k<N){
			
			markTOT=0;  //reset dei marker
			for(i=0; i<N; i++){
				markTOT+=mark[i];
			}
			if(markTOT==N){
				for(i=0; i<N; i++){
					mark[i]=0;
				}
			}
			
			v=(int)(drand48()*N);  //scelta dei genitori
			while(mark[v]==1){
				v=(int)(drand48()*N);
			}
			mark[v]=1;
			w=v;
			while(mark[w]==1){
				w=(int)(drand48()*N);
			}
			mark[w]=1;
			
			meiosi(g[v], gam1);  //riproduzione
			meiosi(g[w], gam2);
			riproduzione_mutazione(gam1, gam2, gF[k]);
			
			fF=0;  //sopravvivenza
			for(i=0; i<L; i++){
				if((*(gF[k]+i)==1)||(*(gF[k]+i+L)==1)){
					fF++;
				}
			}
			sopravvivenza(A[fF], Am, &k);
		}
		
		for(i=0; i<N; i++){
			for(j=0; j<2*L; j++){
				*(*(g+i)+j)=*(*(gF+i)+j);
			}
		}
		
		fenotipi(g, f);
	}
	
	distribuzione(p, f);
	fitness(H0, H, A, &Am, p, g, f, t);
	
	FILE *dati;
	dati=fopen("DD.dat", "w");
	fprintf(dati, "#fenotipo distribuzione:\n");
	for(i=0; i<L+1; i++){
		fprintf(dati, "%d %f\n", i, *(p+i));
	}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness_statica valore:\n");
	//for(i=0; i<L+1; i++){
		//fprintf(dati, "%d %f\n", i, *(H0+i));
	//}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness_tout_court valore:\n");
	//for(i=0; i<L+1; i++){
		//fprintf(dati, "%d %f\n", i, *(H+i));
	//}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness valore:\n");
	//for(i=0; i<L+1; i++){
		//fprintf(dati, "%d %f\n", i, *(A+i));
	//}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fitness_media:\n");
	//fprintf(dati, "%f\n\n", Am);
	//fprintf(dati, "%f\n", mu);
	fclose(dati);
	
	
	return 0;
}

void fenotipi(int **g, int *f){
	
	int i, j;
	
	for(i=0; i<N; i++){
		f[i]=0;
		for(j=0; j<L; j++){
			if((*(g[i]+j)==1)||(*(g[i]+j+L)==1)){
				f[i]++;
			}
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
			*(H+i)-=J*exp(-pow(fabs((j-i)/R),a)/a)*p[j];
		}
		*(A+i)=exp(*(H+i));
	}
	*Am=0;
	for(i=0; i<L+1; i++){
		*Am+=*(A+i)**(p+i);
	}
}

void meiosi(int *g, int *gam){
	
	int i;
	double m=drand48();
	
	for(i=0; i<L; i++){
		gam[i]=(m<0.5)*g[i]+(m>=0.5)*g[i+L];
	}
}

void riproduzione_mutazione(int *gam1, int *gam2, int *gF){
	
	int i;
	double m=drand48();
	int w=(int)(drand48()*L);
	
	for(i=0; i<L; i++){
		gF[i]=gam1[i];
	}
	for(i=L; i<2*L; i++){
		gF[i]=gam2[i-L];
	}
	
	gF[w]=(m<mu)*(1-gF[w])+(m>mu)*gF[w];
}

void sopravvivenza(double A, double Am, int *k){
	
	if(drand48()<A/Am){
		(*k)++;
	}
}

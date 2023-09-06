//Dieckman e Doebeli per organimsi diploidi con solo crattere ecologico e assortatività basata su carattere marker (6 cromosomi)
//riproduzione sessuata ma ermafrodita, in ogni accoppiamento il primo scelto decide con chi accoppiarsi
//meiosi senza crossing over
//forma del panorama adattivo gaussiana
//tasso di mutazione puntuale uguale per tutti i caratteri
//distribuzioni iniziali: delta in f=0 per E, delta in f=L/2 per A, campana centrata in f=L/2 per M

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void genotipi_delta(int **g);
void genotipi_random(int **g);
void fenotipi(int **g, int *f);
void distribuzione(double *p, int *f);
void fitness(double *H0, double *H, double *A, double *Am, double *p, int **g, int *f, int t);
int mating(int v, int fAv, int *fE, double *c);
void meiosi(int *g, int *gam);
void riproduzione_mutazione(int *gam1, int *gam2, int *gF);
void sopravvivenza(double A, double Am, int *k, int *c);

int N, L, T;
double b, r, J, R, a, mu;

int main(){
	
	N=100;  //numero individui
	L=20;  //lunghezza genoma (singolo cromosoma lungo L/2)
	T=1000;  //generazioni
	b=50;
	r=10;
	J=0;
	a=2;
	R=1;
	mu=0.001;
	
	int i, j, k, v, w, control;
	int t;
	int **gE;  //carattere ecologico
	int *fE;
	int **gA;  //carattere assortativo
	int *fA;
	int **gM;  //carattere marker
	int *fM;
	int *gam1, *gam2;
	int **gFE;
	int **gFA;
	int **gFM;
	int fF;
	double *pE;
	double *pA;
	double *pM;
	double *H0;
	double *H;
	double *A;
	double Am;
	int *mark;  //marker per l'accoppiamento
	int markTOT;
	double *comp;  //"compatibilità"
	
	srand48(time(NULL));
	
	gE=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gE[i]=calloc(2*L, sizeof(int));
	}
	fE=calloc(N, sizeof(int));
	gam1=calloc(L, sizeof(int));
	gam2=calloc(L, sizeof(int));
	gFE=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gFE[i]=calloc(2*L, sizeof(int));
	}
	gA=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gA[i]=calloc(2*L, sizeof(int));
	}
	fA=calloc(N, sizeof(int));
	gFA=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gFA[i]=calloc(2*L, sizeof(int));
	}
	gM=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gM[i]=calloc(2*L, sizeof(int));
	}
	fM=calloc(N, sizeof(int));
	gFM=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		gFM[i]=calloc(2*L, sizeof(int));
	}
	pE=calloc(L+1, sizeof(double));
	pA=calloc(L+1, sizeof(double));
	pM=calloc(L+1, sizeof(double));
	H0=calloc(L+1, sizeof(double));
	H=calloc(L+1, sizeof(double));
	A=calloc(L+1, sizeof(double));
	mark=calloc(N, sizeof(int));
	comp=calloc(N, sizeof(double));
	
	for(i=0; i<L+1; i++){  //calcolo fitness
		H0[i]=b*exp(-pow(fabs((i-L/2)/r),a)/a);
	}
	
	genotipi_delta(gA);  //inizializza a delta in f=0 i genotipi assortativi
	genotipi_random(gM);
	fenotipi(gE, fE);
	fenotipi(gA, fA);
	fenotipi(gM, fM);
	
	for(t=0; t<T; t++){
		
		for(i=0; i<N; i++){
			mark[i]=0;
		}
		distribuzione(pE, fE);
		fitness(H0, H, A, &Am, pE, gE, fE, t);
		
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
				w=mating(v, fA[v], fM, comp);
				//w=(int)(drand48()*N);
			}
			mark[w]=1;
			
			meiosi(gE[v], gam1);  //riproduzione per carattere ecologico
			meiosi(gE[w], gam2);
			riproduzione_mutazione(gam1, gam2, gFE[k]);
			
			fF=0;
			for(i=0; i<L; i++){
				if((*(gFE[k]+i)==1)||(*(gFE[k]+i+L)==1)){
					fF++;
				}
			}
			sopravvivenza(A[fF], Am, &k, &control);
			
			if(control!=0){
				meiosi(gA[v], gam1);  //riproduzione per carattere assortativo
				meiosi(gA[w], gam2);
				riproduzione_mutazione(gam1, gam2, gFA[k-1]);
				meiosi(gM[v], gam1);  //riproduzione per carattere marker
				meiosi(gM[w], gam2);
				riproduzione_mutazione(gam1, gam2, gFM[k-1]);
			}
		}
		
		for(i=0; i<N; i++){  //sostituzione della generazione
			for(j=0; j<2*L; j++){
				*(*(gE+i)+j)=*(*(gFE+i)+j);
			}
		}
		for(i=0; i<N; i++){
			for(j=0; j<2*L; j++){
				*(*(gA+i)+j)=*(*(gFA+i)+j);
			}
		}
		for(i=0; i<N; i++){
			for(j=0; j<2*L; j++){
				*(*(gM+i)+j)=*(*(gFM+i)+j);
			}
		}
		fenotipi(gE, fE);
		fenotipi(gA, fA);
		fenotipi(gM, fM);
	}
	
	distribuzione(pE, fE);
	distribuzione(pA, fA);
	distribuzione(pM, fM);
	fitness(H0, H, A, &Am, pE, gE, fE, t);
	
	FILE *dati;
	dati=fopen("DD.dat", "w");
	fprintf(dati, "#fenotipo distribuzione_ecologico distribuzione_assortativo distribuzione_marker:\n");
	for(i=0; i<L+1; i++){
		fprintf(dati, "%d   %f   %f   %f\n", i, *(pE+i), *(pA+i), *(pM+i));
	}
	//fprintf(dati, "\n");
	//fprintf(dati, "#fenotipo distribuzione:\n");
	//for(i=0; i<L+1; i++){
		//fprintf(dati, "%d %f\n", i, *(pA+i));
	//}
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

void genotipi_delta(int **g){
	
	int i, j;
	int *s;
	int mark=0;
	double m, n;
	
	s=calloc(L/2, sizeof(int));
		
	for(i=0; i<L/2; i++){ 
		mark=0;
		while(mark==0){
			s[i]=(int)(drand48()*L/2);
			mark=1;
			for(j=0; j<i; j++){
				if(s[i]==s[j]){
					mark=0;
				}
			}
		}
	}
	
	for(i=0; i<N; i++){  //inizializza in modo casuale metà dei loci dominanti in entrambi o in un solo allele
		for(j=0; j<L/2; j++){
			m=drand48();
			if(m<0.5){
				n=drand48();
				if(n<0.5){
					*(g[i]+s[j])=1;
				}
				else{
					*(g[i]+s[j]+L)=1;			
				}
			}
			else{
				*(g[i]+s[j])=1;
				*(g[i]+s[j]+L)=1;			
			}
		}
	}
	
}

void genotipi_random(int **g){  //
	
	int i, j;
	int m;
	
	for(i=0; i<N; i++){
		for(j=0; j<L; j++){
			*(g[i]+j)=(int)(drand48()*2.0);
			if(*(g[i]+j)==1){
				m=drand48();
				if(m<0.5){
					*(g[i]+j+L)=(int)(drand48()*2.0);
				}
			}
		}	
	}
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

int mating(int v, int fAv, int *fE, double *c){
	
	int i;
	int w;  //individuo da selezionare
	double m;  //fenotipo calcolato come in Dieckman e Doebeli
	double sigma;  //deviazione standard della "curva delle preferenze"
	
	for(i=0; i<N; i++){
		c[i]=0;
	}
	
	if(fAv==L/2){  //calcolo della "compatibilità"
		for(i=0; i<N; i++){
			c[i]=i+1;
		}
		
	}else if(fAv>L/2){
		m=(double)(fAv-L/2)*2.0/(double)L;
		sigma=1.0/(20.0*m*m);
		//sigma=10;
		c[0]=exp(-pow((double)fE[i]-(double)fE[v], 2.0)/(2.0*sigma*sigma));
		for(i=1; i<N; i++){
			c[i]=c[i-1]+exp(-pow((double)fE[i]-(double)fE[v], 2.0)/(2.0*sigma*sigma));
		}
		
	}else{
		m=(double)(fAv-L/2)*2.0/(double)L;
		sigma=1.0/(m*m);
		//sigma=10;
		c[0]=1-exp(-pow((double)fE[i]-(double)fE[v], 2.0)/(2.0*sigma*sigma));
		for(i=1; i<N; i++){
			c[i]=c[i-1]+1-exp(-pow((double)fE[i]-(double)fE[v], 2.0)/(2.0*sigma*sigma));
		}
	}
	
	m=drand48()*c[N-1];  //scelta del partner
	if(m<c[0]){
		w=0;
	}else{
		for(i=1; i<N; i++){
			if((c[i-1]<=m)&&(m<c[i])){
				w=i;
			}
		}
	}

	return w;
}

void meiosi(int *g, int *gam){
	
	int i;
	double m=drand48();
	
	for(i=0; i<L; i++){
		gam[i]=(m<0.5)*g[i]+(m>=0.5)*g[i+L/2];
	}
}

void riproduzione_mutazione(int *gam1, int *gam2, int *gF){
	
	int i;
	double m=drand48();
	int w=(int)(drand48()*L);
	
	for(i=0; i<L; i++){
		gF[i]=gam1[i];
	}
	for(i=L/2; i<2*L; i++){
		gF[i]=gam2[i-L];
	}
	
	gF[w]=(m<mu)*(1-gF[w])+(m>mu)*gF[w];
}

void sopravvivenza(double A, double Am, int *k, int *c){
	
	*c=0;
	if(drand48()<A/Am){
		(*k)++;
		(*c)++;
	}
}

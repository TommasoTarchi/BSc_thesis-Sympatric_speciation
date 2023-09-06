//prima sopravvivenza e poi riproduzione
//sia riproduzione sessuata che asessuata
//specie aploidi
//no monogamia

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* header */

typedef struct _individual{
	
	char **g; // gene sets (0,1)
	int *f;  // phenotype
	
}individual;

// global environment 
struct{
	
	int N; // number of individuals
	int n; // number of traits
	int *lg;   // length of genome per trait
	int rec; // 1 = recombination (sex), 0: nosex
	int mating_trait; // preference trait
	int marker_trait;  // mating trait
	double mu;  // mutation rate
	double rho0; // initial prob genes=1

}env; 


individual * new_individual();
void recombine(individual *dest, individual *parent1, individual *parent2);
void mutate(individual *x);
void copy_individual(individual *dest, individual *src);
void compute_phenotype_distr(int *F, individual **P); 


int main(){
	  
	individual **P, **P1;
	int i, j, k;
	int t, T;
	double *H0, *H;
	double sumH, aveH;
	int *F; // phenotype distribution (ecological trait) 
	int score;// reproduction score 
	int maxscore, maxj;
 
	// definitions
	env.N=1000;
	env.n=3;
	env.lg=calloc(env.n, sizeof(int));
	env.lg[0]=100;
	env.lg[1]=10;
	env.lg[2]=10;  //=0 per random mating
	env.rec=1;  //=0 per riproduzione asessuata, =1 per riproduzione sessuata
	env.mating_trait=2;
	env.marker_trait=0;  //=0 per mating secondo carattere ecologico, =1 per mating secondo carattere non ecologico
	env.mu=0.01;
	env.rho0=0.7; // initial density of genes=1
	T=10000;
	// fitness
	double r=400;
	double b=10;
	double R=20;
	double a=2; 
	double J=0.01;

	P=calloc(env.N, sizeof(individual*));
	P1=calloc(env.N, sizeof(individual*));
	F=calloc(env.lg[0]+1, sizeof(int)); // additive effects of genes
	H0=calloc(env.lg[0]+1, sizeof(double)); 
	H=calloc(env.lg[0]+1, sizeof(double)); 
  
	// inizialization
	srand48(time(NULL));
  
	for(i=0; i<env.N; i++){
		P[i]=new_individual(); // random individuals
		P1[i]=new_individual();
	}
 
	// static fitness
	for(j=0; j<=env.lg[0]; j++){ // 0 is the ecological trait
		//H0[j]=b*(1-fabs(i-env.lg[0]/2.)/r-1/(1+fabs(i-env.lg[0]/2.)/r));
		H0[j]=b*exp(-pow((j-env.lg[0]/2.)/r,2)/2);
		
	}
    
	// compute phenotype distribution 
	compute_phenotype_distr(F, P);

	// main loop 
	for(t=0; t<T; t++){
		
		// compute fitness
		for(j=0; j<=env.lg[0]; j++){ 
			H[j]=H0[j];
			for(k=0; k<=env.lg[0]; k++){
				H[j]-=J*exp(-pow(fabs((j-k)/R),a)/a)*F[k];
			}
		}
		// compute average fitness
		sumH=0;
		for(j=0; j<=env.lg[0]; j++){
			sumH+=F[j]*H[j];
		}
		aveH=sumH/env.N;
   
		// survival 
		j=0;
		for(i=0; i<env.N; i++){
			if(drand48()<H[P[i]->f[0]]/aveH){
				copy_individual(P1[j], P[i]);
				j++;
			}
		}
    
		// complete the population
		while(j<env.N){
			i=(int)(drand48()*env.N); // random individual
			if(drand48()<H[P[i]->f[0]]/aveH){
				copy_individual(P1[j], P[i]);
				j++;
			}
		}
		
		for(i=0; i<env.N; i++){ // for each individual
			if(env.rec){ // preferential mating
				// mating 
				if(env.lg[2]>0){
				// this is quite slow: order N^2 
				// find max score of other individuals (but himself)
					maxscore=0;
					maxj=0;                                                        //NIENTE MONOGAMIA???
					for(j=0; j<env.N; j++){
						if(j==i)continue; // avoid self-mating
						score=(env.lg[env.marker_trait]-abs(P1[j]->f[env.marker_trait]-P1[i]->f[env.marker_trait]))*(P1[i]->f[env.mating_trait]-env.lg[env.mating_trait]/2.); // only assortative
						if(score>maxscore){
							maxscore=score;
							maxj=j;
						}
					}
					// choose mate
					j=maxj;
					//printf("i=%d (%d:%d), j=%d (%d), maxscore=%d\n",i, P1[i]->f[env.mating_trait], P1[i]->f[env.marker_trait],j,P1[j]->f[env.marker_trait],maxscore);
				}else{ // random mating
					j=i;
					while(j==i){
						j=(int)(drand48()*env.N);
					}
				}
				recombine(P[i], P1[i], P1[j]);
			}else{ // nosex, only mutation
				for(i=0; i<env.N; i++){ // for each individual
					copy_individual(P[i], P1[i]);
					mutate(P[i]);
				}
			}
		}
 
		// compute phenotype distribution 
		compute_phenotype_distr(F, P);
	}
	
	FILE *dati;
	dati=fopen("DD.dat", "w");
	fprintf(dati, "#fenotipo distribuzione_ecologico distribuzione_assortativo distribuzione_marker:\n");
	for(i=0; i<env.lg[0]+1; i++){
		fprintf(dati, "%d   %d   %f\n", i, F[i], H[i]);
	}
	fclose(dati);
}
   
void compute_individual_phenotype(individual *x){
	
	int k, j;
  
	for(k=0; k<env.n; k++){
		x->f[k]=0;
		for(j=0; j<env.lg[k]; j++){
			x->f[k]+=x->g[k][j];
		}
	}
}

void recombine(individual *dest, individual *parent1, individual *parent2){

	int k, j;
  
	for(k=0; k<env.n; k++){
		for(j=0; j<env.lg[k]; j++){
			//recombination
			dest->g[k][j]=drand48()<0.5 ? parent1->g[k][j] : parent2->g[k][j];
		}
	}
	mutate(dest);
}
 
void mutate(individual *x){
	
	int k, j;
  
	for(k=0; k<env.n; k++){
		for(j=0; j<env.lg[k]; j++){
			//mutations
			x->g[k][j]^=drand48()<env.mu/env.lg[k];
		}
	}
	compute_individual_phenotype(x);
}
        
individual * new_individual(){
	
	int k, j;
	individual *x;
  
	x=calloc(1, sizeof(individual));
	x->g=calloc(env.n, sizeof(char*));
	x->f=calloc(env.n, sizeof(int));
	for(k=0; k<env.n; k++){
		x->g[k]=calloc(env.lg[k], sizeof(char));
		for(j=0; j<env.lg[k]; j++){
			x->g[k][j]=drand48()<env.rho0;
		}
	}
	compute_individual_phenotype(x);

	return(x);
}

void copy_individual(individual *dest, individual *src){
	
	int k, j;

	for(k=0; k<env.n; k++){
		dest->f[k]=src->f[k];
		for(j=0; j<env.lg[k]; j++){
			dest->g[k][j] = src->g[k][j];
		}
	}
}

void compute_phenotype_distr(int *F, individual **P){
	
	int i, k;
  
	for(k=0; k<=env.lg[0]; k++){
		F[k]=0;
	}
  
	for(i=0; i<env.N; i++){
		F[P[i]->f[0]]++;
	}
}  

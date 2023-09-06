//riferimento teorico: "sympatric speciation trough assortative mating in a long range cellular
//automaton" di Bagnoli e Guardiani

//i fenotipi sono approssimati come somma dei valori booleani associati ai loci e il numero di
//individui rimane costante

//la riproduzione/selezione avviene in questo modo: per N volte scelgo una coppia random di vecchi
//e li faccio accoppiare producendo un cucciolo (tramite ricombinazione uniforme), che poi faccio
//mutare su un gene a caso con probabilità mu; calcolo la fitness e "testo" in ordine casuale i
//cuccioli, eseguendo un altro accoppiamento random tra vecchi ogni volta che un cucciolo muore
//(cioè non supera il test della fitness) per rimpiazzarlo

//IL CODICE FUNZIONA SE SI ESCLUDE LA COMPETIZIONE (SE AD ESEMPIO SI IMPONE UNA CONDIZIONE SEMPRE VERA IN TEST,
//SI OSSERVA LA FORMAZIONE, DOPO MOLTE GENRAZIONI, DI UNA DELTA IN UN VALORE FENOTIPICO PIU' O MENO CENTRALE,
//COME ATTESO DALLA TEORIA): IL PROBLEMA SEMBRA STARE NEI VALORI TROPPO ESTREMI ASSUNTI DALLA FITNESS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define N 100 //scelta parametri
#define L 20
#define T 3e+4
#define DELTA 1
#define beta 1
#define gamma 10
#define J 1
#define alfa 10
#define R 2
#define mu0 1e-1
#define muINF 1e-06
#define delta 3000
#define tau 1e+04

int main(){
	
	int t=0;
	int **vecchi_g;  //genomi
	int **giovani_g;
	int **appoggio_g;
	int *vecchi_f;  //fenotipi
	int *giovani_f;
	int *appoggio_f;
	int *n;  //istogramma
	double *A;  //fitness
	double H;
	double Am;
	double mu; //mutation rate
	int i;
	int j;
	int index1;
	int index2;
	double k;
	double k1;
	int *s;  //indici casuali
	int mark;
	int k2;
	
	s=calloc(N, sizeof(int));
	
	vecchi_g=calloc(N, sizeof(int*));  //allocazione dello spazio genotipico
	for(i=0; i<N; i++){
		*(vecchi_g+i)=calloc(L, sizeof(int));
	}
	giovani_g=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		*(giovani_g+i)=calloc(L, sizeof(int));
	}
	appoggio_g=calloc(N, sizeof(int*));
	for(i=0; i<N; i++){
		*(appoggio_g+i)=calloc(L, sizeof(int));
	}
	
	vecchi_f=calloc(N, sizeof(int));  //allocazione dello spazio fenotipico
	giovani_f=calloc(N, sizeof(int));
	appoggio_f=calloc(N, sizeof(int));
	n=calloc(L+1, sizeof(int));
	
	srand(time(NULL));  //costruzione casuale dei genomi vecchi
	for(i=0; i<N; i++){
		for(j=0; j<L; j++){
			*(*(vecchi_g+i)+j)=(drand48()>0.5);
		}
	}
	
	for(i=0; i<N; i++){  //calcolo dei fenotipi vecchi
		for(j=0; j<L; j++){
			*(vecchi_f+i)=0;
		}
	}
	for(i=0; i<N; i++){
		for(j=0; j<L; j++){
			*(vecchi_f+i)+=*(*(vecchi_g+i)+j);
		}
	}
	
	A=calloc(N, sizeof(double));  //allocazione spazio fitness
	
	while(t<(int)T){
		
		mu=(((double)mu0-(double)muINF)/2)*((double)1-tanh(((double)t-(double)tau)/(double)delta))+(double)muINF;  //probabilità di mutazione di un gene casuale
		
		for(i=0; i<N; i++){  //prima "cucciolata"
			
			index1=(int)(drand48()*(double)N);  //scelta dei genitori
			index2=-1;
			while((index2<0)||(index2==index1)||(index2>N-1)){
				index2=index1-(int)DELTA+(int)(drand48()*(((double)2*(double)DELTA)+(double)1));
			}
			
			for(j=0; j<L; j++){  //ricombinazione uniforme
				k=drand48();
				*(*(giovani_g+i)+j)=(k<0.5)*(*(*(vecchi_g+index1)+j))+(k>=0.5)*(*(*(vecchi_g+index2)+j));
			}
			
			k1=drand48();  //mutazione
			if(drand48()<mu){
				*(*(giovani_g+i)+(int)(k1*L))=-(*(*(giovani_g+i)+(int)(k1*L)))+1;
			}
		}
		
		for(i=0; i<N; i++){  //calcolo dei fenotipi giovani
			for(j=0; j<L; j++){
				*(giovani_f+i)=0;
			}
		}
		for(i=0; i<N; i++){
			for(j=0; j<L; j++){
				*(giovani_f+i)+=*(*(giovani_g+i)+j);
			}
		}
		
		for(i=0; i<L+1; i++){  //costruisco fitness (provvisoria) e istogramma della nuova generazione
			*(n+i)=0;
		}
		for(i=0; i<L+1; i++){
			for(j=0; j<N; j++){
				if(*(giovani_f+j)==i){
					(*(n+i))++;
				}
			}
		}
		for(i=0; i<N; i++){
			*(A+i)=exp(exp(-((double)1/(double)beta)*pow(((*(giovani_f+i))/(double)gamma), (double)beta)));
			if(i==0){
				for(j=1; j<N; j++){
					H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
					*(A+i)*=exp(H);
				}
			}else if(i==N-1){
				for(j=0; j<N-1; j++){
					H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
					*(A+i)*=exp(H);
				}
			}else{
				for(j=0; j<i; j++){
					H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
					*(A+i)*=exp(H);
				}
				for(j=i+1; j<N; j++){
					H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
					*(A+i)*=exp(H);
				}
			}
		}
		Am=0;  //(calcolo fitness media)
		for(i=0; i<N; i++){
			Am+=*(A+i);
		}
		Am=Am/(double)N;
		
		k2=0;

		for(i=0; i<N; i++){  //ordino in modo casuale i numeri da 0 a N-1
			mark=0;
			while(mark==0){
				*(s+i)=(int)(drand48()*N);
				mark=1;
				for(j=0; j<i; j++){
					if(*(s+i)==*(s+j)){
						mark=0;
					}
				}
			}
		}
		
		while(k2<N){  //NON FUNZIONA PERCHE' SI OTTENGONO NUMERI TROPPO PICCOLI PER LA FITNESS: IL COMPUTER NON RIESCE A PROCESSARLI OPPURE LA PROBABILITA' CHE I NUOVI NATI PASSINO IL TEST TENDE A ZERO
				
			if(2<3){  //(test)
				*(appoggio_f+k2)=*(giovani_f+(*(s+k2)));
				for(j=0; j<L; j++){
					*(*(appoggio_g+k2)+j)=*(*(giovani_g+(*(s+k2)))+j);
				}
				k2++;
				              //VA CORRETTO IL FATTO CHE QUANDO UN INDIVIDUO VIENE ELIMINATO, QUELLO CHE LO SOSTITUISCE VIENE SEMPRE "SCELTO" PER IL TEST SUCCESSIVO
			}else{  //sostituzione dell'individuo eliminato
				index1=(int)(drand48()*(double)N);  //(scelta genitori)
				index2=-1;
				while((index2<0)||(index2==index1)||(index2>N-1)){
					index2=index1-(int)DELTA+(int)(drand48()*(double)2*((double)R+(double)1));
				}
				for(j=0; j<L; j++){  //(ricombinazione uniforme)
					k=drand48();
					*(*(giovani_g+(*(s+k2)))+j)=(k<0.5)*(*(*(vecchi_g+index1)+j))+(k>=0.5)*(*(*(vecchi_g+index2)+j));
				}
				k1=drand48();  //(mutazione)
				if(drand48()<mu){
					*(*(giovani_g+(*(s+k2)))+(int)(k1*L))=-(*(*(giovani_g+(*(s+k2)))+(int)(k1*L)))+1;
				}
				for(j=0; j<L; j++){
					*(giovani_f+(*(s+k2)))+=*(*(giovani_g+(*(s+k2)))+j);
				}
				for(i=0; i<L+1; i++){  //(ricalcolo fitness)
					*(n+i)=0;
				}
				for(i=0; i<L+1; i++){
					for(j=0; j<N; j++){
						if(*(giovani_f+j)==i){
							(*(n+i))++;
						}
					}
				}
				for(i=0; i<N; i++){
					*(A+i)=exp(exp(-((double)1/(double)beta)*pow(((*(giovani_f+i))/(double)gamma), (double)beta)));
					if(i==0){
						for(j=1; j<N; j++){
							H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
							*(A+i)*=exp(H);
						}
					}else if(i==N-1){
						for(j=0; j<N-1; j++){
							H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
							*(A+i)*=exp(H);
						}
					}else{
						for(j=0; j<i; j++){
							H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
							*(A+i)*=exp(H);
						}
						for(j=i+1; j<N; j++){
							H=-(double)J*((double)(*(n+j))/(double)N)*exp(-((double)1/(double)alfa)*pow(fabs((*(giovani_f+j)-*(giovani_f+i))/(double)R), (double)alfa));
							*(A+i)*=exp(H);
						}
					}
				}
				Am=0;  //(calcolo fitness media)
				for(i=0; i<N; i++){
					Am+=*(A+i);
				}
				Am=Am/(double)N;
			}
		}
		
		for(i=0; i<N; i++){  //sostituzione della vecchia generazione con la nuova
			*(vecchi_f+i)=*(appoggio_f+i);
			for(j=0; j<L; j++){
				*(*(vecchi_g+i)+j)=*(*(appoggio_g+i)+j);
			}
		}
		
		t++;
	}
	
	for(i=0; i<N; i++){
		*(n+i)=0;
	}
	for(i=0; i<L+1; i++){  //costruzione istogramma
		for(j=0; j<N; j++){
			if(*(vecchi_f+j)==i){
					(*(n+i))++;
			}
		}
	}
	
	FILE *istogramma;
	istogramma=fopen("BG.dat", "w");
	fprintf(istogramma, "#posizione_array_vecchi fenotipo:\n");
	for(i=0; i<N; i++){
		fprintf(istogramma, "%d %d\n", i, *(vecchi_f+i));
	}
	fprintf(istogramma, "#posizione_array_appoggio fenotipo:\n");
	for(i=0; i<N; i++){
		fprintf(istogramma, "%d %d\n", i, *(appoggio_f+i));
	}
	fprintf(istogramma, "#posizione_array_giovani fenotipo:\n");
	for(i=0; i<N; i++){
		fprintf(istogramma, "%d %d\n", i, *(giovani_f+i));
	}
	fprintf(istogramma, "#array_numeri_casuali:\n");
	for(i=0; i<N; i++){
		fprintf(istogramma, "%d\n", *(s+i));
	}
	fprintf(istogramma, "#individuo fitness:\n");
	for(i=0; i<N; i++){
		fprintf(istogramma, "%d  %f\n", i, *(A+i));
	}
	fprintf(istogramma, "fitness_media:  %f\n", Am);

//	FILE *istogramma;
//	istogramma=fopen("BG.dat", "w");
	fprintf(istogramma, "#fenotipo popolazione:\n");
	for(i=0; i<L+1; i++){
		fprintf(istogramma, "%d %d\n", i, *(n+i));
	}
	fclose(istogramma);
	
	return 0;
}	

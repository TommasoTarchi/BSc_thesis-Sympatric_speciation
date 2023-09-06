//modello 1 con modifica: A,-A rinunciano ad accoppiarsi piuttosto che farlo con B,-B

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
	
	int t=0;
	double nA=(double)1/6;
	double Oa=(double)1/6;  //femmine senza preferenze
	double A=(double)1/6;
	double nB=(double)1/6;
	double Ob=(double)1/6;  //femmine senza preferenze
	double B=(double)1/6;
	
	double p_Oa=0;
	double fraz=0;
	
	int t_finale=10;  //numero di generazioni scelte
		
	while(t<t_finale+1){
		
		p_Oa=Oa/3;  //(Oa non hanno preferenze)
		
		if(A+p_Oa<=B){  //B abbondanti
			
			B=A+p_Oa;
			
		}else{  //B scarsi
			
			fraz=A/p_Oa;  //criterio di preferenza proporzionale
			
			if(fraz>=1){
				
				Oa=p_Oa+2*(B/(fraz+1));
				A=(B/(fraz+1))*fraz;
				
			}else{
				
				fraz=p_Oa/A;
				A=B/(fraz+1);
				Oa=p_Oa+2*(B/(fraz+1))*fraz;
				
			}
		}
		
		if(4*p_Oa-Oa<=Ob){  //Ob abbondanti
			
			Ob=4*p_Oa-Oa;
			Oa=3*p_Oa;
			
		}else{  //Ob scarsi
			
			Oa=Oa-p_Oa+Ob;
			
		}
		
		nB=B;  //(la distribuzione è totalmente simmetrica)
		nA=A;
		
		t=t+1;
		
	}
	
	printf("-A: %lf, Oa: %lf, A: %lf\n\n-B: %lf, Ob: %lf, B: %lf\n", nA, Oa, A, nB, Ob, B);
	
	return 0;
}

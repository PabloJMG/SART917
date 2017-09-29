#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <fenv.h>

#include "../headers/init_coe.h"
#include "../headers/poinit_coe.h"
#warning "Compilation Error"

//Hacer la proyección primero y aplicarla a cada píxel y no hacer la proyección para cada píxel.

#define PATH "/home/pablo/local.git/pablo"
#define pi 3.141592
#define DOSPI
//~ #define COEF_COM
#define GUARDARC
#define PARALLEL_BEAM
#define DIMX 100
#define DIMY 100
#define FDIM 1.2
#define NANG 50
#define NPIX 20
#define LIMX (int)(DIMX*FDIM)
#define LIMY (int)(DIMY*FDIM)
#define LIMTOT (int) (LIMX*LIMY)
//~ #define LIMTOT 20

#define D2N (int)(DIMX*DIMY*NANG)

//~ #define D2N 40
#define ind_coe(a, b, c) ((a*DIMY+b)*NANG +c)  //Aqui da problemas, tienes que pensar bien cómo hacerlo
#define ind_coe2(a,b) (a*LIMY+b) //Para n_gamma
#define ind_coe3(a,b) (a*LIMTOT+b) //Para coe y para gamma al parecer
#define coe_gam(a,b) (a*LIMTOT+b)
void *carga_archivo(int NumAng, int NumPix, float *b);
FILE *apertura_archivo(int contador);
FILE *apertura_archivo2(int cont1, int cont2);
void saltar_com(FILE *fin);

//~ typedef struct nopze{
	//~ int nonnum;
	//~ size_t size;
	//~ int ang[];
	
//~ }nopze;

//~ struct nopze2{
	//~ int nonnum;
	//~ int ang[20];
	
//~ };

//~ struct sparse_mat{
	
	//~ float *nozero;
	//~ int *space;
	
//~ }; 
//~ inline
//~ nopze* alloc_nopze(int a, size_t n) {
  //~ nopze * ret = calloc(sizeof(nopze) + n, 1);
  //~ if (ret) memcpy(ret,
                  //~ &(nopze const){ .nonnum=a, .size= n},
                  //~ sizeof(nopze));
  //~ return ret;
//~ } // a = 0; n = NumPix

int main(int argc, char** argv)
{	
	
	int NumPix, gammcon, pro, cont, cont2, cont3, cont4, cont5, cont6, dimx, dimy, limx, limy, pixx, pixy, NumAng, facx, facxt, facy, facyt, i, valu, NumRayos1, ng, dimtotal, nozcoe, suma, ac1, ac2, c1, c2, c3, gc1, gc2, bc1, bc2, d2;
	int aint[3];
	float beta, cosbeta, sinbeta, difgam, factora, deltafantes, deltaang, lambda, deltax, deltay, deltaf, posx, posy, posx2, posy2, diferencia, cant, sumpix, factorpos, lol;
	char filename[100];
	FILE *finp, *pipeplot;
	clock_t inicio, final, com, fin;
	time_t archtime;
	float factordim = FDIM;
	
	printf("Limx = %d, %d", LIMX, ind_coe(2,3,1));
	//~ getchar();
	//~ feenableexcept(FE_DIVBYZERO| FE_INVALID| FE_OVERFLOW);
	//~ exit(1);
	#ifdef DIM_TXT
	finp=fopen("../FP/input/dimSART.txt", "r");
	fscanf(finp, "%d". &dimx);
	fscanf(finp, "%d", &dimy);
	fclose(finp);
	#endif
	
	dimx=DIMX;
	dimy=DIMY;
	d2=D2N;
	lambda=0.01;
	finp=fopen("../FP/input/inputint.txt", "r");
	i=0; // ¿?
	for(cont=0; cont<3; cont++)
		{
			
			saltar_com(finp);
			fscanf(finp, "%d", &aint[cont]);
			printf("%d\n", aint[cont]);
			saltar_com(finp);
		}
	fclose(finp);
	printf("%d x %d\n", dimx, dimy);
	#ifdef DOSPI
	NumAng=aint[0];
	#endif
	
	#ifdef ABANICO
	double angulos[10];
	finp=fopen("NumAngAbanico.txt", "r");
	fscanf(finp, "%d", &NumAng);
	fclose(finp);
	#endif

	NumPix=aint[2];
	
	printf("Antes de la declaracion de structs\n");
	//~ NumAng = aint;
	factora=1000;
	//~ NumRayos=NANG*NumPix;
	NumRayos1=(NANG+1)*NumPix;
	//~ limx=(int)rint(factordim*dimx);
	//~ limy=(int)rint(factordim*dimy);
	limx=LIMX;
	limy=LIMY;
	dimtotal=LIMTOT;
    facx=(int)rint((factordim-1)*dimx/2);
	facy=(int)rint((factordim-1)*dimy/2);
	facxt=(int)rint((factordim+1)*dimx/2);
	facyt=(int)rint((factordim+1)*dimx/2);
	factorpos=(float)dimy/NumPix; 
	deltaang=2*pi/NumAng;
	printf("limx %d limy %d facs %d %d, LIMTOT %d y D2N %d, NumAng = %d\n", limx, limy, facx, facxt, LIMTOT, D2N, NumAng);
	size_t tamdimtotal = dimtotal*sizeof(float);
	  //~ exit(1);
	#ifdef DETORNOT
	if(NumRayos<(dimx*dimy))
		{
			printf("Undetermined system\n");
			return -1;
		}
    #endif
			//~ exit(1);
	//~ nopze n_gamm[LIMTOT];
	int no_gam[LIMTOT];
	//~ int no_gama[LIMTOT*NumPix];
	
	int *no_gama=(int*) calloc(LIMTOT*NumPix, sizeof(int));
	for(cont=0; cont<LIMTOT; cont++)
	{
		no_gam[cont]=0;
		
	}
	//~ for(cont=0; cont<LIMTOT; cont++)
	//~ {
		
		//~ n_gamm[cont]=alloc_nopze(0, NumPix);
		
	//~ }

	//~ nopze *a=alloc_nopze(0, NumPix);

	//~ nopze nzcoe[D2N]; //ind_coe definido para este siempre y cuando a vaya de 0 a dimx*dimy y b vaya de 0 a NAng
	int nzcoe[D2N];
	int *nzcoea=(int *)calloc(NumPix*d2, sizeof(int)); //runtime variable ¿?¿?¿?¿?¿?¿?
	float *gamma=(float *)calloc(NumAng*dimtotal, sizeof(float)); //valgrind echoes error
	float *b=(float *)calloc(NumAng*NumPix, sizeof(float)); //Si const b, quitar esto
	float *betar=(float *)calloc(NumAng*NumPix, sizeof(float)); 
	float *coe=(float *)calloc(dimtotal*NumRayos1,sizeof(float)); //valgrind error
	float *x = (float *) calloc(dimtotal, sizeof(float));
		
	//~ #pragma omp parallel
	//~ {
			//~ #pragma omp for private(cont, cont2)
			//~ for(cont=0; cont<limx;cont++)
			//~ {
				
				//~ for(cont2=0; cont2<limy; cont2++)
				//~ {
					
					
						//~ n_gamm[ind_coe2(cont, cont2)].nonnum=0;
						//~ n_gamm[ind_coe2(cont, cont2)].ang=calloc(NANG, sizeof(int)); //valgrind error
					
				//~ }
			
			
			//~ }
		
		
		
	//~ }
	//~ #pragma omp parallel
	//~ {
		//~ #pragma omp for private(cont, cont2, cont3)
		//~ for(cont=0; cont<dimx; cont++)
		//~ {
			
			//~ for(cont2=0; cont2<dimy; cont2++)
			//~ {
				
				//~ for(cont3=0; cont3<NumAng; cont3++)
				//~ {
					
					//~ nzcoe[ind_coe(cont, cont2, cont3)].ang=calloc(NumPix, sizeof(int));
										
				//~ }				
				
			//~ }
								
		//~ }				
		
	//~ }
	
	#pragma omp parallel
	{
		
		#pragma omp for private(cont, cont2, finp, valu)
			for(cont=0; cont<NumAng; cont++)
			{
				//Nada claro esto de aqui abajo
				//~ gamma[cont]=calloc(dimtotal, sizeof(gamma[0]));
				//~ betar[cont]=calloc(NumPix, sizeof(betar[cont]));
				//~ b[cont]=calloc(NumPix,sizeof(b[0]));
				finp=apertura_archivo(cont);
				
				
				for(cont2=0; cont2<NumPix; cont2++)
				{
					
					if(cont2==0)
					b[cont*NumPix+cont2]=log(1+valu)+log(2);
					else
					{
						
						fscanf(finp, "%d", &valu);
						b[cont*NumPix+cont2]=log(1+valu);
						
					}
					
					
				}
				
				fclose(finp);
			}
				
		
	}
	
	for(cont=0; cont<dimtotal; cont++)
	x[cont]=0.f;
	
	printf("coe = %f", coe[ind_coe3(29,254446)]);
	#ifdef PARALLEL_BEAM
	
	#ifdef DOSPI
	
	#pragma omp parallel 
	{
		
		#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				//Se entra aqui de momento
				pa_initialize_coef_flo_2pi(cont, NumPix, NumAng, factordim, limy, dimx, dimy, LIMTOT, NumRayos1, deltaang, factorpos,  &coe[0]);

			}
		
	}
	
	#endif
	
	#ifdef ABANICO
	
	#pragma omp parallel 
	{
		
			#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				pa_initialize_coef_flo_abanico(cont, NumPix, NumAng, factordim, limy, dimx, dimy, LIMTOT, NumRayos1, deltaang, factorpos, angulos, coe);
			}
		
	}
	
	
	#endif
	
	
	
	#endif
	


	
	
	#ifdef POINT_SOURCE
	
	#ifdef DOSPI
	
	#pragma omp parallel 
	{
		
		#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				
				po_initialize_coef_flo_2pi(cont, NumPix, NumAng, factordim, limy, dimx, dimy, LIMTOT, NumRayos1, deltaang, factorpos,  coe);
			}
		
	}
	
	#endif
	
	#ifdef ABANICO
	
	#pragma omp parallel 
	{
		
		 	#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				po_initialize_coef_flo_abanico(cont, NumPix, NumAng, factordim, limy, dimx, dimy, LIMTOT, NumRayos1, deltaang, factorpos, angulos, coe);
			}
		
	}
	
	
	#endif
	
	
	
	#endif
	
	finp=fopen("../res/Suma_Coef_SART.txt", "w+");
	for(cont=0;	cont<limx; cont++)
	{
		
		for(cont2=0; cont2<limy; cont2++){
			c1=NumAng*NumPix;
			c2=cont*limy+cont2;	
			fprintf(finp, "%d %d %f\n", cont, cont2, coe[ind_coe3(c1, c2)]);
		}
	}
	fclose(finp);
	
	//~ #pragma omp parallel
		//~ {
			
			//~ #pragma omp for private(cont, finp, cont2, cont3, cont4)
				//~ for(cont=0; cont<NumAng; cont++)
					//~ {
						//~ for(cont2=0; cont2<NumPix; cont2++)
							//~ {
								
								//~ finp=apertura_archivo2(cont, cont2);
									//~ for(cont3=0; cont3<limx; cont3++)
										//~ {
											
											//~ for(cont4=0; cont4<limy; cont4++)
												//~ fprintf(finp, "%d %d %f\n", cont3, cont4, coe[cont*NumPix+cont2][cont3*limy+cont4]);
											
											
										//~ }
								
								
							//~ }
						
						
						
					//~ }

			
			
		//~ }
		
	for(cont=0; cont<NumAng; cont++)
	{
		for(cont2=0; cont2<NumPix; cont2++)
		{
			
			for(cont3=facx; cont3<facxt; cont3++)
			{
				for(cont4=facy; cont4<facyt; cont4++)
				{
					c1=cont*NumPix+cont2;
					c2=cont3*limy+cont4;
					betar[cont*NumPix+cont2]+= (coe[ind_coe3(c1,c2)]);
				}
				
				
			}
			
		}
		
		
	}
		
	gammcon=0;	
		
	for(cont=0; cont<NumAng; cont++)
	{
		
		for(cont2=facx; cont2<facxt; cont2++)
		{
			for(cont3=facy; cont3<facyt; cont3++)
			{
				for(cont4=0; cont4<NumPix; cont4++)
				{
			
					c1=cont*NumPix+cont4;
					c2=cont2*limy+cont3;
					gc1=cont;
					gc2=cont2*limy+cont3;
					gamma[ind_coe3(gc1,gc2)]+=(coe[ind_coe3(c1,c2)]);
					
				}
				if(gamma[ind_coe3(gc1,gc2)]<0.0001)
					gammcon++;
			}
			
			
			
		}
		
	}
		
	i=0;
    deltaf=1;
	//~ finp=fopen("../res/debug_data.txt", "at");
	//~ fprintf(finp, "%f %f %f %f %f gammcon %d\n", coe[ind_coe3(3000,445)], b[20*NumPix+100], betar[2*NumPix+2], gamma[ind_coe3(0,665)], gamma[ind_coe3(0,666)], gammcon);
	//~ fclose(finp);
	//~ for(cont=0; cont<limx; cont++)
	//~ {
		//~ for(cont2=0; cont2<limy; cont2++)
			//~ x[cont*limy+cont2]=0.f;
	//~ }
	printf("Dim %d, resta %d\n", dimx, (facxt-facx));
	for(cont2=facx; cont2<facxt; cont2++)
	{
		
		for(cont3=facy; cont3<facyt; cont3++)
		{
			
			for(cont=0; cont<NumAng; cont++)
			{
				nozcoe=0;  
					for(cont4=0; cont4<NumPix; cont4++)
					{	c1=cont*NumPix+cont4;
						c2=cont2*limy+cont3;
						if(coe[ind_coe3(c1,c2)]>0.f)
						{
							//~ nzcoe[(cont2-facx)*dimy+(cont3-facy)][cont].ang[nozcoe]=cont4;
							nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+nozcoe]=cont4;
							nozcoe++;
							
						}
						
					}
				
				nzcoe[ind_coe((cont2-facx),(cont3-facy),cont)]=nozcoe;
				
			}
		
			
			
		}
		
		
	}
	suma=0;
	for(cont2=facx; cont2<facxt; cont2++)
	{
		for(cont3=facy; cont3<facyt; cont3++)
		{
				for(cont4=0; cont4<NumPix; cont4++)
					suma+=nzcoe[ind_coe((cont2-facx),(cont3-facy),cont)];
		}


	}
        printf("Sparcity = %d de %d\n", suma, NumAng*NumPix*dimx*dimy);
        //~ getchar();
//Recorre las gamma (suma de las columnas de la matriz de coef) buscando las no nulas
	for(cont2=facx; cont2<facxt; cont2++)
	{
	
		for(cont3=facy; cont3<facyt; cont3++)
		{

				ng=0;
				for(cont=0; cont<NumAng; cont++)
				{		gc2=cont2*limy+cont3;
						if(gamma[ind_coe3(cont, gc2)] > 0.f)
						{
							no_gama[ind_coe2(cont2, cont3)*NumPix+ng]=cont;
							ng++;
								//~ printf("Se cumple lo de gamma");
						}
				}
				no_gam[ind_coe2(cont2, cont3)]=ng;

		}
	
	}
	betar[0]=0.01;
	inicio=clock();
	deltafantes=100; deltaf=10;
	printf("Comienza iter\n");
	pipeplot=popen("gnuplot -persist","w");
	fprintf(pipeplot, "set palette gray negative\n");
	
	while((deltaf>0.1) && (i<50) && (deltafantes>deltaf))
	{
		deltafantes=deltaf;
		deltaf=0;
		sumpix=0;
		
		for(cont2=facx; cont2<facyt; cont2++)
		{
			
			for(cont3=facy; cont3<facyt; cont3++)
			{
				
				difgam=0;
					//~ com=clock();
				for(cont=0; cont<no_gam[ind_coe2(cont2, cont3)]; cont++)
				{
					diferencia=0;
					
					#pragma omp parallel for private(cont4, cant, cont5, cont6, bc1, bc2, c1, c2, c3, lol, ac1, ac2) reduction(+:diferencia) shared(coe)
					
					for(cont4=0; cont4<nzcoe[ind_coe((cont2-facx),(cont3-facy),cont)]; cont4++)
					{
						
						cant=0;
						c1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont]*NumPix+nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
						c2=cont2*limy+cont3;
						lol=coe[ind_coe3(c1,c2)];
						
						if(lol>0.f){
							
					
						
						//~ printf("Dentro bucle y cont = %d\n", cont);
						
						
						
							//~ printf("Pro es %d\n", pro);
							//~ getchar();
					
						for(cont5=facx; cont5<facxt; cont5++)
						{
							
							for(cont6=facy; cont6<facyt; cont6++)
							{
								ac1=cont*NumPix +nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
								ac2=cont5*limy+cont6;
								
								if(coe[ind_coe3(ac1,ac2)]>0)
								cant+=coe[ind_coe3(ac1,ac2)]*x[ac2]; //Aquí sparse_matrix_multip ¿?
								
								
							}
							
							
						}
		
						
						
						bc1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont];
						bc2=nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
						bc1=bc1*NumPix+bc2;
						c3=no_gama[ind_coe2(cont2, cont3)*NumPix+cont]*NumPix+nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
						diferencia+=(b[c3]-cant)*(coe[ind_coe3(c1,c2)])/(betar[bc1]); 
						
					}
					
					}
					
					
					gc1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont];
					gc2=cont2*limy+cont3;
					difgam+=diferencia/(gamma[ind_coe3(gc1, gc2)]);
				}
				//~ fin=clock();
				//~ printf("Update de pixel tarda %lf\n",(double)(fin-com)/(CLOCKS_PER_SEC));
				//~ getchar();
				x[cont2*limy+cont3]+=lambda*difgam;
				deltaf+=lambda*difgam;
				sumpix+=x[cont2*limy+cont3];
				 
			}  //Final bucle update de un pixel
			
			
			
		}
		
		deltaf=deltaf/sumpix;
		deltaf=fabs(deltaf);
		final=clock();
		printf("Deltaf=%lf y n de iteraciones =%d. Desde el inicio llevamos %lf segundos\n", deltaf, i,(double)(final-inicio)/(CLOCKS_PER_SEC*12));
        sprintf(filename, "../res/Resultados_SART/Resultados_SART_%d.txt", i);

		finp=fopen(filename, "w+");
		
		for(cont=0; cont<limx; cont++)
		{
			for(cont2=limy; cont2>0; cont2--)
				fprintf(finp, "%d %d %f\n", cont, limy-cont2, x[cont*limy+cont2]);
			
		}
		
		fclose(finp);

        fprintf(pipeplot, "plot '%s/res/Resultados_SART/Resultados_SART_%d.txt'w image \n",PATH, i);
        fprintf(pipeplot, "set term postscript \n");

        fprintf(pipeplot, "set output 'img/imSART/NumAng_%d_Numpix_%d,dim_%dx%d__%d.png' \n", NumAng, NumPix, dimx, dimy, archtime);
        fprintf(pipeplot, "replot\n");
        fflush(pipeplot);
        i+=1;

	}

	pclose(pipeplot);
	
	#ifdef GUARDARC
	
	sprintf(filename, "../res/SARTC/NAng_%d.txt", NumAng);
	finp=fopen(filename, "w+");
	
	 for(cont=0; cont<limy; cont++)
	{

		for(cont2=0; cont2<limx; cont2++)
		{
			fprintf(finp, "%d %d %f\n",cont, cont2, x[cont*limx+cont2]);
		}

	}

	fclose(finp);
	
	#endif
	
	//~ #pragma omp parallel
	//~ {
		
		//~ #pragma omp for private(cont)
		//~ for(cont=0; cont<NumRayos1; cont++)
		//~ {
			
			//~ free(coe[cont]); //valgrind error
			
		//~ }
				
	//~ }
	
		//~ #pragma omp parallel
	//~ {
		//~ #pragma omp for private(cont, cont2, cont3)
		//~ for(cont=0; cont<dimx; cont++)
		//~ {
			
			//~ for(cont2=0; cont2<dimy; cont2++)
			//~ {
				
				//~ for(cont3=0; cont3<NumAng; cont3++)
				//~ {
					
					//~ free(nzcoe[ind_coe(cont, cont2, cont3)].ang);
										
				//~ }				
				
			//~ }
								
		//~ }				
		
	//~ }
	
	
	//~ #pragma omp parallel
	//~ {
			//~ #pragma omp for private(cont, cont2)
			//~ for(cont=0; cont<limx;cont++)
			//~ {
				
				//~ for(cont2=0; cont2<limy; cont2++)
				//~ {
					
					
					
						//~ free(n_gamm[ind_coe2(cont, cont2)].ang);
					
				//~ }
			
			
			//~ }
		
		
		
	//~ }
	//~ free(no_gam);
	free(no_gama);
	free(gamma);
	free(coe);
	free(b);
	
	//~ for(cont=0; cont<NANG; cont++)
	//~ free(gamma[cont]);
	
	return 0;
}



void *carga_archivo(int NumAng, int NumPix, float *b) /*Esta función está en el main tal cual*/
{
        int cont, cont2, valu;

        FILE *finp;

 for(cont=0; cont<NumAng; cont++)
 {

         finp=apertura_archivo(cont);

         for(cont2=0; cont2<NumPix; cont2++)
         {
        fscanf(finp, "%d", &valu);

        b[cont*NumPix+cont2]=log(1+valu); /*Esto aquí lo he cambiado*/
        //~ b[cont][cont2]=valu;
        if(cont2==0)
        b[cont*NumPix+cont2]=log(2)+log(1+valu);

         }

         fclose(finp);
 }

        //~ printf("Archivo cargado, presione alguna tecla\n");
        //~ getchar();
}


FILE *apertura_archivo(int contador)
{

char nombrearchivo[250]="Resultados_pixeles.txt";

sprintf(nombrearchivo, "%s/res/angulo_%d.txt",PATH, contador);

if(fopen(nombrearchivo, "r")==NULL)
{
        perror("El archivo no existe o no se puede abrir\n");
}

return fopen(nombrearchivo, "r");



}

                                                              
void saltar_com(FILE *fin) /* En archivo de texto salta a la siguiente línea tras \n*/
{
        char col;
        while(fscanf(fin, "%c", &col))
        {
                if(col=='\n')
                {

                        break;
                }

        }


}

FILE *apertura_archivo2(int cont1, int cont2)
{

        char file_name[200];
        sprintf(file_name, "%s/res/Resultados_SART/Matriz_coef_SART_%d_%d.txt", PATH, cont1, cont2);


        return fopen(file_name, "w+");
}



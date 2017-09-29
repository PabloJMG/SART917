#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <fenv.h>
#include <cuda.h>
#include <cuda_runtime.h>
 #include "device_launch_parameters.h"
//~ #include <helper_cuda.h>
#include "../headers/init_coe.h"
#include "../headers/poinit_coe.h"
#warning "Compilation Error"
//~ #include "../headers/cu_init_coe.h"
#include "../headers/cu_poinit_coe.h"
//Hacer la proyección primero y aplicarla a cada píxel y no hacer la proyección para cada píxel.

//Crear  un array cant con tamaño NumPix*NumAng

//Inicializar coe's con CUDA

// for(NumAng)
// for (NPix)
// for(dimtotal)
//~ {
	
	//~ Suma_proyeccione
	
//~ }

// for(dimtotal)
//~ {
	
	//~ update_pixel
	
//~ }

#define PATH "/home/pablo/local.git/pablo"
#define CUDA 1

#define pi 3.141592
#define DOSPI
//~ #define COEF_COM
#define GUARDARC
#define PARALLEL_BEAM // =! RAND_ANG O SI?
#define DIMX 200
#define DIMY 200
#define FDIM 1.2
#define NANG 30
#define NPIX 600 //Esto sirve para algo ¿?
#define LIMX (int)(DIMX*FDIM)
#define LIMY (int)(DIMY*FDIM)
#define LIMTOT (int) (LIMX*LIMY)
//~ #define TX 10
//~ #define TY 20
//~ #define LIMTOT 20

#define D2N (int)(DIMX*DIMY*NANG)

//Error checking

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

//~ #define D2N 40
#define ind_coe(a, b, c) ((a*DIMY+b)*NANG +c)  //Aqui da problemas, tienes que pensar bien cómo hacerlo
#define ind_coe2(a,b) (a*LIMY+b) //Para n_gamma
#define ind_coe3(a,b) (a*LIMTOT+b) //Para coe y para gamma al parecer
#define coe_gam(a,b) (a*LIMTOT+b)
void *carga_archivo(int NumAng, int NumPix, float *b);
FILE *apertura_archivo(int contador);
FILE *apertura_archivo2(int cont1, int cont2);
void saltar_com(FILE *fin);

 //~ typedef struct dim3{
	
	//~ int x;
	//~ int y;
//~ }

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
	
	int NumPix, gammcon, cont, cont2, cont3, cont4, cont5, cont6, cont7, dimx, dimy, limx, limy, NumAng, facx, facxt, facy, facyt, i, valu, NumRayos1, ng, dimtotal, nozcoe, suma, ac1, ac2, c1, c2, gc1, gc2, bc1,  d2;
	int aint[3];
	float  factora, deltafantes, deltaang, lambda, deltaf, diferencia, cant, sumpix, factorpos, lol;
	char filename[100];
	FILE *finp, *pipeplot;
	clock_t inicio, final, com, fin;
	time_t archtime;
	float factordim = FDIM;
	int TX=32;
	int TY=32;
	int TXX=25;
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
	//~ double angulos[10];
	finp=fopen("../FP/input/NumAngAbanico.txt", "r");
	fscanf(finp, "%d", &NumAng);
	fclose(finp);
	float *angulos=(float*)calloc(NumAng, sizeof(float));
	finp=fopen("../FP/input/AngulosBundle.txt", "r");
	for(cont=0; cont<NumAng; cont++)
	fscanf(finp, "%f", &angulos[cont]);
	fclose(finp);
	#endif

	NumPix=aint[2];
	
	printf("Antes de la declaracion de structs\n");
	//~ NumAng = aint;

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
	int *nzcoea=(int *)calloc(NumPix*d2, sizeof(int)); //runtime variable ¿?¿?¿?¿?¿?¿? Poss lost valgrind
	float *gamma=(float *)calloc(NumAng*dimtotal, sizeof(float)); //valgrind echoes error
	float *b=(float *)calloc(NumAng*NumPix,sizeof(float)); //Si const b, quitar esto. const float *b=read_input_data(); ¿?
	float *betar=(float *)calloc(NumAng*NumPix, sizeof(float)); 
	float *coe=(float *)malloc(dimtotal*NumRayos1*sizeof(float)); //valgrind error
	float *x = (float *) calloc(dimtotal, sizeof(float)); //Def lost
	float *loq=(float *) calloc(NumAng*NumPix, sizeof(float));
	float *difgama=(float *)calloc(dimtotal, sizeof(float));
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
					b[cont*NumPix+cont2]=log(1+valu)+log(2); //No hace mucha diferencia
					//~ b[cont*NumPix+cont2])log(1+valu);
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
	{
		x[cont]=0.f;
		difgama[cont]=0.f;
	}
	//CUDA inits Antes estaba dentro CHANGE 22-9
	float *coe_out;	
	cudaMalloc(&coe_out, dimtotal*NumRayos1*sizeof(float));
	gpuErrchk( cudaPeekAtLastError() );
	dim3 blockSize2(TXX, TXX);
	int axx=(dimtotal+TXX-1)/(TXX);
	int bxx=(NumRayos1+TXX-1)/TXX;
	dim3 gridSize2(axx, bxx);
	printf("Gridsize = %d, blockSize = %d\n", axx,TXX);
	initialize_coe<<<gridSize2, blockSize2>>>(coe_out, NumRayos1); //El segundo valor está bien puesto ¿?
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk(cudaDeviceSynchronize());
	getchar();
	dim3 blockSize(TX, TY);
	int bx=(NumPix+TX-1)/TX; // (NumAng +1)
	int by=(NumAng+TY)/TY; //Antes ponia NumPix
	printf("bx = %d, by=%d\n", bx, by);
	//~ getchar();
	dim3 gridSize(bx,by);
	switch(CUDA){
		
		case(0):
	#ifdef PARALLEL_BEAM
	
	#ifdef DOSPI
	
	#pragma omp parallel 
	{
		
		#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				//Se entra aqui de momento
				pa_initialize_coef_flo_2pi(cont, NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos,  &coe[0]);

			}
		
	}
	break;
	#endif
	
	#ifdef ABANICO
	
	#pragma omp parallel 
	{
		
			#pragma omp for 
			for(cont=0; cont < NumAng; cont++)  //Este for tiene sentido?
			{
				pa_initialize_coef_flo_abanico(cont, NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos, &angulos[0], &coe[0]);
			}
		
	}
	
	break;
	#endif
	
	
	
	#endif
	


	
	
	#ifdef POINT_SOURCE
	
	#ifdef DOSPI
	
	#pragma omp parallel 
	{
		
		#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				
				po_initialize_coef_flo_2pi(cont, NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos,  coe);
			}
		
	}
	break;
	#endif
	
	#ifdef ABANICO
	
	#pragma omp parallel 
	{
		
		 	#pragma omp for 
			for(cont=0; cont < NumAng; cont++)
			{
				po_initialize_coef_flo_abanico(cont, NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos, angulos, coe);
			}
		
	}
	
	break;
	#endif
	
	
	
	#endif
	
	case(1):
	
	//~ cudaMalloc((void**)&coe_out, dimtotal*NumRayos1*sizeof(float));

	//~ dim3 blockSize2(TXX);
	//~ int axx=((dimtotal*NumRayos1)+TXX-1)/TXX;
	//~ dim3 gridSize2(axx);
	//~ initialize_coe<<<gridSize2, blockSize2>>>(coe_out);
	//~ cudaDeviceSynchronize();
	//~ dim3 blockSize(TX, TY);
	//~ int bx=(NumAng+TX)/TX; // (NumAng +1)
	//~ int by=(NumPix+TY-1)/TY; //Antes ponia NumPix
	//~ printf("bx = %d, by=%d\n", bx, by);
	//~ getchar();
	//~ dim3 gridSize(bx,by);
	//~ gpuErrchk(cudaMemcpy(coe_out, coe, dimtotal*NumRayos1*sizeof(float), cudaMemcpyHostToDevice)); //CHANGE 22-9 CHANGE 28-9 Lo he quitado, no sabía qué pintaba aquí
	#ifdef PARALLEL_BEAM
	
	#ifdef DOSPI
	cu_pa_initialize_coef_flo_2pi<<<gridSize, blockSize>>>(NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos, coe_out);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	//~ cu_pa_initialize_coef_flo_2pi<<<(1,1), (NumAng, NumPix)>>>(NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos, coe_out);
	//~ cudaDeviceSynchronize();
	printf("Salio del Kernel\n");

	cudaMemcpy(coe, coe_out, dimtotal*NumRayos1*sizeof(float), cudaMemcpyDeviceToHost);	
	
	for(i=0; i<(dimtotal); i++)
	printf("coe_out[%d] = %f\n", i, coe[i]);
	getchar();
	cudaFree(coe_out);
	break;
	#endif
	
	#ifdef ABANICO
	cu_pa_initialize_coef_flo_abanico<<<gridSize, blockSize>>>(NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, factorpos, angulos, coe_out);
	cudaMemcpy(coe, coe_out, dimtotal*NumRayos1*sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(coe_out);
	break;
	#endif
	#endif
	#ifdef POINT_SOURCE
	
	
	
	#ifdef DOSPI
	cu_po_initialize_coef_flo_2pi<<<gridSize, blockSize>>>(NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, deltaang, factorpos, coe_out);
	cudaMemcpy(coe, coe_out, dimtotal*NumRayos1*sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(coe_out);
	break;
	#endif
	
	#ifdef ABANICO
	cu_po_initialize_coef_flo_abanico<<<gridSize, blockSize>>>(NumPix, NumAng, factordim, limy, dimx, dimy, dimtotal, NumRayos1, factorpos, angulos, coe_out);
	cudaMemcpy(coe, coe_out, dimtotal*NumRayos1*sizeof(float), cudaMemcpyDeviceToHost);
	cudaFree(coe_out);
	break;
	#endif
	
	#endif
	
}
	
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
	
	//while((deltaf>0.1) && (i<50) && (deltafantes>deltaf))
	//{
		//deltafantes=deltaf;
		//deltaf=0;
		//sumpix=0;
		
		//for(cont2=facx; cont2<facyt; cont2++)
		//{
			
			//for(cont3=facy; cont3<facyt; cont3++)
			//{
				
				//difgam=0;
					////~ com=clock();
				//for(cont=0; cont<no_gam[ind_coe2(cont2, cont3)]; cont++)
				//{
					//diferencia=0;
					
					//#pragma omp parallel for private(cont4, cant, cont5, cont6, bc1, bc2, c1, c2, c3, lol, ac1, ac2) reduction(+:diferencia) shared(coe)
					
					//for(cont4=0; cont4<nzcoe[ind_coe((cont2-facx),(cont3-facy),cont)]; cont4++)
					//{
						
						//cant=0;
						//c1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont]*NumPix+nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
						//c2=cont2*limy+cont3;
						//lol=coe[ind_coe3(c1,c2)];
						
						//if(lol>0.f){
							
					
						
						////~ printf("Dentro bucle y cont = %d\n", cont);
						
						
						
							////~ printf("Pro es %d\n", pro);
							////~ getchar();
					
						//for(cont5=facx; cont5<facxt; cont5++)
						//{
							
							//for(cont6=facy; cont6<facyt; cont6++)
							//{
								//ac1=cont*NumPix +nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
								//ac2=cont5*limy+cont6;
								
								//if(coe[ind_coe3(ac1,ac2)]>0)
								//cant+=coe[ind_coe3(ac1,ac2)]*x[ac2]; //Aquí sparse_matrix_multip ¿?
								
								
							//}
							
							
						//}
		
						
						
						//bc1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont];
						//bc2=nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
						//bc1=bc1*NumPix+bc2;
						//c3=no_gama[ind_coe2(cont2, cont3)*NumPix+cont]*NumPix+nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
						//diferencia+=(b[c3]-cant)*(coe[ind_coe3(c1,c2)])/(betar[bc1]);  //Calculo de la diferencia para cada rayo entre la proyección y lo calculado
						
					//}
					
					//}
					
					
					//gc1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont];
					//gc2=cont2*limy+cont3;
					//difgam+=diferencia/(gamma[ind_coe3(gc1, gc2)]); //Esto se podría hacer antes... Se calcula la suma de las diferencias para cada rayo
				//}
				////~ fin=clock();
				////~ printf("Update de pixel tarda %lf\n",(double)(fin-com)/(CLOCKS_PER_SEC));
				////~ getchar();
				//x[cont2*limy+cont3]+=lambda*difgam; //update del pixel
				//deltaf+=lambda*difgam; //variables de control de la convergencia
				//sumpix+=x[cont2*limy+cont3]; // ¿?
				 
			//}  //Final bucle update de un pixel
			
			
			
		//}
		
		//deltaf=deltaf/sumpix;
		//deltaf=fabs(deltaf);
		//final=clock();
		//printf("Deltaf=%lf y n de iteraciones =%d. Desde el inicio llevamos %lf segundos\n", deltaf, i,(double)(final-inicio)/(CLOCKS_PER_SEC*12));
        //sprintf(filename, "../res/Resultados_SART/Resultados_SART_%d.txt", i);

		//finp=fopen(filename, "w+");
		
		//for(cont=0; cont<limx; cont++)
		//{
			//for(cont2=limy; cont2>0; cont2--)
				//fprintf(finp, "%d %d %f\n", cont, limy-cont2, x[cont*limy+cont2]);
			
		//}
		
		//fclose(finp);

        //fprintf(pipeplot, "plot '%s/res/Resultados_SART/Resultados_SART_%d.txt'w image \n",PATH, i);
        //fprintf(pipeplot, "set term postscript \n");

        //fprintf(pipeplot, "set output 'img/imSART/NumAng_%d_Numpix_%d,dim_%dx%d__%d.png' \n", NumAng, NumPix, dimx, dimy, archtime);
        //fprintf(pipeplot, "replot\n");
        //fflush(pipeplot);
        //i+=1;

	//}


while((deltaf>0.01) && (i<100)&& (deltafantes>deltaf) ) //
	{
		deltafantes=deltaf;
		deltaf=0;
		sumpix=0;
		
	//~ com=clock();
				for(cont=0; cont<NumAng; cont++)
				{
					diferencia=0;
					
					#pragma omp parallel for private(cont4, cant, cont5, cont6, bc1,  c1, c2,  lol, ac1, ac2) reduction(+:diferencia) shared(coe)
					
					for(cont4=0; cont4<NumPix; cont4++)
					{
					
					
						
						
						
						
						
							//~ printf("Pro es %d\n", pro);
							//~ getchar();
					
						for(cont5=facx; cont5<facxt; cont5++)
						{
							
							for(cont6=facy; cont6<facyt; cont6++)
							{
								ac1=cont*NumPix+cont4;
								ac2=cont5*limy+cont6;
								
								if(coe[ind_coe3(ac1,ac2)]>0)
								loq[ac1]+=coe[ind_coe3(ac1,ac2)]*x[ac2]; //Aquí sparse_matrix_multip ¿?
								
								
							}
							
							
							
							
						}
						
					
					
				}
				
			}
			cont7=0;
			
			//Calculo de factor de corrección para cada pixel
			for(cont2=facx; cont2<facyt; cont2++)
			{
					for(cont3=facy; cont3<facyt; cont3++)
					{	c2=cont2*limy+cont3;
							for(cont=0; cont<no_gam[ind_coe2(cont2, cont3)]; cont++)
							{
								diferencia=0;
								bc1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont];
								
								
								for(cont4=0; cont4<nzcoe[ind_coe((cont2-facx),(cont3-facy),cont)]; cont4++)
								{
									
									c1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont]*NumPix+nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
								
									lol=coe[ind_coe3(c1,c2)];
									
									
									//~ bc2=nzcoea[ind_coe((cont2-facx),(cont3-facy),cont)*NumPix+cont4];
									bc1++;

									diferencia+=(b[c1]-loq[bc1])*lol/(betar[bc1]); 
									
								
					
								}
					
					
						gc1=no_gama[ind_coe2(cont2, cont3)*NumPix+cont];
						gc2=cont2*limy+cont3;
						difgama[gc2]+=diferencia/(gamma[ind_coe3(gc1, gc2)]);
						}
						
						//~ printf("Hecho para pixel %d\n", cont7);
						cont7++;
					}
				
			}
			
			
			cont7=0;
			
			//update de cada pixel
		for(cont2=facx; cont2<facyt; cont2++)
		{
			
			for(cont3=facy; cont3<facyt; cont3++)
			{
				
				

				//~ printf("Update de pixel %d de %d tarda %lf y es %f\n", cont7, (facyt-facy)*(facxt-facx), (double)(fin-com)/(CLOCKS_PER_SEC), lambda*difgama[cont2*limy+cont3]);
				cont7++;
				//~ getchar();
				x[cont2*limy+cont3]+=lambda*difgama[cont2*limy+cont3];
				deltaf+=lambda*difgama[cont2*limy+cont3];
				sumpix+=x[cont2*limy+cont3];
				 
			}  //Final bucle update de un pixel
			
			
			
		}
		deltaf=deltaf/sumpix;
		deltaf=fabs(deltaf);
		final=clock();
		printf("Deltaf=%lf y n de iteraciones =%d. Desde el inicio llevamos %lf segundos\n", deltaf, i,(double)(final-inicio)/(CLOCKS_PER_SEC*12));
		//~ getchar();
        sprintf(filename, "../res/Resultados_SART/Resultados_SART_%d.txt", i);

		finp=fopen(filename, "w+");
		
		for(cont=0; cont<limx; cont++)
		{
			for(cont2=limy; cont2>0; cont2--)
			
			
				fprintf(finp, "%d %d %f\n", cont, limy-cont2, x[cont*limy+cont2]);  
			
		}
		
		fclose(finp);

       
        i+=1;

	} 
	
	fprintf(pipeplot, "plot '%s/res/Resultados_SART/Resultados_SART_%d.txt'w image \n",PATH, (i-1));
        fprintf(pipeplot, "set term postscript \n");

        fprintf(pipeplot, "set output 'img/imSART/NumAng_%d_Numpix_%d,dim_%dx%d__%d.png' \n", NumAng, NumPix, dimx, dimy, archtime);
        fprintf(pipeplot, "replot\n");
       
	fflush(pipeplot);
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
	free(difgama);
	free(coe);
	free(b);
	free(nzcoea);
	free(x);
	free(betar);
	free(loq);
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



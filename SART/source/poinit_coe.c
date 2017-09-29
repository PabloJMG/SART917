//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>


//#include "../headers/poinit_coe.h"

//void po_initialize_coef_flo_2pi( int cont, int NumPix, int NumAng, float factordim,int limy, int dimx, int dimy, float deltaang, float factorpos, float **coe){


	//int cont2, cont3, pixx, pixy;
	//float cosbeta, sinbeta, beta, deltax, deltay, posx, posy, posx2, posy2;

	////~ #ifdef PARALLEL_BEAM
	////~ #pragma omp parallel 
	////~ {
		////~ #pragma omp for private(cont2, cont3, pixx, pixy, cosbeta, sinbeta, beta, deltax, deltay, posx, posy, posx2, posy2)
			////~ for(cont=0; cont<NumAng; cont++)
			////~ {

				//beta = cont*deltaang;
				//float factora = 1000;
				////~ #ifdef ABANICO
				////~ beta=angulos[cont];
				////~ #endif 
				//cosbeta=cos(beta);
				//sinbeta=sin(beta);
				
				//for(cont2=0; cont2<NumPix; cont2++)
				//{
					
					//for(cont3=0; cont3<dimx; cont3++)
					//{
						
						
						//posx=dimx/2-cont3;
						//posy=dimy/2-cont2*factorpos+(cont2*factorpos-dimy/2)*cont3/dimx; // ¿?¿?¿?
						//posx2=posx;
						//posy2=posy;
						//posx=posx2*cosbeta - posy2*sinbeta; /*inventada a toda prisa*/
						//posy=posy2*cosbeta + posx2*sinbeta;
						 
						//posx+=factordim*dimx/2;
						//posy+=factordim*dimy/2;
						//pixx=trunc(posx);
						//pixy=trunc(posy);
						//deltax=posx-pixx;
						//deltay=posy-pixy;
						
						
						//coe[cont*NumPix+cont2][pixx*limy+pixy]+=(2-deltax-deltay)/factora;
					    //coe[cont*NumPix+cont2][pixx*limy+pixy+1]+=(1-deltax+deltay)/factora;
					    //coe[cont*NumPix+cont2][(pixx+1)*limy+pixy]+=(1-deltay+deltax)/factora;
					    //coe[cont*NumPix+cont2][(pixx+1)*limy+pixy+1]+=(deltax+deltay)/factora;
					    
						//coe[NumAng*NumPix][pixx*limy+pixy]+=coe[cont*NumPix+cont2][pixx*limy+pixy]; /*aqui*/
						//coe[NumAng*NumPix][pixx*limy+pixy+1]+=coe[cont*NumPix+cont2][pixx*limy+pixy+1];
						//coe[NumAng*NumPix][(pixx+1)*limy+pixy]+=coe[cont*NumPix+cont2][(pixx+1)*limy+pixy];
						//coe[NumAng*NumPix][(pixx+1)*limy+pixy+1]+=coe[cont*NumPix+cont2][(pixx+1)*limy+pixy+1];


						
					//}
					
					
					
				//}
				
				
			//}
		
		
		
		
	//void po_initialize_coef_flo_abanico( int cont,int NumPix, int NumAng, float factordim, int limy, int dimx, int dimy, float deltaang, float factorpos, float *angulos, float **coe){	
	

	//int cont2, cont3, pixx, pixy;
	//float cosbeta, sinbeta, beta, deltax, deltay, posx, posy, posx2, posy2;


	//float factora = 1000;
		
				//beta=angulos[cont];
				 
				//cosbeta=cos(beta);
				//sinbeta=sin(beta);
				
				//for(cont2=0; cont2<NumPix; cont2++)
				//{
					
					//for(cont3=0; cont3<dimx; cont3++)
					//{
						
						
							//posx=dimx/2-cont3;
						//posy=dimy/2-cont2*factorpos+(cont2*factorpos-dimy/2)*cont3/dimx; // ¿?¿?¿?
						//posx2=posx;
						//posy2=posy;
						//posx=posx2*cosbeta - posy2*sinbeta; /*inventada a toda prisa*/
						//posy=posy2*cosbeta + posx2*sinbeta;
						////~ posx=dimx/2-cont3;
						////~ posy=dimy/2-cont2*factorpos;
						////~ posx2=posx;
						////~ posy2=posy;
						////~ posx=posx2*cosbeta-posy2*sinbeta;
						////~ posy=posy2*cosbeta+posx2*sinbeta;
						//posx+=factordim*dimx/2;
						//posy+=factordim*dimy/2;
						//pixx=trunc(posx);
						//pixy=trunc(posy);
						//deltax=posx-pixx;
						//deltay=posy-pixy;
						//coe[cont*NumPix+cont2][pixx*limy+pixy]+=(2-deltax-deltay)/factora;
					    //coe[cont*NumPix+cont2][pixx*limy+pixy+1]+=(1-deltax+deltay)/factora;
					    //coe[cont*NumPix+cont2][(pixx+1)*limy+pixy]+=(1-deltay+deltax)/factora;
					    //coe[cont*NumPix+cont2][(pixx+1)*limy+pixy+1]+=(deltax+deltay)/factora;

						//coe[NumAng*NumPix][pixx*limy+pixy]+=coe[cont*NumPix+cont2][pixx*limy+pixy]; /*aqui*/
						//coe[NumAng*NumPix][pixx*limy+pixy+1]+=coe[cont*NumPix+cont2][pixx*limy+pixy+1];
						//coe[NumAng*NumPix][(pixx+1)*limy+pixy]+=coe[cont*NumPix+cont2][(pixx+1)*limy+pixy];
						//coe[NumAng*NumPix][(pixx+1)*limy+pixy+1]+=coe[cont*NumPix+cont2][(pixx+1)*limy+pixy+1];


						
					//}
					
					
					
				//}
				
				
//}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "../headers/poinit_coe.h"

void po_initialize_coef_flo_2pi( int cont, int NumPix, int NumAng, float factordim,int limy, int dimx, int dimy, int dimtot, int NRay, float deltaang, float factorpos, float *coe){


	int cont2, cont3, pixx, pixy, c1, c2a, c2b, c2c, c2d, c1max;
	float cosbeta, sinbeta, beta, deltax, deltay, posx, posy, posx2, posy2;

	//~ #ifdef PARALLEL_BEAM
	//~ #pragma omp parallel 
	//~ {
		//~ #pragma omp for private(cont2, cont3, pixx, pixy, cosbeta, sinbeta, beta, deltax, deltay, posx, posy, posx2, posy2)
			//~ for(cont=0; cont<NumAng; cont++)
			//~ {

				beta = cont*deltaang;
				float factora = 1000;
				//~ #ifdef ABANICO
				//~ beta=angulos[cont];
				//~ #endif 
				cosbeta=cos(beta);
				sinbeta=sin(beta);
				
				for(cont2=0; cont2<NumPix; cont2++)
				{
					
					for(cont3=0; cont3<dimx; cont3++)
					{
						
						posx=dimx/2-cont3;
						posy=dimy/2-cont2*factorpos+(cont2*factorpos-dimy/2)*cont3/dimx; // ¿?¿?¿?
						posx2=posx;
						posy2=posy;
						posx=posx2*cosbeta - posy2*sinbeta; /*inventada a toda prisa*/
						posy=posy2*cosbeta + posx2*sinbeta;
						 
						posx+=factordim*dimx/2;
						posy+=factordim*dimy/2;
						pixx=trunc(posx);
						pixy=trunc(posy);
						deltax=posx-pixx;
						deltay=posy-pixy;
						
						c1=cont*NumPix+cont2;
						c2a=pixx*limy+pixy;
						c2b=pixx*limy+pixy+1;
						c2c=(pixx+1)*limy+pixy;
						c2d=(pixx+1)*limy+pixy+1;
						c1max=NumAng*NumPix;
						
						coe[c1*dimtot+c2a]+=(2-deltax-deltay)/factora;
					    coe[c1*dimtot+c2b]+=(1-deltax+deltay)/factora;
					    coe[c1*dimtot+c2c]+=(1-deltay+deltax)/factora;
					    coe[c1*dimtot+c2d]+=(deltax+deltay)/factora;
					    
						coe[c1max*dimtot+c2a]+=coe[c1*dimtot+c2a]; /*aqui*/
						coe[c1max*dimtot+c2b]+=coe[c1*dimtot+c2b];
						coe[c1max*dimtot+c2c]+=coe[c1*dimtot+c2c];
						coe[c1max*dimtot+c2d]+=coe[c1*dimtot+c2d];


						
					}
					
					
					
				}
				
				
			}
		
		
		
		
	void po_initialize_coef_flo_abanico( int cont,int NumPix, int NumAng, float factordim, int limy, int dimx, int dimy, int dimtot, int NRay,float deltaang, float factorpos, float *angulos, float *coe){	
	

	int cont2, cont3, pixx, pixy,c1, c2a, c2b, c2c, c2d, c1max;
	float cosbeta, sinbeta, beta, deltax, deltay, posx, posy, posx2, posy2;


	float factora = 1000;
		
				beta=angulos[cont];
				 
				cosbeta=cos(beta);
				sinbeta=sin(beta);
				
				for(cont2=0; cont2<NumPix; cont2++)
				{
					
					for(cont3=0; cont3<dimx; cont3++)
					{
						
						c1=cont*NumPix+cont2;
						posx=dimx/2-cont3;
						posy=dimy/2-cont2*factorpos+(cont2*factorpos-dimy/2)*cont3/dimx; // ¿?¿?¿?
						posx2=posx;
						posy2=posy;
						posx=posx2*cosbeta - posy2*sinbeta; /*inventada a toda prisa*/
						posy=posy2*cosbeta + posx2*sinbeta;
						//~ posx=dimx/2-cont3;
						//~ posy=dimy/2-cont2*factorpos;
						//~ posx2=posx;
						//~ posy2=posy;
						//~ posx=posx2*cosbeta-posy2*sinbeta;
						//~ posy=posy2*cosbeta+posx2*sinbeta;
						posx+=factordim*dimx/2;
						posy+=factordim*dimy/2;
						pixx=trunc(posx);
						pixy=trunc(posy);
						deltax=posx-pixx;
						deltay=posy-pixy;
						
						c2a=pixx*limy+pixy;
						c2b=pixx*limy+pixy+1;
						c2c=(pixx+1)*limy+pixy;
						c2d=(pixx+1)*limy+pixy+1;
						c1max=NumAng*NumPix;
						coe[c1*dimtot+c2a]+=(2-deltax-deltay)/factora;
					    coe[c1*dimtot+c2b]+=(1-deltax+deltay)/factora;
					    coe[c1*dimtot+c2c]+=(1-deltay+deltax)/factora;
					    coe[c1*dimtot+c2d]+=(deltax+deltay)/factora;
					    
						coe[c1max*dimtot+c2a]+=coe[c1*dimtot+c2a]; /*aqui*/
						coe[c1max*dimtot+c2b]+=coe[c1*dimtot+c2b];
						coe[c1max*dimtot+c2c]+=coe[c1*dimtot+c2c];
						coe[c1max*dimtot+c2d]+=coe[c1*dimtot+c2d];



						
					}
					
					
					
				}
				
				
}




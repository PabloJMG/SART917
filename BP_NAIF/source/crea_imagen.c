#include "/home/pablo/Resultados_C/LFT/headers/funciones.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
void *crea_imagen( int ray, int numan)
{ 
	FILE *pipeplot=popen("gnuplot -persist", "w");
/*fprintf(pipeplot, "set title /lo que sea/"   */
fprintf(pipeplot, "set pm3d map\n");
fprintf(pipeplot, "set size square 1, 1\n");
fprintf(pipeplot, "splot '/home/pablo/Resultados_C/LFT/res/Resultados_BP.txt' matrix\n");
sleep(5);
fprintf(pipeplot, "set output 'BP_Imagen_Angulos%d_Rayos%d.ps'\n", numan, ray);
/*fprintf(pipeplot, "replot\n");*/

fflush(pipeplot);
fprintf(pipeplot, "quit\n");
pclose(pipeplot);
	
	return 0;
	
} 


//~ FILE *apertura_archivo2(int contador, int contador2, int contador3) 
//~ {

//~ char nombrearchivo[25]="Resultados_pixeles.txt";

//~ sprintf(nombrearchivo, "rayo%d_%d_%d.txt", contador3, contador, contador2);

//~ return fopen(nombrearchivo, "w");



//~ }

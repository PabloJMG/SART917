__global__ void cu_po_initialize_coef_flo_2pi( int NumPix, int NumAng, float factordim, int limy, int dimx, int dimy, int dimtot, int NRay, float deltaang, float factorpos, float *coe);

__global__ void cu_po_initialize_coef_flo_abanico(int NumPix, int NumAng, float factordim, int limy, int dimx, int dimy, int dimtot, int NRay,  float factorpos, float *angulos, float *coe);

__global__ void cu_pa_initialize_coef_flo_2pi(  int NumPix,int NumAng, float factordim, int limy, int dimx, int dimy, int dimtot, int NRay, float deltaang, float factorpos, float *coe);

__global__ void cu_pa_initialize_coef_flo_abanico(  int NumPix,int NumAng, float factordim, int limy, int dimx, int dimy, int dimtot, int NRay, float factorpos, float *angulos, float *coe);

__global__ void initialize_coe(float *coe, int T);

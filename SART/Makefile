NVCC = /usr/local/cuda/bin/nvcc
NVCC_FLAGS = -g -arch sm_35 -G -Xcompiler -Wall

all: exe/main.exe

exe/main.exe: main.o poinit_coe.o  
	$(NVCC) $^ -o $@

main.o: source/camcam.cu headers/cu_poinit_coe.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

poinit_coe.o: source/cu_poinit_coe.cu headers/cu_poinit_coe.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@




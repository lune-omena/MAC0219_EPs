all:    mandelbrot_seq.c mandelbrot_pth.c mandelbrot_omp.c mandelbrot_cuda.cu mandelbrot_ompi.c mandelbrot_ompi_omp.c
	gcc -o mandelbrot_seq mandelbrot_seq.c
	gcc -pthread -o mandelbrot_pth mandelbrot_pth.c
	gcc -o mandelbrot_omp -fopenmp mandelbrot_omp.c
	nvcc -o mandelbrot_cuda mandelbrot_cuda.cu
	mpicc mandelbrot_ompi.c -o mandelbrot_ompi
	mpicc -fopenmp mandelbrot_ompi_omp.c -o mandelbrot_ompi_omp

mandelbrot_seq: mandelbrot_seq.c
	gcc -o mandelbrot_seq mandelbrot_seq.c
	
mandelbrot_pth: mandelbrot_pth.c
	gcc -pthread -o mandelbrot_pth mandelbrot_pth.c

mandelbrot_omp: mandelbrot_omp.c
	gcc -o mandelbrot_omp -fopenmp mandelbrot_omp.c
	
mandelbrot_cuda: mandelbrot_cuda.cu
	nvcc -o mandelbrot_cuda mandelbrot_cuda.cu
	
mandelbrot_ompi: mandelbrot_ompi.c
	mpicc mandelbrot_ompi.c -o mandelbrot_ompi
	
mandelbrot_ompi_omp: mandelbrot_ompi_omp.c
	mpicc -fopenmp mandelbrot_ompi_omp.c -o mandelbrot_ompi_omp
		




mandelbort_seq: mandelbrot_seq.c
	gcc -o mandelbrot_seq mandelbrot_seq.c

mandelbrot_pth: mandelbrot_pth.c
	gcc -pthread -o mandelbrot_pth mandelbrot_pth.c

mandelbrot_omp: mandelbrot_omp.c
	gcc -o mandelbrot_omp -fopenmp mandelbrot_omp.c

mandelbrot_cuda: mandelbrot_cuda.cu
	nvcc -o mandelbrot_cuda mandelbrot_cuda.cu
	./mandelbrot_cuda 0.175 0.375 -0.1 0.1 200 0


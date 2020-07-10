mandelbrot_seq: mandelbrot_seq.c
	gcc -o mandelbrot_seq mandelbrot_seq.c
	./mandelbrot_seq 0.175 0.375 -0.1 0.1 200 0
mandelbrot_pth: mandelbrot_pth.c
	gcc -pthread -o mandelbrot_pth mandelbrot_pth.c

mandelbrot_omp: mandelbrot_omp.c
	gcc -o mandelbrot_omp -fopenmp mandelbrot_omp.c
	
mandelbrot_cuda: mandelbrot_cuda.cu
	nvcc -o mandelbrot_cuda mandelbrot_cuda.cu
	
mandelbrot_ompi: mandelbrot_ompi.c
	mpicc mandelbrot_ompi.c -o mandelbrot_ompi
	mpirun -n 64 --host ./mandelbrot_ompi 0.175 0.375 -0.1 0.1 200 
	
exec:
	./mandelbrot_cuda -2.5 1.5 -2.0 2.0 3096 4 16
clean:
	rm mandelbrot_cuda



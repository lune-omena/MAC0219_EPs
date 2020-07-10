#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include<cuda.h>
#include<assert.h>

struct timer_info {
    clock_t c_start;
    clock_t c_end;
    struct timespec t_start;
    struct timespec t_end;
    struct timeval v_start;
    struct timeval v_end;
};

struct timer_info timer, timerAloc;

void write_to_file(int width,int height, unsigned char * image_buffer){
    FILE * file;
    const char * filename               = "output_cuda.ppm";
    const char * comment                = "# ";
    int max_color_component_value = 255;
    unsigned char* ch=(unsigned char*)malloc(3*sizeof(unsigned char));

    file = fopen(filename,"wb");
    
    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            width, height, max_color_component_value);
    for(int i = 0; i < 3*width*height; i+=3){
        ch[0]=image_buffer[i];
        ch[1]=image_buffer[i+1];
        ch[2]=image_buffer[i+2];
        fwrite(ch, 1 , 3, file);
    };
    fclose(file);
};

__global__ void compute_mandelbrot(double c_x_min, double c_x_max, double c_y_min, double c_y_max, unsigned char * image_buffer,\
     int* colors, double pixel_height,double pixel_width, int i_x_max, int i_y_max,int tasks_per_thread){


    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;
    int iteration_max=200;
    int gradient_size=16;
    double c_x;
    double c_y;

    int iteration;
    unsigned int i_yo = tasks_per_thread*(blockIdx.x * blockDim.x + threadIdx.x);
    unsigned int i_xo = tasks_per_thread*(blockIdx.y * blockDim.y + threadIdx.y);
    unsigned int i_y,i_x;
    for(int i=0;i<tasks_per_thread;i++){
        i_y=i_yo+i;
        for(int j=0;j<tasks_per_thread;j++){
            i_x=i_xo+j;
            if(i_y<i_y_max && i_x<i_x_max){
                c_y = c_y_min + i_y * pixel_height;
                if(c_y < pixel_height / 2 && c_y > -pixel_height / 2){
                    c_y = 0.0;
                };

                c_x         = c_x_min + i_x * pixel_width;

                z_x         = 0.0;
                z_y         = 0.0;

                z_x_squared = 0.0;
                z_y_squared = 0.0;
                
                for(iteration = 0;
                    iteration < iteration_max && \
                    ((z_x_squared + z_y_squared) < escape_radius_squared);
                    iteration++){
                    
                    z_y         = 2 * z_x * z_y + c_y;
                    z_x         = z_x_squared - z_y_squared + c_x;

                    z_x_squared = z_x * z_x;
                    z_y_squared = z_y * z_y;
                };
                if(iteration == iteration_max){
                    image_buffer[((i_y_max * i_y) + i_x)*3+0] = colors[gradient_size*3+0];
                    image_buffer[((i_y_max * i_y) + i_x)*3+1] = colors[gradient_size*3+1];
                    image_buffer[((i_y_max * i_y) + i_x)*3+2] = colors[gradient_size*3+2];

                }
                else{
                    int color = iteration % gradient_size;
            
                    image_buffer[((i_y_max * i_y) + i_x)*3+0] = colors[color*3+0];
                    image_buffer[((i_y_max * i_y) + i_x)*3+1] = colors[color*3+1];
                    image_buffer[((i_y_max * i_y) + i_x)*3+2] = colors[color*3+2];
                };
            }
        }
    }
};


int main(int argc, char *argv[]){

    double c_x_min;
    double c_x_max;
    double c_y_min;
    double c_y_max;
    double pixel_width;
    double pixel_height;
    int image_size;
    int i_x_max;
    int i_y_max;
    int image_buffer_size;
    unsigned char *image_buffer;
    int *colors;
    int threads_per_block;
    int numblock;
    int tasks_per_thread;
    int tot_threads;

    timerAloc.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_start);
    gettimeofday(&timerAloc.v_start, NULL);

    if(argc < 6){
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq  \n");
        exit(0);
    } 
    else{
        sscanf(argv[1], "%lf", &(c_x_min));
        sscanf(argv[2], "%lf", &(c_x_max));
        sscanf(argv[3], "%lf", &(c_y_min));
        sscanf(argv[4], "%lf", &(c_y_max));
        sscanf(argv[5], "%d", &(image_size));
        sscanf(argv[6],"%d",&numblock);
        sscanf(argv[7],"%d", &threads_per_block);
        i_x_max           = image_size;
        i_y_max           = image_size;
        image_buffer_size = image_size * image_size;
    
        pixel_width       = (double)(c_x_max - c_x_min) / i_x_max;
        pixel_height      = (double)(c_y_max - c_y_min) / i_y_max;
    };
    tot_threads=threads_per_block*numblock;
    tasks_per_thread=image_size/tot_threads;
    if(tasks_per_thread<1)tasks_per_thread=1;

    cudaMallocManaged(&(colors),51*sizeof(int));
    int colorvals[51]={
    66, 30, 15,
    25, 7, 26,
    9, 1, 47,
    4, 4, 73,
    0, 7, 100,
    12, 44, 138,
    24, 82, 177,
    57, 125, 209,
    134, 181, 229,
    211, 236, 248,
    241, 233, 191,
    248, 201, 95,
    255, 170, 0,
    204, 128, 0,
    153, 87, 0,
    106, 52, 3,
    16, 16, 16,
    };
    for(int i=0;i<51;i++){
        colors[i]=colorvals[i];
    }
    cudaMalloc((void **) &image_buffer, sizeof(unsigned char)*image_buffer_size*3);


    timer.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
    gettimeofday(&timer.v_start, NULL);
    

    dim3 blockDim(threads_per_block, threads_per_block, 1);
    //dim3 gridDim(i_x_max/ blockDim.x, i_y_max/blockDim.y, 1);
    dim3 gridDim(numblock, numblock, 1);

    compute_mandelbrot<<<gridDim,blockDim,0>>>(c_x_min,c_x_max,c_y_min,c_y_max,image_buffer,colors,pixel_height,pixel_width,i_x_max,i_y_max,tasks_per_thread);    
    

    unsigned char * hostimage_buffer;
    hostimage_buffer=(unsigned char *)malloc(image_buffer_size*3*sizeof(unsigned char));
    cudaMemcpy(hostimage_buffer,image_buffer,sizeof(unsigned char)*image_buffer_size*3,cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();


    timer.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
    gettimeofday(&timer.v_end, NULL);

    timerAloc.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_end);
    gettimeofday(&timerAloc.v_end, NULL);

    write_to_file(i_x_max,i_y_max,hostimage_buffer);
    printf("%f",
        (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
        (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
    
    printf (",");
    printf("%f",
        (double) (timerAloc.t_end.tv_sec - timerAloc.t_start.tv_sec) +
        (double) (timerAloc.t_end.tv_nsec - timerAloc.t_start.tv_nsec) / 1000000000.0);
    printf ("\n");

    cudaFree(colors);
    cudaFree(image_buffer);
    free(hostimage_buffer);

    return 0;
};

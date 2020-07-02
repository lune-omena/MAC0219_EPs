#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include<cuda.h>
#define SIZE 10
#define iteration_max 200
#define gradient_size 16

struct timer_info {
    clock_t c_start;
    clock_t c_end;
    struct timespec t_start;
    struct timespec t_end;
    struct timeval v_start;
    struct timeval v_end;
};

struct timer_info timer, timerAloc;

typedef struct param{
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
} params;

typedef params* Params;

void allocate_image_buffer(Params p){
    int rgb_size = 3;
    cudaMallocManaged(&(p->image_buffer), sizeof(unsigned char) * p->image_buffer_size * rgb_size);
};

Params CreateParams(){
    Params p;
    cudaMallocManaged(&p,sizeof(Params));
    cudaMallocManaged(&(p->colors),51*sizeof(int));
    allocate_image_buffer(p);
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
        p->colors[i]=colorvals[i];
    }
    return p;
}


void init(int argc, char *argv[],Params p){
    if(argc < 6){
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    } 
    else{
        printf("AAAA\n");
        sscanf(argv[1], "%lf", &(p->c_x_min));
        sscanf(argv[2], "%lf", &(p->c_x_max));
        sscanf(argv[3], "%lf", &(p->c_y_min));
        sscanf(argv[4], "%lf", &(p->c_y_max));
        sscanf(argv[5], "%d", &(p->image_size));
        p->i_x_max           = p->image_size;
        p->i_y_max           = p->image_size;
        p->image_buffer_size = p->image_size * p->image_size;
    
        p->pixel_width       = (double)(p->c_x_max - p->c_x_min) / p->i_x_max;
        p->pixel_height      = (double)(p->c_y_max - p->c_y_min) / p->i_y_max;
    };
};
/*
void update_rgb_buffer(int iteration, int x, int y,Params p){
    int color;

    if(iteration == iteration_max){
        image_buffer[(p.i_y_max * y) + x][0] = colors[gradient_size][0];
        image_buffer[(p.i_y_max * y) + x][1] = colors[gradient_size][1];
        image_buffer[(p.i_y_max * y) + x][2] = colors[gradient_size][2];
    }
    else{
        color = iteration % gradient_size;

        image_buffer[(p.i_y_max * y) + x][0] = colors[color][0];
        image_buffer[(p.i_y_max * y) + x][1] = colors[color][1];
        image_buffer[(p.i_y_max * y) + x][2] = colors[color][2];
    };
};
*/
void write_to_file(Params p){
    FILE * file;
    const char * filename               = "output_cuda.ppm";
    const char * comment                = "# ";
    int max_color_component_value = 255;
    char* ch=(char*)malloc(3*sizeof(char));

    file = fopen(filename,"wb");
    
    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            p->i_x_max, p->i_y_max, max_color_component_value);
    printf("CCCCC\n");
    printf("%c\n",p->image_buffer[0]);
    printf("CCCCC\n");
    for(int i = 0; i < 3*p->image_buffer_size; i+=3){
        ch[0]=p->image_buffer[i];
        ch[1]=p->image_buffer[i+1];
        ch[2]=p->image_buffer[i+2];
        fwrite(ch, 1 , 3, file);
    };
    fclose(file);
    printf("CCCCC\n");
};

__global__ void compute_mandelbrot(Params p){


    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;
    double c_x;
    double c_y;

    int iteration;
    uint i_x = (blockIdx.x * blockDim.x) + threadIdx.x;
    uint i_y = (blockIdx.y * blockDim.y) + threadIdx.y;

    printf("%d,%d -",i_x,i_y);

    if(i_y<p->i_y_max){
        c_y = p->c_y_min + i_y * p->pixel_height;
        if(fabs(c_y) < p->pixel_height / 2){
            c_y = 0.0;
        };
        if(i_x<p->i_x_max){

            c_x         = p->c_x_min + i_x * p->pixel_width;

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

            int color;

            if(iteration == iteration_max){ 
                p->image_buffer[((p->i_y_max * i_y) + i_x)*3+0] = p->colors[gradient_size*3+0];
                p->image_buffer[((p->i_y_max * i_y) + i_x)*3+1] = p->colors[gradient_size*3+1];
                p->image_buffer[((p->i_y_max * i_y) + i_x)*3+2] = p->colors[gradient_size*3+2];

            }
            else{
                color = iteration % gradient_size;
        
                p->image_buffer[((p->i_y_max * i_y) + i_x)*3+0] = p->colors[color*3+0];
                p->image_buffer[((p->i_y_max * i_y) + i_x)*3+1] = p->colors[color*3+1];
                p->image_buffer[((p->i_y_max * i_y) + i_x)*3+2] = p->colors[color*3+2];
            };

        }
    }
};

int main(int argc, char *argv[]){
    Params p;
    p=CreateParams();
    timerAloc.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_start);
    gettimeofday(&timerAloc.v_start, NULL);

    init(argc, argv,p);


    timer.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
    gettimeofday(&timer.v_start, NULL);

    dim3 threadsPerBlock(8, 8); 
    dim3 numBlocks(p->image_size/8,p->image_size/8);

    compute_mandelbrot<<<numBlocks,threadsPerBlock>>>(p);    
    cudaDeviceSynchronize();
    printf("BBBBBb\n");
    
    printf("%c ",p->image_buffer[0]);
    printf("BBBBBb\n");
    
    

    timer.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
    gettimeofday(&timer.v_end, NULL);

    timerAloc.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_end);
    gettimeofday(&timerAloc.v_end, NULL);

    write_to_file(p);
    printf("%f",
        (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
        (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
    
    printf (",");
    printf("%f",
        (double) (timerAloc.t_end.tv_sec - timerAloc.t_start.tv_sec) +
        (double) (timerAloc.t_end.tv_nsec - timerAloc.t_start.tv_nsec) / 1000000000.0);
    printf ("\n");

    cudaFree(p->colors);
    cudaFree(p->image_buffer);
    cudaFree(p);

    return 0;
};

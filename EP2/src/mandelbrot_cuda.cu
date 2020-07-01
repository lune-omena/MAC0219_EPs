#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define SIZE 1024

struct timer_info {
    clock_t c_start;
    clock_t c_end;
    struct timespec t_start;
    struct timespec t_end;
    struct timeval v_start;
    struct timeval v_end;
};

struct timer_info timer, timerAloc;

__device__ double c_x_min;
__device__ double c_x_max;
__device__ double c_y_min;
__device__ double c_y_max;

__device__ double pixel_width;
__device__ double pixel_height;

__device__ int iteration_max = 200;

__device__ int image_size;
__device__ unsigned char **image_buffer;

__device__ int i_x_max;
__device__ int i_y_max;
__device__ int image_buffer_size;

__device__ int gradient_size = 16;
__device__ int colors[17][3] = {
                        {66, 30, 15},
                        {25, 7, 26},
                        {9, 1, 47},
                        {4, 4, 73},
                        {0, 7, 100},
                        {12, 44, 138},
                        {24, 82, 177},
                        {57, 125, 209},
                        {134, 181, 229},
                        {211, 236, 248},
                        {241, 233, 191},
                        {248, 201, 95},
                        {255, 170, 0},
                        {204, 128, 0},
                        {153, 87, 0},
                        {106, 52, 3},
                        {16, 16, 16},
                    };

__device__ void allocate_image_buffer(){
    int rgb_size = 3;
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);

    for(int i = 0; i < image_buffer_size; i++){
        image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
    };
};

__host__ void sscanfall(char* argv[]){
    sscanf(argv[1], "%lf", &c_x_min);
    sscanf(argv[2], "%lf", &c_x_max);
    sscanf(argv[3], "%lf", &c_y_min);
    sscanf(argv[4], "%lf", &c_y_max);
    sscanf(argv[5], "%d", &image_size);
}

__device__ void initdevice(){
    i_x_max           = image_size;
    i_y_max           = image_size;
    image_buffer_size = image_size * image_size;

    pixel_width       = (c_x_max - c_x_min) / i_x_max;
    pixel_height      = (c_y_max - c_y_min) / i_y_max;
}

__host__ __device__ void init(int argc, char *argv[]){
    if(argc < 6){
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq -0.188 -0.012 0.554 0.754 11500\n");
        //exit(0);
        return;
    }
    else{
        sscanfall(argv);
        initdevice();
    };
};

__device__ void update_rgb_buffer(int iteration, int x, int y){
    int color;

    if(iteration == iteration_max){
        image_buffer[(i_y_max * y) + x][0] = colors[gradient_size][0];
        image_buffer[(i_y_max * y) + x][1] = colors[gradient_size][1];
        image_buffer[(i_y_max * y) + x][2] = colors[gradient_size][2];
    }
    else{
        color = iteration % gradient_size;

        image_buffer[(i_y_max * y) + x][0] = colors[color][0];
        image_buffer[(i_y_max * y) + x][1] = colors[color][1];
        image_buffer[(i_y_max * y) + x][2] = colors[color][2];
    };
};

__host__ void write_to_file(){
    FILE * file;
    const char * filename               = "output.ppm";
    const char * comment                = "# ";
    int max_color_component_value = 255;
    int i_x_m=i_x_max;

    file = fopen(filename,"wb");
    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            i_x_m, i_y_max, max_color_component_value);

    for(int i = 0; i < image_buffer_size; i++){
        fwrite(image_buffer[i], 1 , 3, file);
    };

    fclose(file);
};

__global__ void compute_mandelbrot(){
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;
    double c_x;
    double c_y;

    int iteration;
    int i_x=threadIdx.x;
    int i_y=threadIdx.y;


    c_y = c_y_min + i_y * pixel_height;
    if(fabs(c_y) < pixel_height / 2){
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

        update_rgb_buffer(iteration, i_x, i_y);
};

int main(int argc, char *argv[]){
    timerAloc.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_start);
    gettimeofday(&timerAloc.v_start, NULL);

    init(argc, argv);

    allocate_image_buffer();

    timer.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
    gettimeofday(&timer.v_start, NULL);

    compute_mandelbrot<<<1,SIZE>>>();
    cudaDeviceSynchronize();

    timer.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
    gettimeofday(&timer.v_end, NULL);

    timerAloc.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_end);
    gettimeofday(&timerAloc.v_end, NULL);

    
    printf("%f",
        (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
        (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
    write_to_file();
    printf (",");
    printf("%f",
        (double) (timerAloc.t_end.tv_sec - timerAloc.t_start.tv_sec) +
        (double) (timerAloc.t_end.tv_nsec - timerAloc.t_start.tv_nsec) / 1000000000.0);
    printf ("\n");
    
    return 0;
};
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

struct timer_info
{
    clock_t c_start;
    clock_t c_end;
    struct timespec t_start;
    struct timespec t_end;
    struct timeval v_start;
    struct timeval v_end;
};

struct timer_info timer, timerAloc;

double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
unsigned char **image_buffer = NULL;

int i_x_max;
int i_y_max;
int image_buffer_size;

int gradient_size = 16;
int colors[17][3] = {
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

void allocate_image_buffer()
{
    int rgb_size = 3;
    image_buffer = (unsigned char **)malloc(sizeof(unsigned char *) * image_buffer_size);

    for (int i = 0; i < image_buffer_size; i++)
    {
        image_buffer[i] = (unsigned char *)malloc(sizeof(unsigned char) * rgb_size);
    };
};

void init(int argc, char *argv[])
{
    if (argc < 6)
    {
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    }
    else
    {
        sscanf(argv[1], "%lf", &c_x_min);
        sscanf(argv[2], "%lf", &c_x_max);
        sscanf(argv[3], "%lf", &c_y_min);
        sscanf(argv[4], "%lf", &c_y_max);
        sscanf(argv[5], "%d", &image_size);

        i_x_max = image_size;
        i_y_max = image_size;
        image_buffer_size = image_size * image_size;

        pixel_width = (c_x_max - c_x_min) / i_x_max;
        pixel_height = (c_y_max - c_y_min) / i_y_max;
    };
};

void update_rgb_buffer(unsigned char **sub_image_buffer, int iteration, int x)
{
    int color;
    if (iteration == iteration_max)
    {

        sub_image_buffer[x][0] = colors[gradient_size][0];
        sub_image_buffer[x][1] = colors[gradient_size][1];
        sub_image_buffer[x][2] = colors[gradient_size][2];
    }
    else
    {
        color = iteration % gradient_size;

        sub_image_buffer[x][0] = colors[color][0];
        sub_image_buffer[x][1] = colors[color][1];
        sub_image_buffer[x][2] = colors[color][2];
    };
};

void write_to_file()
{
    FILE *file;
    char *filename = "output_ompi.ppm";
    char *comment = "# ";

    int max_color_component_value = 255;

    file = fopen(filename, "wb");

    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            i_x_max, i_y_max, max_color_component_value);

    for (int i = 0; i < image_buffer_size; i++)
    {
        fwrite(image_buffer[i], 1, 3, file);
    };

    fclose(file);
};

void compute_mandelbrot(unsigned char **sub_image_buffer, int j)
{
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;
    int i_x;

    int iteration;

    double c_x;
    double c_y;

    c_y = c_y_min + j * pixel_height;
    if (fabs(c_y) < pixel_height / 2)
    {
        c_y = 0.0;
    };
    for (i_x = 0; i_x < i_x_max; i_x++)
    {
        c_x = c_x_min + i_x * pixel_width;

        z_x = 0.0;
        z_y = 0.0;

        z_x_squared = 0.0;
        z_y_squared = 0.0;
        for (iteration = 0;
             iteration < iteration_max &&
             ((z_x_squared + z_y_squared) < escape_radius_squared);
             iteration++)
        {
            z_y = 2 * z_x * z_y + c_y;
            z_x = z_x_squared - z_y_squared + c_x;

            z_x_squared = z_x * z_x;
            z_y_squared = z_y * z_y;
        }
        update_rgb_buffer(sub_image_buffer, iteration, i_x);
    }
}

void RegisterResult(unsigned char **sub_image_buffer, unsigned char **image_buffer, int i)
{
    for (int k = 0; k < i_x_max; k++)
    {
        image_buffer[i_x_max * i + k][0] = sub_image_buffer[k][0];
        image_buffer[i_x_max * i + k][1] = sub_image_buffer[k][1];
        image_buffer[i_x_max * i + k][2] = sub_image_buffer[k][2];
    }
}

unsigned char **alloc_2d_array(int rows, int cols)
{
    unsigned char *data = (unsigned char *)malloc(rows * cols * sizeof(unsigned char));
    unsigned char **array = (unsigned char **)malloc(rows * sizeof(unsigned char *));
    for (int i = 0; i < rows; i++)
        array[i] = &(data[cols * i]);

    return array;
}

int main(int argc, char *argv[])
{
    timerAloc.c_start = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_start);
    gettimeofday(&timerAloc.v_start, NULL);

    init(argc, argv);

    int rank, size, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    int root = 0, i = 0;
    unsigned char **sub_image_buffer;

    MPI_Status status;

    // get rank of this proces
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get total process number
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Get_processor_name(hostname, &h_len);

    //printf("Start! rank:%d size: %d at %s\n", rank, size,hostname);
    sub_image_buffer = alloc_2d_array(i_x_max, 3);

    //o código abaixo roda paralelamenteem multiplos processos
    //variando o rank (numero do processo)

    //processo rank 0 atua como gerenciador
    if (rank == root)
    {
        int i_completed = 0;
        unsigned char *subanswer;
        int CompletedTaskIndex;
        allocate_image_buffer();

        timer.c_start = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
        gettimeofday(&timer.v_start, NULL);

        int *IsWorking = calloc(size, sizeof(int));
        //enquanto todas as rtarefas não forem completas
        while (i_completed < i_y_max)
        {
            //procura por processos ociosos
            for (int j = 1; j < size; j++)
            {
                //se estiver ocioso e ainda houver trabalho para assinalar, assinala
                if (!(IsWorking[j]) && i < i_y_max)
                {
                    //manda o indice i da coluna a ser calculada para o processo j
                    MPI_Send(&i, 1, MPI_INT, j, 0, MPI_COMM_WORLD);
                    IsWorking[j] = 1;
                    i++;
                }
            }
            //espera para receber a resposta de qualquer processo
            MPI_Recv(&(sub_image_buffer[0][0]), 3 * i_x_max, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&CompletedTaskIndex, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //após receber, registra o resultado em image_buffer
            RegisterResult(sub_image_buffer, image_buffer, CompletedTaskIndex);
            //processo que entregou a resposta se torna ocioso
            IsWorking[status.MPI_SOURCE] = 0;
            i_completed++;
        }
        //quando as tarefas acabarem, mandar os processos temrinarem
        for (i = 1; i < size; i++)
        {
            int stop = -1;
            MPI_Send(&stop, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        //cria imagem
        write_to_file();
        for (int i = 0; i < image_buffer_size; i++)
        {
            free(image_buffer[i]);
        }
        free(image_buffer);
        //outros processos além do root
    }
    else
    {
        //processos procurando por trabalho
        while (1)
        {
            //root orienta a calcular a coluna i
            MPI_Recv(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            //se i for -1, o trabalho acabou e pode ser finalizado
            if (i == -1)
                break;
            //calcula
            compute_mandelbrot(sub_image_buffer, i);

            //manda respostas de volta para o root
            MPI_Send(&(sub_image_buffer[0][0]), 3 * i_x_max, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&i, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    free(sub_image_buffer[0]);
    free(sub_image_buffer);

    //printf("End! rank:%d size: %d at %s\n", rank, size,hostname);

    MPI_Finalize();

    timerAloc.c_end = clock();
    clock_gettime(CLOCK_MONOTONIC, &timerAloc.t_end);
    gettimeofday(&timerAloc.v_end, NULL);

    if (rank == root)
    {
        timer.c_end = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
        gettimeofday(&timer.v_end, NULL);
        printf("%f",
               (double)(timer.t_end.tv_sec - timer.t_start.tv_sec) +
                   (double)(timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);

        printf(",");
        printf("%f",
               (double)(timerAloc.t_end.tv_sec - timerAloc.t_start.tv_sec) +
                   (double)(timerAloc.t_end.tv_nsec - timerAloc.t_start.tv_nsec) / 1000000000.0);
        printf("\n");
    }

    return 0;
};

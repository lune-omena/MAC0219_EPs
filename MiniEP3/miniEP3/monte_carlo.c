#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>

#ifndef DEBUG
#define DEBUG 0
#endif

#ifndef VERBOSE
#define VERBOSE 0
#endif

#define FUNCTIONS 1

#define N 10

struct timer_info {
    clock_t c_start;
    clock_t c_end;
    struct timespec t_start;
    struct timespec t_end;
    struct timeval v_start;
    struct timeval v_end;
};

struct timer_info timer;

char *usage_message = "usage: ./monte_carlo SAMPLES FUNCTION_ID N_THREADS\n";

struct function {
    long double (*f)(long double);
    long double interval[2];
};

long double rand_interval[] = {0.0, (long double) RAND_MAX};

long double f1(long double x){
    return 2 / (sqrt(1 - (x * x)));
}

struct function functions[] = {
                               {&f1, {0.0, 1.0}}
};

// Your thread data structures go here

struct thread_data{
    long double *samples;
    long double (*f)(long double s);
    int size;
};

typedef struct thread_data * thread_data_array;

// End of data structures


// começo - funcoes nossas

thread_data_array CriaThread_data(long double (*f)(long double), long double *samples, int size){
    thread_data_array T;
    T=malloc(sizeof(struct thread_data));
    //T->f=malloc(sizeof(long double*));
    T->f=f;
    //T->samples=malloc(sizeof(long double*));
    T->samples=samples;
    T->size=size;
    return T;
}

// termino - funcoes nossas

long double *samples;
long double *results;

long double map_intervals(long double x, long double *interval_from, long double *interval_to){
    x -= interval_from[0];
    x /= (interval_from[1] - interval_from[0]);
    x *= (interval_to[1] - interval_to[0]);
    x += interval_to[0];
    return x;
}

long double *uniform_sample(long double *interval, long double *samples, int size){
    for(int i = 0; i < size; i++){
        samples[i] = map_intervals((long double) rand(),
                                   rand_interval,
                                   interval);
    }
    //printf("%Lf BBBBBBBBB\n",samples[0]);
    return samples;
}

void print_array(long double *sample, int size){
    printf("array of size [%d]: [", size);

    for(int i = 0; i < size; i++){
        printf("%Lf", sample[i]);

        if(i != size - 1){
            printf(", ");
        }
    }

    printf("]\n");
}

long double monte_carlo_integrate(long double (*f)(long double), long double *samples, int size){
    // Your sequential code goes here
    //printf("%Lf",samples[0]);
    //printf("Not implemented yet\n");
    //exit(-1);
    long double soma=0,av;
    for(int i=0;i<size;i++){
        soma+=(*f)(samples[i]);
    }
    av=(long double)soma/size;
    return av; 
}

void *monte_carlo_integrate_thread(void *args){
    // Your pthreads code goes here
    thread_data_array Tdata=(thread_data_array) args;
    printf("%Lf,%Lf, %d, AAAAAAAAAAAAA\n",
        Tdata->f(Tdata->samples[Tdata->size-1]),Tdata->samples[Tdata->size-1],Tdata->size);
    long double* soma=malloc(sizeof(long double)*1);
    *soma=0;
    for (int i = 0; i < Tdata->size; i++)
    {
        (*soma)=(*soma)+Tdata->f(samples[i]);
    }
    return soma;
}

int main(int argc, char **argv){
    if(argc != 4){
        printf("%s",usage_message);
        exit(-1);
    } else if(atoi(argv[2]) >= FUNCTIONS || atoi(argv[2]) < 0){
        printf("Error: FUNCTION_ID must in [0,%d]\n", FUNCTIONS - 1);
        printf("%s",usage_message);
        exit(-1);
    } else if(atoi(argv[3]) < 0){
        printf("Error: I need at least 1 thread\n");
        printf("%s",usage_message);
        exit(-1);
    }

    if(DEBUG){
        printf("Running on: [debug mode]\n");
        printf("Samples: [%s]\n", argv[1]);
        printf("Function id: [%s]\n", argv[2]);
        printf("Threads: [%s]\n", argv[3]);
        printf("Array size on memory: [%.2LFGB]\n", ((long double) atoi(argv[1]) * sizeof(long double)) / 1000000000.0);
    }

    srand(time(NULL));

    int size = atoi(argv[1]);
    struct function target_function = functions[atoi(argv[2])];
    int n_threads = atoi(argv[3]);

    samples = malloc(size * sizeof(long double));

    long double estimate;

    if(n_threads == 1){
        if(DEBUG){
            printf("Running sequential version\n");
        }

        timer.c_start = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
        gettimeofday(&timer.v_start, NULL);

        estimate = monte_carlo_integrate(target_function.f,
                                         uniform_sample(target_function.interval,
                                                        samples,
                                                        size),
                                         size);

        timer.c_end = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
        gettimeofday(&timer.v_end, NULL);
    } else {
        if(DEBUG){
            printf("Running parallel version\n");
        }

        timer.c_start = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_start);
        gettimeofday(&timer.v_start, NULL);

        // Your pthreads code goes here
        int partsize=(int)(size/N);
        long double **parts,*totsample=uniform_sample(target_function.interval,samples,size);
        parts=malloc(N*sizeof(long double*));
        for (int i = 0; i < N; i++)
        {
            parts[i]=malloc(partsize*sizeof(long double));
            for (int j = 0; j < partsize;j++)
            {
                parts[i][j]=totsample[partsize*i+j];
            }
        }
                          
        thread_data_array *Tdatas;
        Tdatas=malloc(N*sizeof(thread_data_array));
        
        for (int i = 0; i < N; i++)
        {
            Tdatas[i]=CriaThread_data(target_function.f,
                                         parts[i],
                                         partsize);
        }
        long double **res;
        pthread_t* threads=(pthread_t* )malloc(N*sizeof(pthread_t));

        res=(double long**)malloc(N*sizeof(long double*));
        for(int i=0;i<N;i++){
            pthread_create(&threads[i],NULL,
                monte_carlo_integrate_thread,Tdatas[i]);
        }
        estimate=0;
        for(int i=0;i<N;i++){
            pthread_join(threads[i],(void*)&res[i]);
            estimate+=*res[i];
        }
        estimate/=size;
        //printf("estimate: %Lf",estimate);

        // Your pthreads code ends here

        timer.c_end = clock();
        clock_gettime(CLOCK_MONOTONIC, &timer.t_end);
        gettimeofday(&timer.v_end, NULL);

        if(DEBUG && VERBOSE){
            print_array(results, n_threads);
        }
    }

    if(DEBUG){
        if(VERBOSE){
            print_array(samples, size);
            printf("Estimate: [%.33LF]\n", estimate);
        }
        printf("%.16LF, [%f, clock], [%f, clock_gettime], [%f, gettimeofday]\n",
               estimate,
               (double) (timer.c_end - timer.c_start) / (double) CLOCKS_PER_SEC,
               (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
               (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0,
               (double) (timer.v_end.tv_sec - timer.v_start.tv_sec) +
               (double) (timer.v_end.tv_usec - timer.v_start.tv_usec) / 1000000.0);
    } else {
        printf("%.16LF, %f\n",
               estimate,
               (double) (timer.t_end.tv_sec - timer.t_start.tv_sec) +
               (double) (timer.t_end.tv_nsec - timer.t_start.tv_nsec) / 1000000000.0);
    }//Membro1
    return 0;
}


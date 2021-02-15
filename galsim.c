#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "graphics/graphics.h"


#define EPSILON 10e-3
#define DELTA_T 10e-5
#define L 1
#define W 1

typedef struct particle
{
    double x;
    double y;
    double mass;
    double v_x;
    double v_y;
    double brightness;
} particle_t;


particle_t** read_particles(int N, FILE* gal_file);
void free_particles(int N, particle_t **particles);
double calculate_r_ij(particle_t *p_i, particle_t *p_j);
double* calculate_F_i(int N, particle_t **particles, double* F_i, double G, int is_x);
void write_particles(particle_t **particles, int N, int n_steps);
void step(int N, particle_t **p, double* F_i_x, double* F_i_y);

float circle_radius=0.0015, circle_color=0;
const int window_width=1000;

int main(int argc, char *argv[])
{
    if(argc != 6)
    {
        printf("Usage: %s <no of particles> <filename> <number of timesteps> <delta t> <graphics [0/1]>\n", argv[0]);
        exit(1);
    }

    // int L = 1, W = 1;
    int N = atoi(argv[1]);
    char* filename = argv[2];
    int n_steps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int is_graphics = atoi(argv[5]);
    double G = 100/N;
    // double epsilon = 1e-3;
    // delta_t = 10e-5;
    printf(" === Warning: delta_t, epsilon and G are set as constants.\n");
    printf("%d %s %d %f %d\n", N, filename, n_steps, delta_t, is_graphics);

    // Open the file for reading.
    FILE *gal_file = fopen(filename, "rb");

    if(!gal_file)
    {
        // Pointer is null, something went wrong!
        printf("File was not opened correctly.\n");
        exit(EXIT_FAILURE);
    }

    // Read the initial state of the particles.
    particle_t **p = read_particles(N, gal_file);
    // And close now that the particles are stored in memory.
    if(fclose(gal_file))
    {
        printf("\nFile was not closed correctly.\n");
        exit(EXIT_FAILURE);
    }

    double *F_i_x = calloc(N, sizeof(particle_t));
    double *F_i_y = calloc(N, sizeof(particle_t));
    if (is_graphics)
    {
        InitializeGraphics(argv[0], window_width, window_width);
        SetCAxes(0,1);
        int counter = 0;
        do
        {
            counter++;
            F_i_x = calculate_F_i(N, p, F_i_x, G, 1);
            F_i_y = calculate_F_i(N, p, F_i_y, G, 0);

            ClearScreen();
            // One step.
            step(N, p, F_i_x, F_i_y);
            for (size_t i = 0; i < N; i++)
            {
                DrawCircle(p[i]->x, p[i]->y, L, W, circle_radius, circle_color);
            }
            Refresh();
            usleep(30000);
            if(counter + 1 == n_steps) printf("Passed n_steps\n");

        } while (!CheckForQuit() && is_graphics);
        FlushDisplay();
        CloseDisplay();
    }
    if(!is_graphics)
    {
        // When graphics are turned off, no while-loop.
        for (size_t n = 0; n < n_steps; n++)
        {
            F_i_x = calculate_F_i(N, p, F_i_x, G, 1);
            F_i_y = calculate_F_i(N, p, F_i_y, G, 0);
            step(N, p, F_i_x, F_i_y);
        }
        write_particles(p, N, n_steps);
    /*
    printf("\n");
    printf("|||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
    printf("VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n");
        for (size_t j = 0; j < N; j++) {
            printf("\nPARTICLE %ld\n\tx: %0.2f\t\ty: %0.2f\n\tx_vel: %0.2f\ty_vel: %0.2f\n\tmass: %0.2f\n\tbrightness %0.2f\n", j+1, p[j]->x, p[j]->y,  p[j]->v_x, p[j]->v_y, p[j]->mass, p[j]->brightness);
        }
    }
    */
    // Put back when doing real sims?
    // particle_t **particles_check = read_particles(N, strcat("ref_output_data/", gal_file));
    }
    free(F_i_x);
    free(F_i_y);
    free_particles(N, p);
    //free_particles(N, particles_check);
    return 0;
}

void step(int N, particle_t **p, double* F_i_x, double* F_i_y)
{
    for (size_t i = 0; i < N; i++)
    {
        /*
        a_x = (F_i_x[i] / particles[i]->mass);
        a_y = (F_i_y[i] / particles[i]->mass);
        */
        /*
        for (size_t j = 0; j < N - 1; j++)
        {
            if(i != j)
            {
                F_x += (((-G* (p[i]->mass) * (p[j]->mass)) * (p[i]->x - p[j]->x) ) / (pow(sqrt((pow((p[i]->x) - (p[j]->x), 2) + pow((p[i]->y) - (p[j]->y), 2))) + EPSILON, 3)));
                F_y += (((-G* (p[i]->mass) * (p[j]->mass)) * (p[i]->y - p[j]->y) ) / (pow(sqrt((pow((p[i]->x) - (p[j]->x), 2) + pow((p[i]->y) - (p[j]->y), 2))) + EPSILON, 3))) ;
            }
        }
        */
        p[i]->v_x += DELTA_T*(F_i_x[i] / p[i]->mass);
        p[i]->v_y += DELTA_T*(F_i_y[i] / p[i]->mass);

        p[i]->x += DELTA_T*(p[i]->v_x);
        p[i]->y += DELTA_T*(p[i]->v_y);
    }
}

double* calculate_F_i(int N, particle_t **particles, double* F_i, double G, int is_x)
{
    particle_t *p_i, *p_j;
    double r_ij, denominator, Gm_i, r_xy;
    /*
    if(is_x)
        printf("X:\n");
    else
        printf("Y:\n");
    */
    for (size_t i = 0; i < N; i++)
    {
        // For each particle.
        p_i = particles[i];
        Gm_i = (-G*(p_i->mass));
        for (size_t j = 0; j < N; j++)
        {
            // Calculate the force from the others, not counting oneself.
            if(i != j)
            {
                p_j = particles[j];
                // calculate r_ij.
                r_ij = calculate_r_ij(p_i, p_j);
                // printf("\tr %ld -> %ld is %f\n", i + 1, j + 1, r_ij);
                denominator = pow(r_ij + EPSILON, 3);
                // Which component to compute is the input of is_x.
                if(is_x)
                {
                    r_xy = p_i->x - p_j->x;
                }
                else
                {
                    r_xy = p_i->y - p_j->y;
                }
                F_i[i] += (Gm_i * ((p_j->mass) / (denominator)))*r_xy;
                // printf("P%ld: %0.2f N [P%ld]\n", i + 1, F_i[i], j + 1);
            }
        }
    }
    return F_i;
}

double calculate_r_ij(particle_t *p_i, particle_t *p_j)
{
    return sqrt((pow((p_i->x) - (p_j->x), 2) + pow((p_i->y) - (p_j->y), 2)));
}

particle_t** read_particles(int N, FILE* gal_file)
{
    // Allocate an array of particles.
    particle_t **particles = malloc(N*sizeof(particle_t*));

    // while(!feof(gal_file))
    // Read N initial particles and their data.
    for (size_t i = 0; i < N; i++)
    {
        particle_t *p = (particle_t *)malloc(sizeof(particle_t));
        fread(p, sizeof(particle_t), 1, gal_file);
        if(!feof(gal_file))
        {
            // printf("\nPARTICLE %ld\n\tx: %0.2f\t\ty: %0.2f\n\tx_vel: %0.2f\ty_vel: %0.2f\n\tmass: %0.2f\n\tbrightness %0.2f\n", i+1, p->x, p->y,  p->v_x, p->v_y, p->mass, p->brightness);
            particles[i] = p;
        }
    }
    return particles;
}

void write_particles(particle_t **particles, int N, int n_steps)
{
    FILE * out_file;

    out_file = fopen("result.gal", "wb");
    for (size_t i = 0; i < N; i++)
    {
        fwrite(particles[i], sizeof(particle_t), 1, out_file);
    }
    int closed = fclose(out_file);
    if(closed) printf("didn't close write correctly\n");
}

void free_particles(int N, particle_t **particles)
{
    for (size_t i = 0; i < N; i++)
    {
        free(particles[i]);
    }
    free(particles);
}

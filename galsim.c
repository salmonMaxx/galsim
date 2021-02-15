#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "graphics/graphics.h"


#define EPSILON 0.01
#define DELTA_T 0.00001
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
    // Zero-allocate two arrays for the forces on each particle.
    // These will be updated for each time step of the simulation.
    double *F_i_x = calloc(N, sizeof(particle_t));
    double *F_i_y = calloc(N, sizeof(particle_t));

    // Branching off for the graphics, this will not save output to file.
    if (is_graphics)
    {
        InitializeGraphics(argv[0], window_width, window_width);
        SetCAxes(0,1);
        // int counter = 0;
        do
        {
            // counter++;
            // For each step in the while loop, update the forces on each particle.
            F_i_x = calculate_F_i(N, p, F_i_x, G, 1);
            F_i_y = calculate_F_i(N, p, F_i_y, G, 0);

            ClearScreen();
            // Calculate one time step with the updated forces.
            step(N, p, F_i_x, F_i_y);
            for (size_t i = 0; i < N; i++)
            {
                // Draw all the particles on the screen.
                DrawCircle(p[i]->x, p[i]->y, L, W, circle_radius, circle_color);
            }
            Refresh();
            usleep(3000);
        } while (!CheckForQuit() && is_graphics);
        FlushDisplay();
        CloseDisplay();
    }

    // Branching for when no graphics are needed and when output needs to be
    // saved to a file.
    if(!is_graphics)
    {
        // When graphics are turned off, no while-loop.
        for (size_t n = 0; n < n_steps; n++)
        {
            // Calculate forces for time step n.
            F_i_x = calculate_F_i(N, p, F_i_x, G, 1);
            F_i_y = calculate_F_i(N, p, F_i_y, G, 0);

            // Take a step.
            step(N, p, F_i_x, F_i_y);
        }
        // Save particles to file.
        write_particles(p, N, n_steps);
    }
    // Freeing all the used memory.
    free(F_i_x);
    free(F_i_y);
    free_particles(N, p);
    return 0;
}

void step(int N, particle_t **p, double* F_i_x, double* F_i_y)
{
    // Updates each particle component-wise with the forces given by F_i_x and F_i_y.
    for (size_t i = 0; i < N; i++)
    {
        // Updating velocities first, so that the position is then calculated
        // with v_{i+1}.

        // a_i = F_i / m
        // => v_{i+1} = v_i + dt * a_i
        p[i]->v_x += DELTA_T*(F_i_x[i] / p[i]->mass);
        p[i]->v_y += DELTA_T*(F_i_y[i] / p[i]->mass);
        // => x_{i+1} = x_i + dt * v_{i+1}
        p[i]->x += DELTA_T*(p[i]->v_x);
        p[i]->y += DELTA_T*(p[i]->v_y);
    }
}

double* calculate_F_i(int N, particle_t **particles, double* F_i, double G, int is_x)
{
    particle_t *p_i, *p_j;
    double r_ij, denominator, Gm_i, r_xy;

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
                F_i[i] = (Gm_i * ((p_j->mass) / (denominator)))*r_xy;
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

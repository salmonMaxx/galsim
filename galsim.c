#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "graphics/graphics.h"
#include "assert.h"

#define EPSILON 0.001
#define L 1
#define W 1

double delta_t;

typedef struct vector2
{
    double x;
    double y;
} vector2_t;

typedef struct particle
{
    vector2_t position;
    double mass;
    vector2_t velocity;
    double brightness;
} particle_t;

typedef struct save_format
{
    double x;
    double y;
    double mass;
    double v_x;
    double v_y;
    double brightness;
} save_format_t;

particle_t** read_particles(int N, FILE* gal_file);
void free_particles(int N, particle_t **particles);
vector2_t* calculate_F_i(int N, particle_t **p, vector2_t *F, double G);
void step(int N, particle_t **p, vector2_t *F);
void write_particles(particle_t **particles, int N, int n_steps);

float circle_radius=1.5*0.00128, circle_color=0;
const int window_width=1000;

int main(int argc, char *argv[])
{
    if(argc != 6)
    {
        printf("Usage: %s <no of particles> <filename> <number of timesteps> <delta t> <graphics [0/1]>\n", argv[0]);
        exit(1);
    }
    int N = atoi(argv[1]);
    char* filename = argv[2];
    int n_steps = atoi(argv[3]);
    delta_t = atof(argv[4]);
    int is_graphics = atoi(argv[5]);

    particle_t **p;
    vector2_t *F;
    double G = 100.0/N;

    // Open the file for reading.
    FILE *gal_file = fopen(filename, "rb");

    if(!gal_file)
    {
        // Pointer is null, something went wrong!
        printf("File was not opened correctly.\n");
        exit(EXIT_FAILURE);
    }
    // Read the initial state of the particles.
    p = read_particles(N, gal_file);

    // And close now that the particles are stored in memory.
    if(fclose(gal_file))
    {
        printf("\nFile was not closed correctly.\n");
        exit(EXIT_FAILURE);
    }

    F = calloc(N, sizeof(vector2_t));

    if (is_graphics)
    {
        InitializeGraphics(argv[0], window_width, window_width);
        SetCAxes(0,1);
        do
        {
            F = calculate_F_i(N, p, F, G);
            ClearScreen();
            step(N, p, F);
            for (size_t i = 0; i < N; i++)
            {
                DrawCircle(p[i]->position.x, p[i]->position.y, L, W, circle_radius, circle_color);
            }
            Refresh();
            usleep(3000);
        } while (!CheckForQuit() && is_graphics);
        FlushDisplay();
        CloseDisplay();
    }
    else
    {
        for (size_t n = 0; n < n_steps; n++)
        {
            F = calculate_F_i(N, p, F, G);
            step(N, p, F);
        }
        write_particles(p, N, n_steps);
    }
    free(F);
    free_particles(N, p);
    return 0;
}

void step(int N, particle_t **p, vector2_t *F)
{
    for (size_t i = 0; i < N; i++)
    {
        p[i]->velocity.x += delta_t*(F[i].x / p[i]->mass);
        p[i]->velocity.y += delta_t*(F[i].y / p[i]->mass);
        p[i]->position.x += delta_t*(p[i]->velocity.x);
        p[i]->position.y += delta_t*(p[i]->velocity.y);
    }
}

vector2_t* calculate_F_i(int N, particle_t **p, vector2_t *F, double G)
{
    double r_ij, denominator, Gm_i, temp;

    for (size_t i = 0; i < N; i++)
    {
        F[i].x = 0;
        F[i].y = 0;
        Gm_i = (-G*(p[i]->mass));
        for (size_t j = 0; j < N; j++)
        {
            if (i != j)
            {
                r_ij = sqrt((((p[i]->position.x) - (p[j]->position.x))*((p[i]->position.x) - (p[j]->position.x))) + (((p[i]->position.y) - (p[j]->position.y))*((p[i]->position.y) - (p[j]->position.y))));
                denominator = (r_ij + EPSILON)*(r_ij + EPSILON)*(r_ij + EPSILON);
                temp = (Gm_i * ((p[j]->mass) / (denominator)));
                F[i].x += temp*(p[i]->position.x - p[j]->position.x);
                F[i].y += temp*(p[i]->position.y - p[j]->position.y);
            }
        }
    }
    return F;
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

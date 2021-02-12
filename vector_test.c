#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct vector2
{
    double x;
    double y;
} vector2_t;

typedef struct particle
{
    vector2_t position;
    vector2_t velocity;
    double mass;
    double brightness;
} particle_t;

int main(int argc, char const *argv[]) {
    particle_t *p;
    p = malloc(sizeof(particle_t));
    return 0;
}

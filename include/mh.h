#ifndef _MH
#define _MH

#define MAX_ARRAY_INT 1000  // Tamaño máximo suficiente para los problemas concretos

    typedef struct {
        int array_int[MAX_ARRAY_INT];  // Cambio a array de longitud estática
        double fitness;
    } Individuo;

    void cruzar(Individuo *, Individuo *, Individuo *, Individuo *, int, int);
    void mutar(Individuo *, int, int, double);
    void fitness(const double *, Individuo *, int, int);
    double distancia_ij(const double *, int, int, int);
#endif
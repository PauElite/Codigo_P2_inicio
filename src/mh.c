#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#include "../include/mh.h"

#define MUTATION_RATE 0.15
#define PRINT 1

int aleatorio(int n) {
    return (rand() % n);  // genera un numero aleatorio entre 0 y n-1
}

int find_element(int *array, int end, int element)
{
    int i = 0;
    int found = 0;

    // comprueba que un elemento no está incluido en el individuo (el cual no admite enteros repetidos)
    while((i < end) && !found) {
        if(array[i] == element) {
            found = 1;
        }
        i++;
    }
    return found;
}

void crear_individuo(int n, int m, Individuo *individuo)
{
    int i = 0, value;

    // inicializa array de elementos
    memset(individuo->array_int, -1, m * sizeof(int));

    while(i < m) {
        value = aleatorio(n);
        // si el nuevo elemento no está en el array...
        if(!find_element(individuo->array_int, i, value)) {
            individuo->array_int[i] = value;  // lo incluimos
            i++;
        }
    }
}

int comp_array_int(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

int comp_fitness(const void *a, const void *b) {
    /* qsort pasa un puntero al elemento que está ordenando */
    return (*(Individuo **)b)->fitness - (*(Individuo **)a)->fitness;
}

double aplicar_mh(const double *d, int n, int m, int n_gen, int tam_pob, double m_rate, int *sol)
{
    // para reiniciar la secuencia pseudoaleatoria en cada ejecución
    srand(time(NULL) + getpid());
    int i, g, mutation_start;

    // crea poblacion inicial (array de individuos)
    Individuo **poblacion = (Individuo **) malloc(tam_pob * sizeof(Individuo *));
    assert(poblacion);

    // crea cada individuo (array de enteros aleatorios)
    for(i = 0; i < tam_pob; i++) {
        poblacion[i] = (Individuo *) malloc(sizeof(Individuo));
        crear_individuo(n, m, poblacion[i]);

        // calcula el fitness del individuo
        fitness(d, poblacion[i], n, m);
    }

    // evoluciona la poblacion durante un numero de generaciones
    for(g = 0; g < n_gen; g++)
    {    

        // los hijos de los ascendientes mas aptos sustituyen a la ultima mitad de los individuos menos aptos
        for(i = 0; i < (tam_pob / 2) - 1; i += 2) {
            cruzar(poblacion[i], poblacion[i + 1], poblacion[tam_pob / 2 + i], poblacion[tam_pob / 2 + i + 1], n, m);
        }

        // inicia la mutacion a partir de 1/4 de la poblacion
        mutation_start = tam_pob / 4; // Indice por donde empieza a realizar la mutación

        // muta 3/4 partes de la poblacion
        for(i = mutation_start; i < tam_pob; i++) {
            mutar(poblacion[i], n, m, m_rate);
        }

        // recalcula el fitness del individuo
        for(i = 0; i < tam_pob; i++) {
            fitness(d, poblacion[i], n, m);
        }

        // ordena individuos segun la funcion de bondad (mayor "fitness" --> mas aptos)
        qsort(poblacion, tam_pob, sizeof(Individuo *), comp_fitness);

        if (PRINT) {
            printf("Generacion %d - ", g);
            printf("Fitness = %.0lf\n", poblacion[0]->fitness);
        }
    }
    // ordena el array solucion
    qsort(poblacion[0]->array_int, m, sizeof(int), comp_array_int);

    // y lo mueve a sol para escribirlo
    memmove(sol, poblacion[0]->array_int, m * sizeof(int));

    // almacena el mejor valor obtenido para el fitness
    double value = poblacion[0]->fitness;

    // se libera la memoria reservada
    for(i = 0; i < tam_pob; i++) {
        free(poblacion[i]);  // Liberamos cada individuo
    }
    free(poblacion);  // Finalmente, liberamos el array de punteros a individuos

    // devuelve el valor obtenido para el fitness
    return value;
}

void cruzar(Individuo *padre1, Individuo *padre2, Individuo *hijo1, Individuo *hijo2, int n, int m)
{
    // Elegir un "punto" de corte aleatorio a partir del que se realiza el intercambio de los genes
    int punto_corte = aleatorio(m);  // Elegimos un punto de corte aleatorio
    int i, j;

    // Los primeros genes del padre1 van al hijo1. Idem para el padre2 e hijo2.
    for(i = 0; i < punto_corte; i++) {
        hijo1->array_int[i] = padre1->array_int[i];
        hijo2->array_int[i] = padre2->array_int[i];
    }

    // Y los restantes son del otro padre, respectivamente.
    for(i = punto_corte, j = punto_corte; i < m && j < m; j++) {
        // Asignamos genes del padre2 al hijo1 si no están ya presentes
        hijo1->array_int[i] = padre2->array_int[j];
        hijo2->array_int[i] = padre1->array_int[j];
        i++;
    }

    // Comprobamos que los genes añadidos después del punto de corte, no estuviesen ya presentes. Si se da el caso, modificarlos por otros aleatorios que no estén
    for(int i = punto_corte; i < m; i++) {
        if(find_element(hijo1->array_int, i, hijo1->array_int[i])) {
            int flag = 0;
            int num_aleatorio;
            while(!flag) {
                num_aleatorio = aleatorio(n);
                if(!find_element(hijo1->array_int, i, num_aleatorio)) {
                    hijo1->array_int[i] = num_aleatorio;
                    flag = 1;
                }
            }
        }

        if(find_element(hijo2->array_int, i, hijo2->array_int[i])) {
            int flag = 0;
            int num_aleatorio;
            while(!flag) {
                num_aleatorio = aleatorio(n);
                if(!find_element(hijo2->array_int, i, num_aleatorio)) {
                    hijo2->array_int[i] = num_aleatorio;
                    flag = 1;
                }
            }
        }
    }
}

void mutar(Individuo *actual, int n, int m, double m_rate)
{
    // Decidir cuantos elementos mutar:
    // Si el valor es demasiado pequeño la convergencia es muy pequeña y si es demasiado alto diverge
    int num_aleatorio;
    int num_mutaciones = m * m_rate;
    for (int i = 0; i < num_mutaciones; i++) {

        // Encontrar un valor aleatorio que no esté presente en el individuo
        num_aleatorio = aleatorio(n);
        while (find_element(actual->array_int, m, num_aleatorio)) {
            num_aleatorio = aleatorio(n);
        }

        // Asignamos el nuevo valor a la posición seleccionada
        actual->array_int[aleatorio(m)] = num_aleatorio;
    }
}

double distancia_ij(const double *d, int i, int j, int n) {
    // Asumiendo que d representa una matriz de distancias en forma compacta (como se menciona en el código)
    if (i == j) {
        return 0.0;
    } else if (i > j) {
        int aux = j;
        j = i;
        i = aux;
    }

    int k = (((n * n) - n) / 2) - ((((n - i) * (n - i)) - (n - i)) / 2) + j - i - 1;
    return d[k];
}


void fitness(const double *d, Individuo *individuo, int n, int m)
{
    // Determina la calidad del individuo calculando la suma de la distancia entre cada par de enteros
    double suma = 0.0;
    for(int i = 0; i < m; i++) {
        for(int j = i + 1; j < m; j++) {
            suma = suma + distancia_ij(d, individuo->array_int[i], individuo->array_int[j], n);
        }
    }
    individuo->fitness = suma;
}

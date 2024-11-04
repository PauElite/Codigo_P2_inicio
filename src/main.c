#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>
#include <stddef.h>
#include <string.h>
#include "../include/mh.h"

#include "../include/io.h"
#define TIME 1
#define NGM 20

extern void aplicar_mh(const double *d, int n, int m, int g, int tam_pob, double m_rate, Individuo *poblacion, int rank);
int aleatorio(int n);
int find_element(int *array, int end, int element);
void crear_individuo(int n, int m, Individuo *individuo);
int comp_fitness(const void *a, const void *b);
int comp_fitness_menorAMayor(const void *a, const void *b);
void crear_tipo_datos(int m, MPI_Datatype *individuo_type)
{
    int blocklen[2] = {m, 1};
    MPI_Datatype dtype[2] = {MPI_INT, MPI_DOUBLE};

    MPI_Aint disp[2];
    disp[0] = offsetof(Individuo, array_int);
    disp[1] = offsetof(Individuo, fitness);

    MPI_Type_create_struct(2, blocklen, disp, dtype, individuo_type);
    MPI_Type_commit(individuo_type);
}

void realizarMigracion(MPI_Datatype individuo_type, Individuo *poblacion, int tam_pob, int rank, int size, int NEM)
{
    // Cada subpoblación se encuentra ya ordenada por el valor fitness
    Individuo *mejoresIndividuos = (Individuo *)malloc(NEM * sizeof(Individuo));
    for (int i = 0; i < NEM; i++)
    {
        mejoresIndividuos[i] = poblacion[i];
    }

    Individuo *mejoresIndividuos_total = NULL;
    if (rank == 0)
    {
        mejoresIndividuos_total = (Individuo *)malloc(size * NEM * sizeof(Individuo));
    }

    MPI_Gather(mejoresIndividuos, NEM, individuo_type, mejoresIndividuos_total, NEM, individuo_type, 0, MPI_COMM_WORLD);

    free(mejoresIndividuos);
    if (rank == 0)
    {
        qsort(mejoresIndividuos_total, size * NEM, sizeof(Individuo), comp_fitness);
    }
    qsort(poblacion, tam_pob, sizeof(Individuo), comp_fitness_menorAMayor); // ordeno de menor a mayor para sustituir los peores individuos

    MPI_Scatter(mejoresIndividuos_total, NEM, individuo_type, poblacion, NEM, individuo_type, 0, MPI_COMM_WORLD);
    qsort(poblacion, tam_pob, sizeof(Individuo), comp_fitness);
}

static double mseconds()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec * 1000 + t.tv_usec / 1000;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv); // Inicializar entorno MPI

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 0;
    int m = 0;
    int n_gen = 0;
    int tam_pob = 0;
    double m_rate = 0.0;

    // Check Number of Input Args
    if (rank == 0)
    {
        if (argc < 5)
        {
            fprintf(stderr, "Ayuda:\n");
            fprintf(stderr, "  mpirun -np <n_procesos> ./programa n m nGen tamPob mRate\n");
            MPI_Finalize();
            return (EXIT_FAILURE);
        }

        n = atoi(argv[1]);
        m = atoi(argv[2]);
        n_gen = atoi(argv[3]);
        tam_pob = atoi(argv[4]);
        m_rate = atof(argv[5]);

        // Check that 'm' is less than 'n'
        assert(m < n);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tam_pob, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m_rate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *d = NULL;
    if (rank == 0)
    {
        d = read_distances(n);
    }
    else
    {
        d = (double *)malloc(((n * n - n) / 2) * sizeof(double));
    }

    MPI_Bcast(d, (n * n - n) / 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Crear tipo de datos MPI para Individuo
    MPI_Datatype individuo_type;
    crear_tipo_datos(m, &individuo_type);

    // Inicializar población local
    int tam_pob_local = tam_pob / size;
    int resto = tam_pob % size;
    if (rank < resto)
    {
        tam_pob_local++; // Distribuir el resto de manera equitativa entre los primeros procesos
    }
    Individuo *poblacion_total = NULL;
    if (rank == 0)
    {
        poblacion_total = (Individuo *)malloc(tam_pob * sizeof(Individuo));
        for (int i = 0; i < tam_pob; i++)
        {
            crear_individuo(n, m, &poblacion_total[i]);
        }
    }

    Individuo *poblacion = (Individuo *)malloc(tam_pob_local * sizeof(Individuo));
    int *elementosADifundir = (int *)malloc(size * sizeof(int));
    int *posicionesIniciales = (int *)malloc(size * sizeof(int));
    int desplazamiento = 0;
    for (int i = 0; i < size; i++)
    {
        elementosADifundir[i] = (tam_pob / size) + (i < resto ? 1 : 0);
        posicionesIniciales[i] = desplazamiento;
        desplazamiento += elementosADifundir[i];
    }
    int NEM = tam_pob / size / 2;

    MPI_Scatterv(poblacion_total, elementosADifundir, posicionesIniciales, individuo_type, poblacion, tam_pob_local, individuo_type, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate memory for output data
    int *sol = (int *)malloc(m * sizeof(int));

#ifdef TIME
    MPI_Barrier(MPI_COMM_WORLD);
    double ti = mseconds();
#endif
    int ngm = 0;
    for (int g = 0; g < n_gen; g++)
    {
        // Cada proceso evoluciona su subpoblación
        aplicar_mh(d, n, m, g, elementosADifundir[rank], m_rate, poblacion, rank);
        ngm++;
        if (ngm == NGM)
        {   
            // Al llegar a NGM evoluciones, se realiza la migración de los mejores NEM individuos y se resetea el contador NGM
            realizarMigracion(individuo_type, poblacion, elementosADifundir[rank], rank, size, NEM);
            ngm = 0;
        }
    }

    // Nos quedamos con el mejor valor fitness de entre los calculados por cada proceso
    double global_value = 0.0;
    MPI_Reduce(&poblacion->fitness, &global_value, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
#ifdef TIME
        double tf = mseconds();
        printf("Execution Time: %.2lf sec\n", (tf - ti) / 1000);
#endif
        printf("Best Fitness Value: %.2lf\n", global_value);
    }

    // Free Allocated Memory
    if (rank == 0)
    {
        free(poblacion_total);
    }
    free(poblacion);
    free(sol);
    free(d);
    free(elementosADifundir);
    free(posicionesIniciales);

    MPI_Type_free(&individuo_type);
    MPI_Finalize(); // Finalizar entorno MPI

    return (EXIT_SUCCESS);
}

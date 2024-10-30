#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>
#include <stddef.h>
#include "../include/mh.h"

#include "../include/io.h"
#define TIME 1

extern double aplicar_mh(const double *, int, int, int, int, double, int *);
void crear_tipo_datos(int m, MPI_Datatype *individuo_type) {
    int blocklen[2] = {m, 1};
    MPI_Datatype dtype[2] = { MPI_INT, MPI_DOUBLE };
    
    MPI_Aint disp[2];
    disp[0] = offsetof(Individuo, array_int);
    disp[1] = offsetof(Individuo, fitness);
    
    MPI_Type_create_struct(2, blocklen, disp, dtype, individuo_type); 
    MPI_Type_commit(individuo_type);
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

        // Generate matrix D with distance values among elements
        // double *d = read_distances(n);
    }
    // printf("\nSoy el proceso %d y estos son mis datos: n=%d\tm=%d\tn_gen=%d\ttam_pob=%d\tm_rate=%.2lf", rank, n, m, n_gen, tam_pob, m_rate);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tam_pob, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m_rate, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // printf("\nSoy el proceso %d y estos son mis datos: n=%d\tm=%d\tn_gen=%d\ttam_pob=%d\tm_rate=%.2lf\n", rank, n, m, n_gen, tam_pob, m_rate);


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

    MPI_Scatterv(poblacion_total, elementosADifundir, posicionesIniciales, individuo_type, poblacion, tam_pob_local, individuo_type, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < tam_pob_local; i++) {
        printf("Proceso %d: %d\t%d\t%d\t%d con fitness %.2lf\n", rank, poblacion[i].array_int[0], poblacion[i].array_int[2], poblacion[i].array_int[3], poblacion[i].array_int[4], poblacion[i].fitness);
    }

    if (rank == 0)
    {
        free(poblacion_total);
    }
    free(elementosADifundir);
    free(posicionesIniciales);

    // Allocate memory for output data
    int *sol = (int *)malloc(m * sizeof(int));

#ifdef TIME
    MPI_Barrier(MPI_COMM_WORLD);
    double ti = mseconds();
#endif

    // Call Metaheuristic
    //double value = aplicar_mh(d, n, m, n_gen, tam_pob / size, m_rate, sol);

    // Gather results at root process
    double global_value;
    //MPI_Reduce(&value, &global_value, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
#ifdef TIME
        double tf = mseconds();
        printf("Execution Time: %.2lf sec\n", (tf - ti) / 1000);
#endif
        printf("Best Fitness Value: %.2lf\n", global_value);
    }

    // Free Allocated Memory
    free(sol);
    free(d);
    free(poblacion);

    MPI_Type_free(&individuo_type);
    MPI_Finalize(); // Finalizar entorno MPI

    return (EXIT_SUCCESS);
}

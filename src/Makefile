# Definicion del compilador y las banderas
CC = mpicc
CFLAGS = -O3 -Wall -std=c99 -g

# Definicion de los parámetros de entrada
N = 500
M = 200
N_GEN = 100
T_POB = 50
M_RATE = 0.15

# Nombre del ejecutable
EXEC = sec
OUTFILE = ../output/out
C_FILES = main.c io.c mh.c

# Regla para compilar el ejecutable con MPI
sec: $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC) -DTIME -DDEBUG -lm

# Regla para ejecutar el programa con MPI, usando el archivo de entrada y parámetros adicionales
run_sec: sec
	mpirun -np 1 ./$(EXEC) < ../input/input_500_200_100_50.txt $(N) $(M) $(N_GEN) $(T_POB) $(M_RATE)

# Regla para limpiar el ejecutable generado
clean:
	rm -f $(EXEC)

.PHONY: all clean run_sec

CC = g++
CFLAGS = -O3 -Wall -g

N = 50
M = 5
N_GEN = 4
T_POB = 8
M_RATE = 0.15

EXEC = igen
OUTFILE = ../input
C_FILES = igen.cpp

igen: $(C_FILES)
	$(CC) $(CFLAGS) $(C_FILES) -o $(EXEC)

exe_igen:
	./$(EXEC) $(N) $(M) $(N_GEN) $(T_POB) $(M_RATE) > $(OUTFILE)_$(N)_$(M)_$(N_GEN)_$(T_POB).txt

clean:
	rm -f $(EXEC)

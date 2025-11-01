CC = g++
CFLAGS = -c -O3 -I/usr/include/eigen3 -fopenmp -lgmp -lmpfr -std=c++17 -lgsl -lgslcblas -lm
LDFLAGS = -fopenmp -lgmp -lmpfr -std=c++17 -lgsl -lgslcblas -lm
SRC = WSA.cpp Trasser.cpp PFSS_Data.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = hello

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXEC)

.PHONY: all clean
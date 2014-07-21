CFLAGS = -Wall -O3 -g -pthread #-pg

CC = g++ $(CFLAGS)

RM = rm -f

GP_OBJS = main.o io.o graph.o graph_compression.o getopt_pp.o perf.o

all: graph_compression

graph_compression: $(GP_OBJS)
	@echo ""
	@echo " --- graph_compression ---"
	@$(CC) $(GP_OBJS) -o graph_compression
	@echo ""

%.o: %.cc %.h
	@echo " --- COMPILING OBJECT \"$@\""
	@$(CC) $< -c

clean:
	$(RM) graph_compression *.o
	clear


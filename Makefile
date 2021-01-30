all: mst

mst: mst-skeleton.c mst-solution.c mst-kruskal.c
	mpicc -O3 mst-skeleton.c -lm -o mst

clean:
	rm mst

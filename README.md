# How to run
* Open the project in terminal (or `cd` to the project folder).
* Type `make all`.
* The `.cpp` and `.h` files will be compiled and the new file called `mst` will appear. This is the executable file.
* Generate data using `create-graph.py`: `python create-graph.py 1000 4000 5000 graph_test`. The words `1000`, `4000`, `5000`, `graph_test` here are the number of vertexes, the number of edges, the maximum edge weight and the test file name, respectively. These information can be adjusted. The generated data should be placed in a separate folder called `data`.
* Run the `mst` file by the following command: `mpirun --hostfile /etc/hosts --oversubscribe -np 8 ./mst data/medium_test kruskal_par`.


GraphProject4990_2023F
graph generator can be used to generate file for graph of any size.\
g++ GraphGenerator.cpp -o ggen -O3\

./ggen 0 0 10000 100 0 499 OutputGraph.txt\
This command generates a graph of size |V| = 1 000 000, |E| = 49 519 998\
argument 1 is the number of complete graphs which are subgraphs to the output graph\
argument 2 is size of those complete graphs. for example 2 5 will mean that 2 K_5 graphs will be subgraphs to the output graph\
argument 3 is the number of clique graphs which are subgraphs to the output\
argument 4 is the size of those cliques. for example 10 000 100 will generate 10 000 cliques each with size 100 vertices\
argument 5 is whether or not to allow "3 edges", that is an edge with 3 vertices instead of 2. I don't understand this so I always leave it 0\
argument 6 is a random seed, I use 499\
argument 7 is the name of the file to output the graph into.

Compile C-Style dfs:\
g++ dfs.cpp -o dfs -O3 -std=c++20

run:\
./dfs 32V.txt //any graph can be substituted for 32V.txt, but it must be in the same folder as the compiled program

Compile C++ style dfs:\
g++ dfs2-using-2d-vector.cpp -o dfs2 -O3 -std=c++20

run:\
./dfs2 32V.txt  //Program does not yet produce an output file.

C style is slightly more consistent in having fast speeds, but both are almost equally fast.

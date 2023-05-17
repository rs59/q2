# Project "q2"

## Algorithm summary

A p-way partitioning of graph G is a division of G into p sub-graphs in which the vertices in each subset do not overlap and satisfy specific properties:
- The sum of the weights of the nodes in each subgraph is balanced.
- The sum of the weights of the edges crossing between subsets is minimized.

## Project outputs

1. The source C++ files, including: 
     - the sequential version of the algorithm
     - the parallel version of the algorithm
2. A readme USER MANUAL text file (written in plain ASCII) describing how to:
    - compile and run the program,
    - under which system,
    - which API.
3. A DOCUMENTATION text file (written in markdown) including:
    - How the reading part has been performed
    - How the data structure has been organized
    - How the parallelism has been designed
    - Experimental evaluation: tables or graphics reporting a reasonable set of experimental results.
    - Experimental evaluation: comparison (in terms of memory and of elapsed time) between the original sequential version of the tool and the 2-3 parallel versions with different parallelization levels (1, 2, 4, 8, threads) and increasing complexity of the input graph (up to millions of nodes).
    - Optional: Compare the result of the written tool with publicly available applications such as hMetis, or the one available in Python in scikit-learn.
4. An OVERHEAD set (organized in PowerPoint or similar) to be used during the project discussion in the project evaluation phase.


To read the full project details, view the [instructions](https://github.com/rs59/q2/blob/main/instructions.pdf) PDF. Additional relevant research can be found in the [resources](https://github.com/rs59/q2/tree/main/resources) folder.

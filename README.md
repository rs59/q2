# Project "q2"

## Algorithm summary

A p-way partitioning of graph G is a division of G into p sub-graphs in which the vertices in each
subset do not overlap and satisfy specific properties:
- The sum of the weights of the nodes in each subgraph is balanced.
- The sum of the weights of the edges crossing between subsets is minimized.

## Project outputs

- The source C++ files.
    - Including the sequential and the parallel version of the algorithm.
- A README text file (written in plain ASCII) including a short “user manual”, i.e., a document describing how to compile and run the program, under which system, which API, etc.
- A DOCUMENTATION text file (written in Word, Latex, Mark-down, etc.) including a short “designer manual”, i.e., a document including:
-  All main design choices (how the reading part has been performed, how the data structure has been organized, how the parallelism has been designed, etc.)
-  The experimental evaluation of the tools, i.e., tables or graphics reporting a reasonable set of experimental results. The experimental evidence should include a comparison (in terms of memory and of elapsed time) between the original sequential version of the tool and the 2-3 parallel versions with different parallelization levels (i.e., with 1, 2, 4, 8, etc., threads) and different size of the input graph (up to millions of nodes).
-  An OVERHEAD set (organized in PowerPoint or similar) to be used during the project discussion in the project evaluation phase.

The comparison must consider the computation time and the memory usage of the main phases (at least, graph loading and path computation) and compare the sequential version with the parallel one with an increasing number of threads (i.e., 1, 2, 3, 4, etc.) on graphs of increasing size and complexity.

Optional: Compare the result of the written tool with publicly available applications such as hMetis, or the one available in Python in scikit-learn.



To read the full project details, view the [instructions](https://github.com/rs59/q2/blob/main/instructions.pdf) PDF. Additional relevant research can be found in the [resources](https://github.com/rs59/q2/tree/main/resources) folder.

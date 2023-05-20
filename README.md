# Project "q2"

## Algorithm summary

A p-way partitioning of graph G is a division of G into p sub-graphs in which the vertices in each subset do not overlap and satisfy specific properties:
- The sum of the weights of the nodes in each subgraph is balanced.
- The sum of the weights of the edges crossing between subsets is minimized.

## Project outputs

1. The [source C++ files](https://github.com/rs59/q2/blob/main/src/2023partition.cpp), including: 
     - the sequential version of the algorithm
     - the parallel version of the algorithm
2. A readme [USER MANUAL](https://github.com/rs59/q2/blob/main/output/USER_MANUAL.md) text file (written in plain ASCII) describing how to:
    - compile and run the program, under which system, and which API.
3. A DOCUMENTATION text file (written in Word/Markdown/LaTeX/etc) including:
  
[Part 1](https://github.com/rs59/q2/blob/main/output/DOCUMENTATION_1.md):
 - How the reading part has been performed
 - How the data structure has been organized
 - How the parallelism has been designed
 - Selection of algorithm: determine what heuristics make sense (i.e. level of imbalance between the partitions and level of failure of edges)
  
[Part 2](https://github.com/rs59/q2/blob/main/output/DOCUMENTATION_2.ipynb):
 - Experimental evaluation: comparison (in terms of [memory and of elapsed time](https://unix.stackexchange.com/questions/207209/how-to-calculate-the-memory-consumed-by-a-c-program-in-linux)) between the original sequential version of the tool and the parallel version with different parallelization levels (1, 2, 4, 8, threads) (Y axis) and increasing complexity of the input graph (up to millions of nodes) (X axis).
 - Tables or graphics reporting a reasonable set of experimental results.
  
[Part 3](https://github.com/rs59/q2/blob/main/output/DOCUMENTATION_3.ipynb): 
 - Evaluate a standard graph sample dataset that can be used with the tool (real-world dataset)
  
[Part 4](https://github.com/rs59/q2/blob/main/output/DOCUMENTATION_4.ipynb):
 - Optional: Compare the result of the written tool with publicly available applications such as hMetis, or the one available in Python in scikit-learn.

4. An OVERHEAD set (organized in PowerPoint or similar) to be used during the project discussion in the project evaluation phase.


To read the full project details, view the [instructions](https://github.com/rs59/q2/blob/main/instructions.pdf) PDF. Additional relevant research can be found in the [resources](https://github.com/rs59/q2/tree/main/resources) folder.

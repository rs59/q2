---
USER MANUAL: 2023partition
---
A readme USER MANUAL text file (written in plain ASCII) describing how to run the program

### Running the program:

`2023partition nthreads npartitions maxdeviation inputfile outputfile`
  ex. `2023partition 2 3 1.05 input.gv output.gv`

- `nthreads`: number of threads, 1 or greater (1 if sequential, more if parallelizable)
- `npartitions`: number of partitions desired in the output graph
- `maxdeviation`: maximum deviation of the output partition size, expressed as a float value greater than 1.00. This value is multplied by the ideal partition size (sum_weights/npartitions) to create an upper bound for the partition size
- `inputfile`: input file for graph input (see format info below)
- `outputfile`: output file for graph output (see format info below)

### Program output:

`[301,305,296]`

If the program run is successful, the size (sum of vertex weights) of each partition.

`Unable to split graph: ERROR`

If the program run is unsuccessful, a message indicating failure and the error reason.

### Input file

The program reads in source files with inputs in the format `V(n1,4) V(n2,3) ....  E(n1,n2,6)` where:
- All non-newline whitespace is removed
- Any line starting with a # is disregarded
- Vertex information: The format V(X,Y) is matched, with X being any string of characters before the comma, and Y being any value that can be can be cast to a positive int
- Edge information: The format E(X,Y,Z) is matched, with X being any string of characters before the comma, Y being any string of characters before the comma, and Z being any value that can be can be cast to a positive int
- Any other information is disregarded

Example input file:

    #define _STR(x) #x
    #define STR(x) _STR(x)
    
    #define V(A,B) A [label=B];
    #define E(A,B,C) A -- B [label="  "+STR(C)];
    
    graph G {
      
        V(n1,4)
        V(n2,3)
        V(n3,2)
        V(n4,1)
    
        E(n1,n3,6)
        E(n3,n2,7)
        E(n1,n4,9)
        E(n2,n4,4)
    }

Simple rendering (node weight in brackets, edge weight between dashes):

    [1] --4-- [3] --7-- [2] --6-- [4]
     \_____________9______________/

### Output file

The program outputs source files with in the format `x1 subgraph c1 { V(n1,4) V(n4,1) } subgraph c2 { V(n2,3) V(n3,2) } ....  E(n1,n2,6) x2` where:
- x1 consists of any text before the first instance of the term `subgraph`
- Partition information: marked by `subgraph X {` and closed by `}`, with X indicating any unique string of characters not shared by other subgraphs. Between the curly braces, the number of matches of `V(` indicates the number of vertices in the partition.
- Vertex information: The format `V(X,Y)` is matched, with X being any string of characters before the comma, and Y being any value that can be can be cast to a positive int. Vertices are only inside partitions.
- Edge information: The format `E(X,Y,Z)` is matched, with X being any string of characters before the comma, Y being any string of characters before the comma, and Z being any value that can be can be cast to a positive int. Edges are not inside partitions.
- Any other information (x2) is disregarded.


Example output file:
    
    #define _STR(x) #x
    #define STR(x) _STR(x)
    
    #define V(A,B) A [label=B];
    #define E(A,B,C) A -- B [label="  "+STR(C)];
    
    graph G {
      
        subgraph cluster1 {
          V(n1,4)
          V(n4,1)
        }
    
        subgraph cluster2 {
          V(n2,3)
          V(n3,2)
        }
    
        E(n1,n3,6)
        E(n3,n2,7)
        E(n1,n4,9)
        E(n2,n4,4)
    }

Simple rendering (node weight in brackets, edge weight between dashes, split marked using `/` character):

         [3] --7-- [2]
         |           |
    //// 4 ///////// 6 ////
         |           |
         [1] --9-- [4]
 
### Compile and run the program

The program uses the standard C++ libraries (no Boost) and can be compiled by running on the latest version of Ubuntu (22.04):
- `mkdir bin`
- `cd bin && g++ -std=c++17 ../src/2023partition.cpp -o 2023partition`
- `cd bin && chmod +x 2023partition`
- `cd bin && ./2023partition`

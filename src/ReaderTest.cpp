#include "Reader.cpp"
#include "graph.h"

#define FILENAME "resources/metismodels/x15y30m20q20.metis"
#define NUMTHREADS 4

int main(){

    Graph graph = metisRead(FILENAME, NUMTHREADS);
    graph.print();
}
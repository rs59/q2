#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include "graph.h"

std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


void vertexRead(std::ifstream& file, Graph& graph, unsigned int nodeId, std::mutex& mtx){
    int nodeWeight;
    std::string line;
    while (std::getline(file, line)) {
            std::stringstream sLine(line);
            if(sLine >> nodeWeight){
                mtx.lock();
                graph.addVertex(nodeId, static_cast<double>(nodeWeight) );
                mtx.unlock();
            } else {
                std::cerr << "Invalid File Data Format" << std::endl;
            }
        nodeId++;
        }
}

void edgeRead(std::ifstream& file, Graph& graph, unsigned int nodeId, std::mutex& mtx){
    int nodeEdge, edgeWeight;
    std::string line;
    while (std::getline(file, line)) {
            std::stringstream sLine(line);
            if(sLine >> nodeEdge){ // first value already read -> dump
                while(sLine >> nodeEdge){
                    sLine >> edgeWeight;
                    mtx.lock();
                    graph.addEdge(nodeId, nodeEdge, static_cast<double>(edgeWeight));
                    mtx.unlock();
                }
            } else {
                std::cerr << "Invalid File Data Format" << std::endl;
            }
        nodeId++;
        }
}


void readFromFile(std::ifstream& file, Graph& graph, int threadId, int numThreads, unsigned int startingLine, std::mutex& mtx) {
    if (!file.is_open()) {
        GotoLine(file, startingLine);
        vertexRead(file, graph, startingLine, mtx);
        //barrier to implement
        GotoLine(file, startingLine);
        edgeRead(file, graph, startingLine, mtx);
    } else {
        std::cerr << "Error: Could not open file " << std::endl;
    }
}



int main() {
    const std::string filename = "your_file.txt";
    const int numThreads = 4; // You can change the number of threadPool

    Graph graph;
    std::vector<std::thread> threadPool;

    for (int i = 0; i < numThreads; ++i) {
        threadPool.emplace_back(readFromFile, filename, std::ref(graph), i, numThreads);
    }

    for (auto& thread : threadPool) {
        thread.join();
    }

    return 0;
}
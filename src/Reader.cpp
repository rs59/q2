#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "graph.h"

int current_count;

std::ifstream& GotoLine(std::ifstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


void vertexRead(std::ifstream& file, Graph& graph, const unsigned int& nodeId, const unsigned int& linesToRead,  std::mutex& mtx){
    int nodeWeight;
    std::string line;
    for(int i = 0; std::getline(file, line) && i < linesToRead; i++) {
            std::stringstream sLine(line);
            if(sLine >> nodeWeight){
                std::unique_lock<std::mutex> lock(mtx);
                graph.addVertex(nodeId+i, static_cast<double>(nodeWeight));
                lock.unlock();
            } else {
                std::cerr << "Invalid File Data Format" << std::endl;
            }
        }
}

void edgeRead(std::ifstream& file, Graph& graph, const unsigned int& nodeId, const unsigned int& linesToRead, std::mutex& mtx){
    int nodeEdge, edgeWeight;
    std::string line;
    for(int i = 0; std::getline(file, line) && i < linesToRead; i++) {
            std::stringstream sLine(line);
            if(sLine >> nodeEdge){ // first value already read -> dump
                while(sLine >> nodeEdge){
                    sLine >> edgeWeight;
                    std::unique_lock<std::mutex> lock(mtx);
                    graph.addEdge(nodeId + i, nodeEdge, static_cast<double>(edgeWeight));
                    lock.unlock();
                }
            } else {
                std::cerr << "Invalid File Data Format" << std::endl;
            }
        }
}


void readFromFile(const std::string& filename, Graph& graph, const int& numThreads, const unsigned int& startNodeId, const unsigned int& nodeToRead,  std::mutex& mtx, std::mutex& barMtx, std::condition_variable& cv, std::mutex& creationMtx) {

    creationMtx.unlock();
    std::ifstream file(filename);
    if (file.is_open()) {
        
        GotoLine(file, startNodeId);
        vertexRead(file, graph, startNodeId, nodeToRead, mtx);

        // barrier implementation with cv
        std::unique_lock<std::mutex> lock(barMtx);
        current_count++;
        if (current_count == numThreads)
        {
                // All threads have arrived, notify all waiting threads
                current_count = 0;
                cv.notify_all();
        }
        else
        {
                // Not all threads have arrived, wait
                cv.wait(lock, []
                        { return current_count == 0; });
        }

        GotoLine(file, startNodeId);
        edgeRead(file, graph, startNodeId, nodeToRead, mtx);
        file.close();
    } else {
        std::cerr << "Error: Could not open file " << std::endl;
    }
}

Graph metisRead(const std::string& filename, const int& numThreads){
    
    if(numThreads < 1) {
        std::cerr << "Null or negative number of threads" << std::endl;
        exit(1);
    }

    Graph graph;
    std::string line;
    unsigned int numLines;

    std::vector<std::thread> threadPool;

    std::mutex mtx, creationMtx;

    std::condition_variable cv;
    std::mutex barMtx;
    current_count = 0;
    
    std::ifstream file(filename);
    if(!file.is_open()){
        std::cerr << "File not Valid" << std::endl;
        exit(1);
    }
    std::getline(file, line);
    file.close();
    std::stringstream sLine(line);
    sLine >> numLines;
    int linesPerThread = (numThreads>1) ? numLines/(numThreads-1) : numLines;
    //std::cout << "Number of Lines: " << numLines << std::endl;
    //std::cout << "Number of Lines per Thread: " << linesPerThread << std::endl;
    for (int i = 0; i < numThreads; ++i) {
        creationMtx.lock();
        int firstNodeId = linesPerThread*i + 1;
        //std::cout << "First line for thread: "<< i+1 << ": " << firstNodeId << std::endl;
        threadPool.emplace_back([&] {
            readFromFile(filename, graph, numThreads, firstNodeId, linesPerThread, std::ref(mtx), std::ref(barMtx), std::ref(cv), std::ref(creationMtx));
            });
    }

    for (auto& thread : threadPool) {
        thread.join();
    }

    return graph;
}


using namespace std;

void printGraphVertexes(Graph& graph){
    auto vertices = graph.getVertices();
    for (auto& vertex : vertices){
        cout << "Id: " << vertex.first << "\tWeight: " << vertex.second << endl;
    }
}

void printGraphEdges(Graph& graph){
    auto vertices = graph.getEdgeWeights();
    for (auto& edge : vertices){
        cout << "Node1: " << edge.first.first << "\tNode2: "<< edge.first.second << "\tWeight: " << edge.second << endl;
    }
}

int main() {
    //resources\metismodels\x200000y440000m20q20.metis
    const string filename = "resources/metismodels/x200000y440000m20q20.metis";
    const int numThreads = 7; // You can change the number of threadPool
    Graph graph = metisRead(filename, numThreads);
    /*printGraphVertexes(graph);
    printGraphEdges(graph);*/
    return 0;
}
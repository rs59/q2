#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "graph.h"


int current_count;

std::ifstream& gotoLine(std::ifstream& file, unsigned int num){
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
                //std::cout << nodeId + i << ": " << nodeWeight << std::endl;
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


void readFromFile(const std::string& filename, Graph& graph, const int& i, const int& numThreads, const int startNodeId, const int nodeToRead,  std::mutex& mtx, std::mutex& barMtx, std::condition_variable& cv, std::mutex& creation_mtx) {

    creation_mtx.unlock();
    std::ifstream file(filename); // create new reference to file
    if (file.is_open()) {
        
        DEBUG_STDOUT("Unlocked creation_mtx, file open");
        DEBUG_STDOUT("Will read i "+std::to_string(i)+" startnode "+std::to_string(startNodeId)+" numbertoread "+std::to_string(nodeToRead));
        gotoLine(file, startNodeId);
        DEBUG_STDOUT("Thread now reading "+std::to_string(startNodeId)+" "+std::to_string(nodeToRead));
        vertexRead(file, graph, startNodeId, nodeToRead, mtx);
        DEBUG_STDOUT("Completed "+std::to_string(startNodeId)+" "+std::to_string(nodeToRead));

        // barrier implementation with cv
        std::unique_lock<std::mutex> lock(barMtx);
        current_count++;
        if (current_count == numThreads)
        {
                // All threads have arrived, notify all waiting threads
                DEBUG_STDOUT("All threads arrived.");
                current_count = 0;
                cv.notify_all();
        }
        else
        {
                // Not all threads have arrived, wait
                DEBUG_STDOUT("Thread now waiting....");
                cv.wait(lock, []
                        { return current_count == 0; });
                DEBUG_STDOUT("Thread active....");
        }

        gotoLine(file, startNodeId);
        edgeRead(file, graph, startNodeId, nodeToRead, mtx);
        file.close();
    } else {
        std::cerr << "Error: Could not open file " << std::endl;
    }
}

Graph metisRead(const std::string& filename, const int& numThreads){
    
    DEBUG_STDOUT("numThreads: "+std::to_string(numThreads));
    if(numThreads < 1) {
        std::cerr << "Null or negative number of threads" << std::endl;
        exit(1);
    }

    Graph graph;
    std::string line;
    unsigned int numLines;

    std::vector<std::thread> threadPool;

    std::mutex mtx, creation_mtx;

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
    int linesPerThread = numLines/numThreads;
    int lineStartsAt = 1-(numLines/numThreads);

    DEBUG_STDOUT("Number of Lines: "+std::to_string(numLines));


    for (int i = 0; i < numThreads; ++i) {
        creation_mtx.lock();

        lineStartsAt+=(numLines/numThreads);
       
        
        DEBUG_STDOUT("Number of Lines assigned to thread #"+std::to_string(i)+": "+std::to_string(linesPerThread));
        DEBUG_STDOUT("thread #"+std::to_string(i)+" lines ["+std::to_string((numLines/numThreads)*i + 1)+"--"+std::to_string(((numLines/numThreads)*i + 1) + linesPerThread-1)+"]");


        // Assigns last thread all remaining lines
        if(i == numThreads - 1) {
          threadPool.emplace_back([&] {                                  // all remaining lines
              readFromFile(filename, graph, i, numThreads, lineStartsAt, numLines-(linesPerThread*(numThreads-1)), std::ref(mtx), std::ref(barMtx), std::ref(cv), std::ref(creation_mtx));
              });
        } else {
          threadPool.emplace_back([&] {
              readFromFile(filename, graph, i, numThreads, lineStartsAt, linesPerThread, std::ref(mtx), std::ref(barMtx), std::ref(cv), std::ref(creation_mtx));
              });
        }
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

/*
int main() {
    //resources\metismodels\x200000y440000m20q20.metis
    const string filename = "resources/metismodels/x15y30m20q20.metis";
    const int numThreads = 7; // You can change the number of threadPool
    Graph graph = metisRead(filename, numThreads);
    //graph.print();
    return 0;
}*/

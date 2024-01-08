#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include "graph.h"


#define TRUE 1
#define FALSE 0

struct edge{
    int neighbourNode;
    double weight;

    edge(const int _node, const double _weight) : neighbourNode(_node), weight(_weight) {}
};

struct node{
    int id;
    double weight;
    std::vector<edge> neighbours;

    node(const int _id, const double _weight, const std::vector<edge> _neighbours) : id(_id),  weight(_weight), neighbours(_neighbours) {}
};


struct mysync {
    std::mutex creationMtx, barrierMtx, writingMtx, sectionMtx;
    std::shared_mutex eofMtx;
    std::condition_variable barrier_cv;
};

int section;
int current_count;
char endOfFile = FALSE;



int getEOF(std::shared_mutex& mtx){
    std::shared_lock lock(mtx);
    return endOfFile;
}

void setEOF(std::shared_mutex& mtx, char value){
    std::shared_lock lock(mtx);
    endOfFile = value;
}

std::ifstream& gotoLine(std::ifstream& file, unsigned int num){
    file.seekg(0, std::ios::beg);
    file.clear();
    for(unsigned int i=0; i < num; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        if(file.eof()){
            return file;
        }
    }
    return file;
}


void readLines(std::ifstream& file, std::vector<node>& nodes, const int& startId, const unsigned int& linesToRead){
    double nodeWeight, edgeWeight;
    int nodeNeighbour;
    std::string line;
    for(unsigned int i = 0; std::getline(file, line) && i < linesToRead; i++) {
        std::stringstream sLine(line);
        std::vector<edge> neighbours;
        if(sLine >> nodeWeight){
            while(sLine >> nodeNeighbour){
                if(sLine >> edgeWeight){
                    neighbours.emplace_back(nodeNeighbour, edgeWeight);
                } else {
                    std::cerr << "Invalid File Data Format" << std::endl;
                }
            }
        } else {
            std::cerr << "Invalid File Data Format" << std::endl;
        }
            nodes.emplace_back(startId+i, nodeWeight, neighbours);
    }
    if(file.eof()){
        DEBUG_STDOUT("END OF FILE REACHED -----------------------------------------------------------------");
    }
}

void addVerticesToGraph(Graph& graph, const std::vector<node>& nodes, std::mutex& writeMtx)
{
    std::unique_lock<std::mutex> lock(writeMtx);
    for(const auto& node: nodes){
        graph.addVertex(node.id, node.weight);
    }
}

void addEdgesToGraph(Graph& graph, const std::vector<node>& nodes, std::mutex& writeMtx){
    std::unique_lock<std::mutex> lock(writeMtx);
    for(const auto& node: nodes){
        for(const auto& nNode: node.neighbours){
            graph.addEdge(node.id, nNode.neighbourNode, nNode.weight);
        }
    }
}

int sectionUpdate(std::mutex& mtx) {
    std::unique_lock lock(mtx);
    int temp = section;
    section++;
    return temp;
}



void readFromFile(int threadNum, const std::string& filename, Graph& graph, const int& numThreads, const int nodeToRead,  mysync& sincro) {

    sincro.creationMtx.unlock();
    
    std::ifstream file(filename); // create new reference to file
    std::vector<node> nodes;
    if (file.is_open()) {
        while(!getEOF(sincro.eofMtx)){
            DEBUG_STDOUT("Unlocked creation_mtx, file open");
            int thisSection = sectionUpdate(sincro.sectionMtx);
            int startId = thisSection * nodeToRead + 1;
            DEBUG_STDOUT("Thread #" + std::to_string(threadNum) + " Section:" + std::to_string(thisSection));
            DEBUG_STDOUT("Thread #" + std::to_string(threadNum) + " Lines: ["+std::to_string(startId) + " - " + std::to_string(startId + nodeToRead) + "]");
            gotoLine(file, startId);
            readLines(file, nodes, startId, nodeToRead);
            DEBUG_STDOUT("Completed "+std::to_string(startId)+" - "+std::to_string(startId+nodeToRead));
            DEBUG_STDOUT(file.rdstate());
            std::string temp;
            if(!getline(file, temp)){
                DEBUG_STDOUT(file.rdstate());
                DEBUG_STDOUT("Called END FOR ALL");
                setEOF(sincro.eofMtx, TRUE);
            }
            DEBUG_STDOUT("Thread #" + std::to_string(threadNum) + " EOF: " +std::to_string(getEOF(sincro.eofMtx)));
        }
    } else {
        std::cerr << "Error: Could not open file " << std::endl;
    }
    file.close();
    addVerticesToGraph(graph, nodes, sincro.writingMtx);

    // barrier implementation with cv
    std::unique_lock<std::mutex> lock(sincro.barrierMtx);
    current_count++;
    if (current_count == numThreads)
    {
            // All threads have arrived, notify all waiting threads
            DEBUG_STDOUT("All threads arrived.");
            current_count = 0;
            sincro.barrier_cv.notify_all();
    }
    else
    {
            // Not all threads have arrived, wait
            DEBUG_STDOUT("Thread now waiting....");
            sincro.barrier_cv.wait(lock, []
                    { return current_count == 0; });
            DEBUG_STDOUT("Thread active....");
    }
        addEdgesToGraph(graph, nodes, sincro.writingMtx);
}

unsigned int getTotLines(const std::string& fileName){
    std::string line;
    unsigned int numLines;
    std::ifstream file(fileName);
    if(!file.is_open()){
        std::cerr << "File not Valid" << std::endl;
        exit(1);
    }
    std::getline(file, line);
    file.close();
    std::stringstream sLine(line);
    sLine >> numLines;
    return numLines;
}

Graph metisRead(const std::string& filename, const int& numThreads){
    
    DEBUG_STDOUT("numThreads: "+std::to_string(numThreads));
    if(numThreads < 1) {
        std::cerr << "Null or negative number of threads" << std::endl;
        exit(1);
    }

    Graph graph;

    std::vector<std::thread> threadPool;

    mysync sincro;
    current_count = 0;
    section = 0;
    setEOF(sincro.eofMtx, FALSE);
    
    unsigned int numLines = getTotLines(filename);
    int linesPerSection = numLines/(numThreads*1.5);
    
    DEBUG_STDOUT("Number of Lines: "+std::to_string(linesPerSection));

    for (int i = 0; i < numThreads; i++) {
        sincro.creationMtx.lock();
        
        DEBUG_STDOUT("Thread #"+std::to_string(i)+" created");

        threadPool.emplace_back([&] {                                  // all remaining lines
                readFromFile(i, filename, graph, numThreads, linesPerSection, std::ref(sincro));
            });
    }

    for (auto& thread : threadPool) {
        thread.join();
    }

    return graph;
}


void writeToFile(const std::vector<std::vector<int>>& partitions, const std::string& filename) {
    // Open the file
    std::ofstream of(filename);

    // Check if the file is opened successfully
    if (!of.is_open()) {
        std::cerr << "Error in Opening the file in output"<< std::endl;
        exit(-1);
    }

    // Write partitions to the file
    for (const auto& partition : partitions) {
        for (const auto& node : partition) {
            of << node << ' ';
        }
        of << std::endl; 
    }

    // Close the file
    of.close();
}


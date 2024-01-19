#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <condition_variable>
#include "graph.h"

//Boolean compiler definitions
#define TRUE 1
#define FALSE 0

// Strct to save a single Neighbout with the weight of the edge between the 2 nodes
struct edge{
    int neighbourNode;
    double weight;

    edge(const int _node, const double _weight) : neighbourNode(_node), weight(_weight) {}
};

// Struct to save the node and its neighbours
struct node{
    int id;
    double weight;
    std::vector<edge> neighbours;

    node(const int _id, const double _weight, const std::vector<edge> _neighbours) : id(_id),  weight(_weight), neighbours(_neighbours) {}
};

// Struct to compress the sincronizartion variables
struct mysync {
    std::mutex creationMtx, barrierMtx, writingMtx, sectionMtx;
    std::shared_mutex eofMtx;
    std::condition_variable barrier_cv;
};

// Global variables for sincronization
int section;
int current_count;
char endOfFile;

// Thread safe for reading and writing eof variable
int getEOF(std::shared_mutex& mtx){
    std::shared_lock lock(mtx);
    return endOfFile;
}

void setEOF(std::shared_mutex& mtx, char value){
    std::shared_lock lock(mtx);
    endOfFile = value;
}


// Function to move the file pointer to a specific line
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

// Function to insert all edges read into the graph data class
void addToGraph(Graph& graph, const std::vector<node>& nodes, std::mutex& writeMtx){
    std::unique_lock<std::mutex> lock(writeMtx);
    for(const auto& node: nodes){
        graph.setEdgeWeight(node.id, node.weight);
        for(const auto& nNode: node.neighbours){
            graph.addEdge(node.id, nNode.neighbourNode, nNode.weight);
        }
    }
}

// Function to read the file and save in the thread local variable the value read
void readLines(std::ifstream& file, Graph& graph, const int& startId, const unsigned int& linesToRead, std::mutex& writeMtx){
    std::vector<node> nodes;
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
    addToGraph(graph, nodes, writeMtx);
}

// Function to insert the value read to the Graph data-class
void addVerticesToGraph(Graph& graph, const std::vector<node>& nodes, std::mutex& writeMtx){
    std::unique_lock<std::mutex> lock(writeMtx);
    for(const auto& node: nodes){
        graph.addVertex(node.id, node.weight);
    }
}

// Function to insert all edges read into the graph data class
void addEdgesToGraph(Graph& graph, const std::vector<node>& nodes, std::mutex& writeMtx){
    std::unique_lock<std::mutex> lock(writeMtx);
    for(const auto& node: nodes){
        for(const auto& nNode: node.neighbours){
            graph.addEdge(node.id, nNode.neighbourNode, nNode.weight);
        }
    }
}

// Thread safe read and update the section variable
int sectionUpdate(std::mutex& mtx) {
    std::unique_lock lock(mtx);
    int temp = section;
    section++;
    return temp;
}


// Thread Function to read multiple sequential sections of lines
void readFromFile(int threadNum, const std::string& filename, Graph& graph, const int& numThreads, const int nodeToRead,  mysync& sincro) {
    // Unlock the mutex to proceed the creation of the next one
    sincro.creationMtx.unlock();
    // Create new reference to file
    std::ifstream file(filename); 
    // Check file successfully open
    if (file.is_open()) {
        // Keep reading until any thread hits the EOF
        while(!getEOF(sincro.eofMtx)){
            DEBUG_STDOUT("Unlocked creation_mtx, file open");
            // Update of the section of the file to Read
            int thisSection = sectionUpdate(sincro.sectionMtx);
            // Update of the start node id to give
            int startId = thisSection * nodeToRead + 1;
            DEBUG_STDOUT("Thread #" + std::to_string(threadNum) + " Section:" + std::to_string(thisSection));
            DEBUG_STDOUT("Thread #" + std::to_string(threadNum) + " Lines: ["+std::to_string(startId) + " - " + std::to_string(startId + nodeToRead) + "]");
            // Update File Pointer
            gotoLine(file, startId);
            // Read Data from file and save into @nodes
            readLines(file, graph, startId, nodeToRead, sincro.writingMtx);
            DEBUG_STDOUT("Completed "+std::to_string(startId)+" - "+std::to_string(startId+nodeToRead));
            // Check if the file reacher the end, file.eof() not working: flags not risen
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
    // End of Reading, close the File
    file.close();
}
// Function that reads the first line of the file 
// containing the total number of nodes and edges in the file
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
    // Data
    Graph graph;

    std::vector<std::thread> threadPool;

    // Syncronization Data
    mysync sincro;
    current_count = 0;
    section = 0;
    setEOF(sincro.eofMtx, FALSE);
    // Calculate the lines per thread to read
    unsigned int numLines = getTotLines(filename);
    DEBUG_STDOUT("Number of Lines: "+std::to_string(numLines));
    for(unsigned int i = 1; i <= numLines; i++){
        graph.addVertex(i, 1);
        //DEBUG_STDOUT("Added vertex :" + std::to_string(i));
    }
    DEBUG_STDOUT("Added nodes");
    int linesPerSection = numLines/(numThreads*1.5);
    
    
    // Creations of the threads 
    for (int i = 0; i < numThreads; i++) {
        sincro.creationMtx.lock();
        
        DEBUG_STDOUT("Thread #"+std::to_string(i)+" created");

        threadPool.emplace_back([&] {                                  // all remaining lines
                readFromFile(i, filename, graph, numThreads, linesPerSection, std::ref(sincro));
            });
    }
    // Waiting for all threads to terminate
    for (auto& thread : threadPool) {
        thread.join();
    }

    return graph;
}

// Function to write graph partitions to a file
void writeToFile(const std::vector<std::vector<int>>& partitions, const std::string& filename) {
    // Open the file
    std::ofstream of(filename);

    // Check if the file is opened successfully
    if (!of.is_open()) {
        std::cerr << "Error in Opening the file in output"<< std::endl;
        exit(-1);
    }
    std::cout << "WRITING TO FILE" << std::endl;

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


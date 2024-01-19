#define DEBUG 1

#ifdef DEBUG
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif

#include "ReaderWriter2.cpp"

// Function to calculate the cut size of partitions
double calculateCutSize(const Graph& G, const std::vector<std::vector<int>>& partitions){
    double cutSize = 0;
    for(const auto& partition : partitions)
        for(const auto& node : partition){
            auto nNodes = G.getNeighborsKL(node);
            for(const auto& nNode : nNodes)
                if(std::find(partition.begin(), partition.end(), nNode) == partition.end())
                    cutSize += G.getEdgeWeightKL(node, nNode);
        }
    return cutSize / 2;  // Dividing by 2 as each edge is counted twice
}
/*
std::vector<int> splitLine(std::string line){
    int temp;
    std::stringstream sLine(line);
    std::vector<int> partition;
    while(sLine >> temp){
        partition.push_back(temp);
    }
    return partition;
}*/

std::vector<std::vector<int>> readEndFile(std::string filename, int&  numNodes){
    // Create new reference to file
    std::ifstream file(filename);
    std::vector<std::vector<int>> partitions;
    if (file.is_open()) {
        if(file >> numNodes){
            int numPartitions;
            file >> numPartitions;
            for(int i = 0; i< numPartitions; i++){
                partitions.push_back(std::vector<int>());
            }
            int partitionIndex;
            int i = 1;
            while(file >> partitionIndex){
                partitions[partitionIndex].push_back(i);
                i++;
            }
        }else{
        DEBUG_STDOUT("Number of nodes not read");
        std::cerr << "Error: Could not read the File " << std::endl;
        exit(-1);
        }
    }
    else{
        std::cerr << "Error: Could not open file " << std::endl;
        exit(-1);
    }
    return partitions;
}


int main(int argc, char* argv[]){

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " input-filename output-filename numThreads" << std::endl;
        return 1;
    }

    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];
    int numThreads = std::stoi(argv[3]);

    if(numThreads < 1){
        std::cerr << "Usage: " << argv[0] << " input-filename output-filename numThreads" << std::endl;
        return 1;
    }
    DEBUG_STDOUT("Start Reading Input");
    Graph graph = metisRead(inputFilename, numThreads);
    DEBUG_STDOUT("Finished Reading");
    int numNodes;
    DEBUG_STDOUT("Started Reading output");
    auto partitions = readEndFile(outputFilename, numNodes);
    DEBUG_STDOUT("Finished Reading output");
    if(numNodes != graph.numVertices()){
        DEBUG_STDOUT("Inconguent Number of nodes ");
        std::cerr << "Inconguent Number of nodes " << std::endl;
        exit(-1);
    }
    // // Display the generated vectors
    // for (const auto& part : partitions) {
    //     std::cout << "Vector: ";
    //     for (const int& node : part) {
    //         std::cout << node << " ";
    //     }
    //     std::cout << std::endl;
    // }
    DEBUG_STDOUT("Started Calculating Cutsize");
    double cutSize = calculateCutSize(graph, partitions);
    DEBUG_STDOUT("Finished Calculating Cutsize");
    std::cout << "Cutsize: " << cutSize << std::endl;
}
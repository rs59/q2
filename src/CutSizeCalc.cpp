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

std::vector<int> splitLine(std::string line){
    int temp;
    std::stringstream sLine(line);
    std::vector<int> partition;
    while(sLine >> temp){
        partition.push_back(temp);
    }
    return partition;
}

std::vector<std::vector<int>> readEndFile(std::string filename){
    // Create new reference to file
    std::ifstream file(filename);
    std::vector<std::vector<int>> partitions;
    if (file.is_open()) {
        std::string line;
        while(std::getline(file, line)){
            partitions.push_back(splitLine(line));
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
    
    Graph graph = metisRead(inputFilename, numThreads);
    auto partitions = readEndFile(outputFilename);
    double cutSize = calculateCutSize(graph, partitions);
    std::cout << "Cutsize: " << cutSize << std::endl;
}
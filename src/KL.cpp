// #define DEBUG 1

#ifdef DEBUG 
#define DEBUG_STDERR(x) (std::cerr << (x) << std::endl)
#define DEBUG_STDOUT(x) (std::cout << (x) << std::endl)
#else 
#define DEBUG_STDERR(x)
#define DEBUG_STDOUT(x)
#endif

#include "multipartitionKL.cpp"
#include <iostream>
#include "RdWr.cpp"
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " nthreads npartitions inputfile outputfile mode: -rr for round robin or -b  for blob" << std::endl;
        return 1;
    }

    // Start the clock
    auto start_time = std::chrono::high_resolution_clock::now();

    int nthreads = std::stoi(argv[1]);
    int npartitions = std::stoi(argv[2]);
    std::string inputfile = argv[3];
    std::string outputfile = argv[4];
    std::string mode = argv[5];

    if(mode != "-rr" && mode != "-b"){
        std::cerr << "invalid mode or mode not selected" << std::endl;
        exit(1);
    }

    // initialize the Graph from reading the file
    Graph graph = metisRead(inputfile, nthreads);
    auto reading_end_time = std::chrono::high_resolution_clock::now();
    auto reading_duration = std::chrono::duration_cast<std::chrono::microseconds>(reading_end_time - start_time);
    double reading_seconds = reading_duration.count() / 1e6;
    std::cout << "Reading Execution time: " << reading_seconds << " seconds" << std::endl;

    // graph.print();

    std::vector<std::vector<int>> partitions;
    // Call the algorithm function based on the mode selected
    if(mode == "-rr") partitions = multipartitionKL_round_robin(graph, npartitions);
    else  partitions = multipartitionKL_blob(graph, npartitions);

#ifdef DEBUG
    printPartitions(partitions, graph.getVertices());
#endif


    auto partitioning_end_time = std::chrono::high_resolution_clock::now();
    auto partitioning_duration = std::chrono::duration_cast<std::chrono::microseconds>(partitioning_end_time - reading_end_time);\
    double partitioning_seconds = partitioning_duration.count() / 1e6;
    std::cout << "Partition Execution time: " << partitioning_seconds << " seconds" << std::endl;
    
    // Write the result on file
    writeToFile(partitions, outputfile);
    auto writing_end_time = std::chrono::high_resolution_clock::now();
    auto writing_duration = std::chrono::duration_cast<std::chrono::microseconds>(writing_end_time - partitioning_end_time);
    double writing_seconds = writing_duration.count() / 1e6;
    std::cout << "Writing Execution time: " << writing_seconds << " seconds" << std::endl;
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(writing_end_time - start_time);
    double seconds = duration.count() / 1e6;
    // Print the execution time
    std::cout << "Overall Execution time: " << seconds << " seconds" << std::endl;

    return 0;
}


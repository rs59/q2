#include <iostream>
#include <string>

using namespace std;


// Define the usage instructions for the executable.
string USAGE_INSTRUCTIONS = "\n" "\n" "Usage:" "\n"
"2023partition nthreads npartitions maxdeviation inputfile outputfile" "\n"
"ex. 2023partition 2 3 1.05 input.gv output.gv" "\n"
 "\n"
"nthreads: number of threads, 1 or greater (1 if sequential, more if parallelizable)" "\n"
"npartitions: number of partitions desired in the output graph" "\n"
"maxdeviation: maximum deviation of the output partition size, expressed as a float value greater than 1.00." "\n"
"This value is multplied by the ideal partition size (sum_weights/npartitions) to create an upper bound for the partition size" "\n"
"inputfile: input file for graph input (see format info below)" "\n"
"outputfile: output file for graph output (see format info below)";


// Utility function to test if a string is a positive int.
bool isPositiveInt(string str) {
    return !str.empty() && str.find_first_not_of("0123456789") == std::string::npos && std::stoi(str) > 0;
}

// Custom exception to throw up to terminal output.
class CustomException : public std::exception {
    private:
    string message;

    public:
    CustomException(string msg) : message(msg) {}
    string what () {
        return message;
    }
};


int main(int argc, char** argv)
{
    // Check arguments
    try {
        if(argc!=6) { // There are 6 arguments (1 program name + 5 parameters)
            throw CustomException("Wrong number of arguments provided.");
        } else if(!isPositiveInt(argv[1]) || !isPositiveInt(argv[2])) {
            // The first two arguments are positive ints
            throw CustomException("The first two arguments (nthreads and npartitions) must be positive, nonzero integers."+USAGE_INSTRUCTIONS);
        }



    } catch (CustomException ce) {
        cout << ce.what() << endl;
        exit(1);
    }
}
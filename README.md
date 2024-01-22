------------------------------------------------------------
                     Q2 project README
------------------------------------------------------------

**Compatibility:**
The program can be compiled on Ubuntu 22.04.3 LTS with g++ 11.4.0. Alternatively, it can be compiled on Windows 10 with g++ 13.2.0 (MSYS2 project). The implementation relies solely on standard libraries, without using Boost.

**Compile METIS Implementation:**
Compile:
    g++ -Wall -o ./bin/mtmetis.exe ./src/mtmetis.cpp
Run:
    ./bin/mtmetis.exe nThreads nPartitions maxDeviation inputFile outputFile cutoffSize partitioningAlg
    - nThreads: number of threads for partitioning
    - nPartitions: the final number of partitions desired
    - maxDeviation: tolerance for partition weights
    - inputFile: address of the input file
    - outputFile: address to save the result
    - cutoffSize: size below which partitioning can occur
    - partitioningAlg: -g for greedy or -kl for KL

**Compile and Run KL Implementation:**
Compile:
    g++ -Wall -o ./bin/KL.exe ./src/KL.cpp
Run:
    ./bin/KL.exe nThreads nPartitions input_filename output_filename mode
    - nThreads: number of threads (affecting only reading in KL sequential mode)
    - nPartitions: the final number of partitions desired
    - input_filename: address of the input file
    - output_filename: address to save the result
    - mode: -rr for round robin or -b for blob

**METIS Models and Testing:**
./resources/metismodels contains various Metis-type model files, generated using genmetis.py (provided in ./resources). Each file structure is xaaaaybbbbmccqdd.metis:
- xaaaaybbbb: number of nodes and edges
- mcc: maximum value of node weights
- qdd: maximum value of edge weights

**CutSize and Cut Quality Calculation:**
Use CutSizeCalc.cpp to calculate cutsize and cut quality.
Compile:
    g++ -Wall -o ./bin/CutSizeCalc.exe ./src/CutSizeCalc.cpp
Run:
    ./bin/CutSizeCalc.exe input_filename output_filename nThreads
    - input_filename: address of the input file
    - output_filename: address to save the result
    - nThreads: number of threads for the file reader

**Installing and Running MT-METIS:**
Run the following commands in bash:
curl https://dlasalle.github.io/mt-metis/releases/mt-metis-0.7.2.tar.gz | tar -xz
cd mt-metis-0.7.2; ./configure; make; make install
apt-get install sysstat
sudo apt install zsh -y
g++ ./q2/src/mtmetis.cpp -o ./q2/bin/mtmetis.exe; chmod +x ./q2/bin/mtmetis.exe
./mt-metis-0.7.2/build/Linux-x86_64/bin/mtmetis -T 4 ./q2/resources/metismodels/x10000y20000m20q20.metis 10 test.part

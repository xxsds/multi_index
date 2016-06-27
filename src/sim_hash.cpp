#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "sim_hash/sim_hash.hpp"

/* Compute SimHash and OddSketch for each row read from stdin. Hashes are stored in two binary files, 8 bytes each */

int main(int argc, char* argv[]) {

    if( argc < 2 ) {
       std::cout << "Usage: ./" << argv[0] << " output_hash_filename [window_size]" << std::endl;
       std::cout << "\twindow_size: Size of the window sliding over the text. Default window_size=2" << std::endl;
       return 1;
    }
    size_t window_size = 2;
    if ( argc > 2 ) {
        window_size = std::stoull(argv[2]);    
    }

    std::string outfilename1 = std::string(argv[1]) + std::string(".SimHash");
    std::string outfilename2 = std::string(argv[1]) + std::string(".OddSketch");
    std::ofstream file1 (outfilename1.c_str());
    std::ofstream file2 (outfilename2.c_str());
    
    std::string line;
    uint64_t simhash, oddsketch;
    uint64_t tot = 0; 
    while (std::getline(std::cin, line)) {
         simhash = simhash64(line, window_size);
         oddsketch = oddsketch64(line, window_size);
         file1.write(reinterpret_cast<const char *>(&simhash), sizeof(simhash));
         file2.write(reinterpret_cast<const char *>(&oddsketch), sizeof(oddsketch));
         tot++;
         if(tot % 10000000 == 0) std::cout << "Processed " << tot/1000000 << " Million rows." << std::endl;
    }
    return 0;
}

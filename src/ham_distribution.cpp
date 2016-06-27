#include <iostream>
#include <fstream>
#include <array>
#include <cstdint>
#include <sdsl/bits.hpp>

using namespace std;

int main(int argc, char* argv[]){
    if ( argc < 2 ){
        cout << "Usage: ./ham_distribution file " << endl;
        return 1;
    }
    ifstream file(argv[1], std::ifstream::binary);
    if ( file.is_open()){
        array<uint64_t, 65> cnt = {0};
        uint64_t x;
        while ( file.read((char*)&x, sizeof(x) ) ){
            ++cnt[sdsl::bits::cnt(x)];
        }
        for(size_t i=0; i<cnt.size(); ++i){
            cout << i << "," << cnt[i] << endl;
        } 
    } else {
        cout << "Unable to open file" << endl;
        return 1;
    }

}

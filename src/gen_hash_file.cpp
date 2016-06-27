#include <iostream>
#include <sdsl/int_vector.hpp>

using namespace std;

int main(int argc, char* argv[]){
    if ( argc < 2 ) {
        cout << "Usage: ./" << argv[0] << "file" << endl;
        return 1;
    }
    vector<uint64_t> keys = {0x00000000000003FF,
                             0xFF800000000003FF};
    ofstream out(argv[1], std::ofstream::binary);
    for (size_t i=0; i < keys.size(); ++i){
        out.write((char*)&keys[i], sizeof(keys[i]));
    }
}

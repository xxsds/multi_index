#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

using namespace std;

int main(){
    vector<uint64_t> data = {0x1111222233334444,
                             0x1011222233334444,
                             0x1111222233334454,
                             };
    ofstream out("debug.SimHash");
    out.write((char*)data.data(), sizeof(data[0])*data.size());
}

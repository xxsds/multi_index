#include "chi/chi_idx.hpp"
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace phi;
using namespace sdsl;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

const string index_name = INDEX_NAME;


vector<uint64_t> unique_vec(vector<uint64_t> v){
    if(v.size()>0){
        std::sort(v.begin(), v.end());
        auto end = std::unique(v.begin(), v.end());
        v.erase(end, v.end()); 
    }
    return v;
}

vector<uint64_t> load_keys(string hash_file)
{
    vector<uint64_t> keys;
    int_vector_buffer<64> key_buf(hash_file, ios::in, 1<<20, 64, true);
//    cout << "# hashes: " << key_buf.size() << endl;
    for (size_t i=0; i<key_buf.size(); ++i){
        keys.push_back(key_buf[i]);
    }
    return unique_vec(keys);
}


int main(int argc, char* argv[]){
    constexpr uint8_t t_n = N;
    constexpr uint8_t t_k = K;
    typedef INDEX_TYPE      index_type;
    if ( argc < 1 ) {
        cout << "Usage: ./" << argv[0] << " hash_file" << endl;
        return 1;
    }
    string hash_file = argv[1];
    string idx_file = hash_file + "."+to_string(t_n)+"_"+to_string(t_k)+"."+index_name;
    index_type pi;
    if ( !load_from_file(pi, idx_file) ) {
        vector<uint64_t> keys = load_keys(hash_file);
        {   
//            my_timer<> t("index construction");
            auto temp = index_type(keys);
//            std::cout<<"temp.size()="<<temp.size()<<std::endl;
            pi = std::move(temp);
        }
        store_to_file(pi, idx_file);
    }
    pi.info(); 
}

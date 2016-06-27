#include "multi_idx/multi_idx.hpp"
#include "multi_idx/multi_idx_red.hpp"
#include "multi_idx/multi_idx_helper.hpp"
#include "chi/chi_idx.hpp"
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;
using namespace phi;
using namespace sdsl;
using namespace multi_index;

//#define N 4
//#define K 7
//#define INDEX_NAME "triangle_clusters_bv_split_threshold100_red"
//#define INDEX_TYPE multi_idx_red<triangle_clusters_binvector_split_threshold<100>,t_k>

const string index_name = INDEX_NAME;
string filtered_index_name;

constexpr uint8_t t_n = N;
constexpr uint8_t t_k = K;
typedef INDEX_TYPE      index_type;

template std::map<uint64_t,uint64_t> get_dist<bucket_accumulator>(const index_type&);
template std::map<uint64_t,uint64_t> get_dist<cluster_accumulator>(const index_type&);
template std::map<uint64_t,uint64_t> get_dist<error_accumulator>(const index_type&);

template<uint8_t t_n, uint8_t t_k, typename t_idx>
struct idx_file_trait{
    static std::string value(){
        return to_string(t_n)+"_"+to_string(t_k)+"."+filtered_index_name;
    }
};


template<uint8_t t_n, uint8_t t_k, typename t_strat>
struct idx_file_trait<t_n, t_k, multi_idx_red<t_strat, t_k>>{
    static std::string value(){
        return to_string(t_n)+"."+filtered_index_name;
    }
};


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
    for (size_t i=0; i<key_buf.size(); ++i){
        keys.push_back(key_buf[i]);
    }
    return unique_vec(keys);
}

/*
template<typename t_acc, typename t_index>
void output(t_index& pi, string id, string hash_file){
    auto res = t_acc(pi);
    for(auto x : res){
        cout << index_name << " " << hash_file << " " << N << " "<< id <<" " << x.first << " " << x.second << endl; 
    }
}
*/


int main(int argc, char* argv[]){

    filtered_index_name = index_name;

    if ( argc < 2 ) {
        cout << "Usage: ./" << argv[0] << " hash_file" << endl;
        return 1;
    }

    string hash_file = argv[1];
    string idx_file = hash_file + "." + idx_file_trait<t_n,t_k,index_type>::value();
    index_type pi;

    if ( !load_from_file(pi, idx_file) ) {
        vector<uint64_t> keys = load_keys(hash_file);
        {   
            pi = index_type(keys,true);
        }
        store_to_file(pi, idx_file);
    }

    {
        auto res = get_dist<bucket_accumulator>(pi);
        for(auto x : res){
            cout << index_name << " " << hash_file << " " << N << " b " << x.first << " " << x.second << endl; 
        }
    }

    {
        auto res = get_dist<cluster_accumulator>(pi);
        for(auto x : res){
            cout << index_name << " " << hash_file << " " << N << " c " << x.first << " " << x.second << endl; 
        }
    }

    {
        auto res = get_dist<error_accumulator>(pi);
        for(auto x : res){
            cout << index_name << " " << hash_file << " " << N << " e " << x.first << " " << x.second << endl; 
        }
    }

//    output<dist<bucket_accumulator>>(pi, string("b"), hash_file);
//    output<dist<cluster_accumulator>>(pi, string("c"), hash_file);

    return 0;
}

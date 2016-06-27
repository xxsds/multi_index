#include "multi_idx/multi_idx.hpp"
#include "multi_idx/multi_idx_red.hpp"
#include "multi_idx/linear_scan.hpp"
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <vector>
#include <sstream>

using namespace std;
using namespace sdsl;
using namespace multi_index;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

const string index_name = INDEX_NAME;
string filtered_index_name;

template<class duration_type=std::chrono::seconds>
struct my_timer{
    string phase;
    time_point<timer> start;

    my_timer() = delete;
    my_timer(string _phase) : phase(_phase) {
        std::cout << "Start phase ``" << phase << "''" << std::endl;
        start = timer::now();
    };
    ~my_timer(){
        auto stop = timer::now();
        std::cout << "End phase ``" << phase << "'' ";
        std::cout << " ("<< duration_cast<duration_type>(stop-start).count() << " ";
        std::cout << " ";
        std::cout << " elapsed)" << std::endl;
    }
};

vector<uint64_t> unique_vec(vector<uint64_t> v){
    if(v.size()>0){
        std::sort(v.begin(), v.end());
        auto end = std::unique(v.begin(), v.end());
        v.resize(end-v.begin());
    }
    return v;
}

template<uint8_t t_b, uint8_t t_k, typename t_idx>
struct idx_file_trait{
    static std::string value(std::string hash_file){
        return hash_file + "." + to_string(t_b)+"_"+to_string(t_k)+"."+filtered_index_name;
    }
};


template<uint8_t t_b, uint8_t t_k, typename t_strat>
struct idx_file_trait<t_b, t_k, multi_idx_red<t_strat, t_k>>{
    static std::string value(std::string hash_file){
        return hash_file + "." + to_string(t_b)+"."+filtered_index_name;
    }
};

template<uint8_t t_b, uint8_t t_k>
struct idx_file_trait<t_b, t_k, linear_scan<t_k>>{
    static std::string value(std::string hash_file){
        return hash_file;
    }
};

vector<uint64_t> load_keys(string hash_file)
{
    vector<uint64_t> keys;
    int_vector_buffer<64> key_buf(hash_file, ios::in, 1<<20, 64, true);
    for (size_t i=0; i<key_buf.size(); ++i){
        keys.push_back(key_buf[i]);
    }
    return unique_vec(keys);
}


void warmup_core_and_cache(){
    std::vector<uint64_t> v(1ULL<<23, 0xABCDABCDABCDABCDULL); // generate 8*4M = 64 MB data
    for(size_t i=0; i<v.size(); ++i){
        v[i] = v[i]-7*i;
    }
    uint64_t cnt = 0;
    for(size_t i=0; i<v.size()/2; ++i){
        cnt += v[i]^v[v.size()-1-i];
    }
    std::cout<<"warmup "<<cnt<<std::endl;
}

int main(int argc, char* argv[]){
    constexpr uint8_t t_b = BLOCKS;
    constexpr uint8_t t_k = K;
    typedef INDEX_TYPE      index_type;
    filtered_index_name = index_name;
    std::string in_fix = "_simd";
    auto pos = index_name.find(in_fix);
    if ( pos != std::string::npos ) {
        filtered_index_name = index_name.substr(0, pos) + 
                              index_name.substr(pos+in_fix.size());
    }

    if ( argc < 2 ) {
        cout << "Usage: ./" << argv[0] << " hash_file [query_file] [search_only] [check_mode] [print_header_for_search_only] [parallel_construction]" << endl;
        cout << " search_only: 0=No (default); 1=Yes" << endl;
        cout << " check_mode: 0=No (default); 1=Yes" << endl;
        cout << " print_header_for_search_only: 0=No (default); 1=Yes" << endl;
        cout << " parallel construction: 0=No (default); 1=Yes" << endl;
        return 1;
    }

    string hash_file = argv[1];
    string idx_file = idx_file_trait<t_b,t_k,index_type>::value(hash_file);
    index_type pi;

    {
        ifstream idx_ifs(idx_file);
        if ( !idx_ifs.good() ){
             vector<uint64_t> keys = load_keys(hash_file);
            {   
    //            my_timer<> t("index construction");
//                auto temp = index_type(keys, async);
                auto temp = index_type(keys, false);
    //            std::cout<<"temp.size()="<<temp.size()<<std::endl;
                pi = std::move(temp);
            }
            store_to_file(pi, idx_file);
            write_structure<HTML_FORMAT>(pi, idx_file+".html");
        }
    }

    string line;
    while ( argc >=3 or std::getline (std::cin, line, ':') ) {
        string qry_file = "";
        size_t search_only = 0;
        size_t check_mode = 0;
        bool print_info = false;
        bool async = false;
        
        if ( argc >= 3 ) {
            qry_file = argv[2];
            if ( argc > 3 ) { search_only = stoull(argv[3]); }
            if ( argc > 4 ) { check_mode  = stoull(argv[4]); }
            if ( argc > 5 ) { print_info  = stoull(argv[5]); }
            if ( argc > 6 ) { async       = stoull(argv[6]); }
        } else {
            stringstream linestream(line); 
            linestream >> qry_file;
            if ( qry_file.size() == 0 )
                continue;
            if ( linestream ) linestream >> search_only;
            if ( linestream ) linestream >> check_mode;
            if ( linestream ) linestream >> print_info;
            if ( linestream ) linestream >> async;
        }

        

        if ( pi.size() == 0 and !load_from_file(pi, idx_file) ) {
            std::cout<<"ERROR. Index size == 0 and index could not be loaded from disk."<<std::endl;
            return 1;
        }

        warmup_core_and_cache();

        int_vector<64> qry;
        if ( !load_vector_from_file(qry, qry_file, 8) ){
            cout << "Error: Could not load query file " << qry_file << "." << endl;
            return 1;
        } else {
            cout << "Loaded " << qry.size()<< " queries form file "<<qry_file << endl;
        }

        if(!search_only or (search_only and print_info)) {
            cout << "# hash_file = " << hash_file << endl;
            cout << "# hashes = " << pi.size() << endl;
            cout << "# qry_file = " << qry_file << endl;
            cout << "# index = " << index_name << endl;
            cout << "# b = " << (size_t)t_b << endl;
            cout << "# k = " << (size_t)t_k << endl;
            cout << "# index_size_in_bytes = " << size_in_bytes(pi) << endl;

        //    vector<uint64_t> pat;
        //    {        
        //        int_vector_buffer<64> key_buf(hash_file, ios::in, 1<<20, 64, true);
        //        pat.push_back( 17615424867151718988ULL );
        //        for (size_t i=0; i < 10000 and i < key_buf.size(); ++i) {
        //            pat.push_back(key_buf[i]);
        //        }
        //    }
            cout << "# queries = " << qry.size() << endl;
        }
        if ( !check_mode ) {
            size_t check_cnt = 0;
            size_t match_cnt = 0;
            size_t unique_cnt = 0;
            
            if(!search_only) {
                {
                  auto start = timer::now();
                  for (size_t i=0; i<qry.size(); ++i){
                      auto result = pi.match(qry[i]);
                      check_cnt += get<1>(result);
                      match_cnt += get<0>(result).size();
                      unique_cnt += unique_vec(get<0>(result)).size();
#ifdef STATS
                      // cum_clusters
                      cout << "@ "<< (int)(std::tuple_element<0,decltype(pi.m_idx)>::type::threshold);
                      cout << " " << stat_counter["cum_clusters"];
                      cout << " " << stat_counter["cum_sub_queries"];
                      cout << " " << stat_counter["cum_clusters_checked"];
                      cout << " " << stat_counter["cum_clusters_visited"];
                      cout << " " << stat_counter["cum_cluster_sizes"];
                      cout << " " << stat_counter["cum_cluster_survivors"];
                      cout << " " << stat_counter["early_exits"];
                      cout << endl;
                      stat_counter.clear();
#endif
                    }
                  auto stop = timer::now();
                  cout << "# time_per_full_query_in_us = " << duration_cast<chrono::microseconds>(stop-start).count()/(double)qry.size() << endl;
                }
                cout << "# check_cnt_full = " << check_cnt << endl;
                cout << "# check_match_full = " << match_cnt << endl;
                cout << "# check_unique_matches_full = " << unique_cnt << endl;
                cout << "# matches_per_query = " << ((double)match_cnt)/qry.size() << endl;
                cout << "# unique_matches_per_query = " << ((double)unique_cnt)/qry.size() << endl;
                cout << "# candidates_per_query = " << ((double)check_cnt)/qry.size() << endl;
                cout << "# check_cnt_search = " << check_cnt << endl;
            } else {
                check_cnt = 0;
                {
                  auto start = timer::now();
                  for (size_t i=0; i<qry.size(); ++i){
                      check_cnt += get<1>(pi.match(qry[i], true));
                    }
                  auto stop = timer::now();
                  cout << "# time_per_search_query_in_us = " << duration_cast<chrono::microseconds>(stop-start).count()/(double)qry.size() << endl;
              }
            }
        }
        else
        {
             auto keys = load_keys(hash_file);
             if (argc > 5){
                qry = int_vector<64>(1, stoull(argv[5]));
                cout << "Check only one key" << endl;
             }
             cout << "Checking results "<< endl;
             for (size_t i=0; i<qry.size(); ++i){

                auto res = get<0>(pi.match(qry[i]));
                res = unique_vec(res);
/*                
                std::cout<<"+++++++++++++++"<<std::endl;
                for(size_t i=0; i<res.size(); ++i){
                    if(res[i]==9649569770737409418ULL){
                        std::cout<<"res["<<i<<"]="<<res[i]<<std::endl;
                    }
                }
                cout<<"3="<<bitset<64>(3)<<endl;
                std::cout<<"+++++++++++++++"<<std::endl;
*/
                vector<uint64_t> res_check;
                for (size_t j=0; j<keys.size(); ++j){
                    if ( bits::cnt(keys[j] ^ qry[i]) <= t_k ) {
                        res_check.push_back(keys[j]);
                    }
                }
                res_check = unique_vec(res_check);

                bool error = false;
                if (res.size() != res_check.size()){
                    error = true;
                    cout<<"ERROR: res.size()="<<res.size()<<"!="<<res_check.size()<<"=res_check.size()"<<endl;
                    cout<<"key="<<bitset<64>(qry[i])<<" ("<<qry[i]<<")"<<endl;

                    size_t j0=0, j1=0;
                    while ( j0+j1 < res.size() + res_check.size() ) {
                        if ( j0 < res.size() and j1 < res_check.size() ) {
                            if ( res[j0] == res_check[j1] ) {
                                ++j0; ++j1;
                            } else {
                                cout<<"res="<<bitset<64>(res[j0])<<" ("<<res[j0]<<") d="<< bits::cnt(res[j0]^qry[i]) <<" j0="<<j0<<endl;
                                cout<<"cck="<<bitset<64>(res_check[j1])<<" ("<<res_check[j1]<<") d="<< bits::cnt(res_check[j1]^qry[i])<<" j1="<<j1<<endl;
                                if ( res[j0] < res[j1] )
                                    ++j0;
                                else
                                    ++j1;
                            }
                        } else if ( j0 == res.size() ){
                             cout<<"res processed"<<endl;
                             cout<<"cck="<<bitset<64>(res_check[j1])<<" ("<<res_check[j1]<<") d="<< bits::cnt(res_check[j1]^qry[i])<<endl;
                            ++j1;
                        } else { // j1  == res.size()
                            cout<<"res="<<bitset<64>(res[j0])<<" ("<<res[j0]<<") d="<< bits::cnt(res[j0]^qry[i]) <<endl;
                             cout<<"cck processed"<<endl;
                            ++j0;
                        }
                    }
                    break;
/*
                    for(size_t j=0; j<res.size(); ++j){
                        cout<<"res="<<bitset<64>(res[j])<<" ("<<res[j]<<") d="<< bits::cnt(res[j]^qry[i]) <<endl;
                    }
                    for(size_t j=0; j<res_check.size(); ++j){
                        cout<<"cck="<<bitset<64>(res_check[j])<<" ("<<res_check[j]<<") d="<< bits::cnt(res_check[j]^qry[i])<<endl;
                    }
*/                    
                } 
                for (size_t j=0; j<res.size(); ++j){
                    if (res[j]!=res_check[j]){
                        cout<<"ERROR: res["<<j<<"]="<<res[j]<<"!="<<res_check[j]<<"=res_check["<<j<<"]"<<endl;
                        cout<<"key="<<qry[i]<<endl;
                        error = true;
                        break;
                    }
                }
                if (!error){
                    if ( res.size() > 1 )
                        cout << "*";
                    else
                        cout << ".";
                    cout.flush();
                } else {
                    return 1;
                }
            }       
            cout << endl;
        }


        if ( argc >= 3 ) {
            break;
        }
    }
    return 0;
}

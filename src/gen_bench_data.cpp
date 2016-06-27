#include <iostream>
#include <algorithm>  
#include <vector>
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;


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


void store_keys(vector<uint64_t> v, string &filename) {
  ofstream out(filename, std::ofstream::binary);
  for (size_t i=0; i < v.size(); ++i){
      out.write((char*)&v[i], sizeof(v[i]));
  }
  out.close();
}


int main(int argc, char* argv[]){
  
  uint64_t n_queries = 10000;
  if ( argc < 2 ) {
      cout << "Usage: ./" << argv[0] << " hash_file [num_queries]" << endl;
      return 1;
  }
  
  string hash_file = argv[1];
  if ( argc > 2 ) {
      n_queries = stoull(argv[2]);    
  }
  
  vector<uint64_t> keys = load_keys(hash_file);
  random_shuffle (keys.begin(), keys.end());
  
  /* Select and write real but not existing keys as queries */
  vector<uint64_t> queries (n_queries, 0);
  copy(keys.end()-n_queries, keys.end(), queries.begin());
  string query_file = hash_file + string(".real.query");
  store_keys(queries, query_file);
  keys.erase(keys.end()-n_queries, keys.end());
  
  /* Select and write existing keys as queries */
  copy(keys.begin(), keys.begin()+n_queries, queries.begin()); 
  query_file = hash_file + string(".existing.query");
  store_keys(queries, query_file);
/*
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis;
  
  // Select and random keys as queries
  for(size_t i = 0; i < n_queries; i++) 
    queries[i] = dis(gen);
  query_file = hash_file+string(".random.query");
  store_keys(queries, query_file);
*/  
  /* Write remaining keys as data */
  hash_file = hash_file + string(".data");
  store_keys(keys, hash_file);
 
}

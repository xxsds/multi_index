#include <iostream>
#include <vector>
#include "sdsl/int_vector.hpp"

using timer = std::chrono::high_resolution_clock;

template<size_t BITS>
  struct scanning {
    uint64_t sum;
    
    scanning(size_t n, size_t m, size_t runs) {
      sdsl::int_vector<BITS> keys(n);
      
      for(size_t i = 0; i < n; ++i) 
        keys[i] = i;
      
    
      std::vector<uint64_t> queries;
      queries.reserve(runs);
      srand(42); // make everything dterministic
      for(size_t i = 0; i < runs; ++i) 
        queries.push_back(rand() % (n-m-1));
      
      auto start = timer::now();
      sum = 0;
      for(auto x: queries) {
        auto begin = keys.begin() + x;
        auto end = keys.begin() + x + m;
        for (auto it = begin; it != end; ++it)
          sum += *it;
      }
      auto stop = timer::now();
      std::cout << BITS << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count()/(double)(runs*m) << " " << sum << std::endl;   
    }
  };

template<uint8_t B0, uint8_t B1>
void scanning_range(size_t n, size_t m, size_t runs){
    if ( B0 <= B1 ) {
        { scanning<B0<B1 ? B0 : B1>(n, m, runs); }
        scanning_range<B0+(B0<=B1), B1>(n, m, runs);
    }
}

int main() {//int argc, char** argv) {
  
    size_t runs = 100000;
    size_t n = (1ULL << 30);
    size_t m = (1ULL << 14);
    
    std::cout << "# number of keys " << n << " run length " << m << " number of runs " << runs << std::endl;
    std::cout << "# bits,nanosecs_per_key" << std::endl; 
    scanning_range<8, 64>(n, m, runs);
}

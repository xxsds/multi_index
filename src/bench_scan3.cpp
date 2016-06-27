#include <iostream>
#include <vector>
#include "sdsl/int_vector.hpp"
#include "multi_idx/perm.hpp"

using timer = std::chrono::high_resolution_clock;

template<size_t STEP>
  struct scanning {
    uint64_t sum;
    
    scanning(size_t n, size_t m, size_t runs) {
      sdsl::int_vector<32> keys(n);
      uint64_t mask = n/2;
      
      for(size_t i = 0; i < n; ++i) 
        keys[i] = i;
          
      std::vector<uint64_t> queries;
      queries.reserve(runs);
      srand(42); // make everything deterministic
      for(size_t i = 0; i < runs; ++i) 
        queries.push_back(rand() % (n-m-1));
      
      auto start = timer::now();
      sum = 0;
      size_t iter = 0;
      for(auto x: queries) {
        auto it = keys.begin() + x;
        auto end = keys.begin() + x + m;
        
        if(STEP == 1) {
          for (; it != end; ++it) {
            if ( sdsl::bits::cnt(*it^mask) <= 3 ) 
              sum += perm<4,1>::mi_rev_permute[0](*it);
            iter += 1;
          }      
        } else {  // Skip STEP elements every STEP elements
          while(it < end) {
            auto end2 = it + STEP;
            for (; it != end2; ++it) {
              if ( sdsl::bits::cnt(*it^mask) <= 3 ) 
                sum += perm<4,1>::mi_rev_permute[0](*it);
              iter += 1;
            }
            it += STEP;
          }
        }
      }
      auto stop = timer::now();
      size_t v = STEP == 1 ? 0 : STEP;
      std::cout << v << " " << std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count()/(double)(runs) << " " << sum << " " << iter << std::endl;   
    }
  };

template<size_t S0, size_t S1>
void scanning_range(size_t n, size_t m, size_t runs){
    if ( S0 <= S1 ) {
        { scanning<S0<S1 ? S0 : S1>(n, m, runs); }
        scanning_range<S0*2, S1>(n, m, runs);
    }
}

int main() {//int argc, char** argv) {
  
    size_t runs = 100000;
    size_t n = (1ULL << 30);
    size_t m = (1ULL << 14);
    
    std::cout << "# number of keys " << n << " run length " << m << " number of runs " << runs << std::endl;
    std::cout << "# step,nanosecs_per_key" << std::endl; 
    scanning_range<1, 256>(n, m, runs);
}

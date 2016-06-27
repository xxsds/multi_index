import itertools
import sys
import random
from subprocess import call

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def TAoCPpermutation(n,k):
  """
  Implements solution presented in The Art of Computer Programming III, Exercise 6.5-1 
  to generate a set of (n choose k) perms such that, 
  for any t-subset of [n], there exists a permutation prefixed by the subset, 
  where t \leq k or t \geq n-k 
  """
  perms = []
  for subset in itertools.combinations(range(n), k):
     A = []; B = []; C = []; min = 0; j = 0; up = 0
     for i in xrange(n):
        if(j>=k or i != subset[j]):
            B.append(i)
            up +=1
        else:
            up -=1
            j += 1
            if(up < min):
                min = up
                B.append(i)
            else:
                A.append(i)
                C.append(B.pop())
     perms.append(A+B+C)
  return perms
  

def get_k_sets(perm, k):
    curr_set = []
    for i in xrange(n):
        l = (perm*2)[i: i+k]
        curr_set.append(tuple(sorted(l)))
    return curr_set

def is_intersetion_empty(big, small):
    for x in small:
        if x in big:
            return False
    return True

print "Usage: python", sys.argv[0], "n k"

n = int(sys.argv[1]) # number of blocks
k = int(sys.argv[2]) # number of blocks that have to be consecutive at least once in some permutation

sys.stderr.write(str(n)+"\n")
sys.stderr.write(str(k)+"\n")
sys.stderr.write(str(TAoCPpermutation(n,k))+"\n")

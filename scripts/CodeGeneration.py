import sys
import ast
import subprocess

def write_function(perm, function_name, orig_perm, hpp_file):
    print orig_perm
    hpp_file.write("\t// " + str(orig_perm)+ "\n")
    hpp_file.write("\tstatic constexpr uint64_t " + function_name +"(uint64_t x) {\n")
    proc = subprocess.Popen(["./calcperm", " ".join([str(x) for x in perm])], stderr=subprocess.PIPE)
    code = proc.stderr.read().replace("\n","\n\t\t")
    hpp_file.write("\t\t" + code)
    hpp_file.write("return x;\n")
    hpp_file.write("\t}\n\n") 


def write_code_set(selected_perms, permute_name, rev_permute_name, struct_name, hpp_file, cpp_file):
    # write the set of perms perms
    fun_vec_name = permute_name + "_permute"
    rev_fun_vec_name = rev_permute_name + "_permute"
    permute_array_name = permute_name + "_perms"
    perm_function_name = permute_name + "_{}" # permutation_id
    revperm_function_name = rev_permute_name + "_{}" # permutation_id
    perm_cnt = str(len(selected_perms))
    block_cnt = str(len(sizes))
    array_str = "std::array<uint8_t,"+block_cnt+">";

    id = 0
    for perm in selected_perms:
        permutation = [] # it says where each of the 64 bits has to go accordingly to perm
        for x in perm:
            permutation += bits[x]
        write_function(permutation, perm_function_name.format(id), perm, hpp_file)
  
        inverse = [0] * len(permutation)
        for i, p in enumerate(permutation):
            inverse[p] = i
        write_function(inverse, revperm_function_name.format(id), perm, hpp_file) 
        id += 1
    
    hpp_file.write("\tstatic constexpr t_fun_vec_"+permute_name+ " " + fun_vec_name + " = {\n")
    for id in xrange(len(selected_perms)-1):
        hpp_file.write("\t\t&" + struct_name + "::" + perm_function_name.format(id) + ",\n")
    hpp_file.write("\t\t&" + struct_name + "::"+ perm_function_name.format(len(selected_perms)-1) + "\n};\n\n")


    hpp_file.write("\tstatic constexpr t_fun_vec_"+permute_name+ " " + rev_fun_vec_name + " = {\n")
    for id in xrange(len(selected_perms)-1):
        hpp_file.write("\t\t&" + struct_name + "::"+ revperm_function_name.format(id)+",\n")
    hpp_file.write("\t\t&" + struct_name + "::"+ revperm_function_name.format(len(selected_perms)-1)+"\n};\n\n")

    hpp_file.write("\tstatic constexpr std::array<"+array_str+","+perm_cnt+"> " + permute_array_name + "={{\n")
    for id in xrange(len(selected_perms)-1):
        hpp_file.write("\t\t{" + ", ".join(str(x) for x in selected_perms[id]) + "},\n")
    hpp_file.write("\t\t{" + ", ".join(str(x) for x in selected_perms[-1]) + "}\n}};\n\n")

    hpp_file.write("\tstatic constexpr "+array_str+ " " + fun_vec_name + "_block_sizes = {{\n\t" + ", ".join(str(x) for x in sizes) + "\n}};\n")

    widths = []
    widths_sums =[]
    for perm in selected_perms:
        l = []
        for el in perm:
            l.append(sizes[el])
        p = [0]*len(l)
        for i in xrange(1, len(l)):
            p[i] = p[i-1] + l[i-1]

        widths.append("\t{" + ",".join([str(x) for x in l]) + "}")
        widths_sums.append("\t{" + ",".join([str(x) for x in p]) + "}")

    # write width of each block in each permutation    
    hpp_file.write("\tstatic constexpr std::array<"+array_str+","+perm_cnt+"> " + fun_vec_name + "_block_widths = {{ // width of the blocks in the 64-bit word \n")
    hpp_file.write(",\n".join(widths) + "\n\t}};\n\n")

    # write prefix sums of blocks' width in each permutation
    hpp_file.write("\tstatic constexpr std::array<"+array_str+","+perm_cnt+"> " + fun_vec_name + "_block_widths_sums = {{ // prefix sums of widths of the blocks in the 64-bit word \n")
    hpp_file.write(",\n".join(widths_sums) + "\n\t}};\n\n")
    
    # write vectors of functions pointers
    cpp_file.write("constexpr "+ struct_name + "::t_fun_vec_"+permute_name+" " + struct_name + "::" + fun_vec_name + ";")
    cpp_file.write("constexpr "+ struct_name + "::t_fun_vec_"+permute_name+" " + struct_name + "::" + rev_fun_vec_name + ";")     
    cpp_file.write("constexpr std::array<"+array_str+","+perm_cnt+"> " + struct_name + "::" + permute_array_name + ";\n")

    cpp_file.write("\nconstexpr "+array_str+" " + struct_name + "::" + fun_vec_name + "_block_sizes;\n");
    cpp_file.write("\nconstexpr std::array<"+array_str+","+perm_cnt+"> " + struct_name + "::" + fun_vec_name + "_block_widths;\n")
    cpp_file.write("\nconstexpr std::array<"+array_str+","+perm_cnt+"> " + struct_name + "::" + fun_vec_name + "_block_widths_sums;\n")

if len(sys.argv) < 3:
    print "Usage: python", sys.argv[0], "permutation_filename code_filename"
    print "Make sure you compiled calcperm with g++ -o calcperm calcperm.cpp"
    sys.exit(1)

# Read permutation file 
# Input format:
#   - n 
#   - k
#   - Multi_index permutations

perm_file = open(sys.argv[1], "r")
hpp_file = open(sys.argv[2] + ".hpp", "w")
cpp_file = open(sys.argv[2] + ".cpp", "w")

n = int(perm_file.readline())
k = int(perm_file.readline())
mi_perms = ast.literal_eval(perm_file.readline())

print "Multi Index permutations:"
for perm in mi_perms:
  print perm

# Generate permuting/inverting code

# block sizes are w_b or w_b +1. There are exactly w_b_rem of the latter size. 
w_b = 64/n  
w_b_rem = 64%n
 
#sizes = [w_b+1]*w_b_rem + [w_b]*(n-w_b_rem) 
sizes = [w_b]*(n-w_b_rem) + [w_b+1]*w_b_rem
print "Sizes of the blocks are:"
print sizes

if (sum(sizes) != 64):
    print "ERROR: sizes do not sum to 64"

bits = []
start = 0
for x in sizes:
    end = start + x
    bits.append([x for x in xrange(start, end)])
    start = end
print bits

# generate code for the perms and their inverted.

# Header
hpp_file.write("#pragma once\n\n#include <vector>\n#include \"multi_idx/perm.hpp\"\n#include \"multi_idx/aux_bits.hpp\"\t\n")

struct_name = "perm<" +str(n) + "," + str(k)+">"
hpp_file.write("template<>\nstruct " + struct_name + " {\n\n")
hpp_file.write("\ttypedef uint64_t (*perm_fun_t)(uint64_t);\n")
hpp_file.write("\ttypedef std::array<perm_fun_t,"+str(len(mi_perms))+"> t_fun_vec_mi;\n\n")
hpp_file.write("\tstatic constexpr uint8_t max_dist = " + str(n-k) + ";\n")
hpp_file.write("\tstatic constexpr uint8_t match_len = " + str(k) + ";\n\n")

cpp_file.write("#include \"multi_idx/"+sys.argv[2]+".hpp\"\n\n")
cpp_file.write("// initialize static const members\n\n")

# generate functions 
write_code_set(mi_perms, "mi", "mi_rev", struct_name, hpp_file, cpp_file)

# Footer
hpp_file.write("};")

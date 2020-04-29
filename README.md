# MinCNE
To compile:
    g++ -o MinCNE MinCNE.cpp ./edlib/src/edlib.cpp -std=c++11 -O3 -march=native -lz -I ./edlib/include/
To run:
    ./MinCNE <m> <input.txt> <q> <b> <n> <th> <minLen>
    Where m = minhash size (default: 200)
          q = q-gram size (default: 13)
          b = bucket size (default: 5)
          n = number of sequences
          th = threshold similarity 
          minLen = minimum length of CNE
    input.txt contains the lines in the following format 
        species1.chr.fa,chr#,start,end
        species2.chr.fa,chr#,start,end
          

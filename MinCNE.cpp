/*
 
*/
#include <sstream>
#include <stack> 
#include <omp.h>
#include "edlib.h"
#include <unordered_map>
#include "metrohash64.cpp"
#include <string>
#include <cstring>
#include <stdint.h>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <openssl/md5.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <stdlib.h>
#include "ntHashIterator.hpp"
KSEQ_INIT(gzFile, gzread)

using namespace std;

void usage(const string s)
{
    fprintf(stderr, "usage: ./minhash <k>  <fasta> <kmerLen> <bucketsize> <min_element_cluster> <threshold_similarity>\n");
    if (s != "") fprintf(stderr, "%s\n", s.c_str()); 
    exit(1);
}

int getPos(int x, vector<int> y, int *new_pos){
  int target = 0;   
  *new_pos = x;
  for(int i=0; i<y.size(); i++){
    target += y.at(i);
    if(i != 0) *new_pos = *new_pos - y.at(i-1); 
    if(x < target) return i;
  }
}

int main(int argc, char **argv)
{
  vector <string> seqnames;             // Sequence names
  vector <string> seqns;                // Sequences
  vector <unsigned char *> min_hashes;  // The minimum hashes for each file.
  int k;                                // The number of hashes
  int bbh;                              // Bytes per hash
  int hash_buf_size;                    // Size of the hash buffers (k*bbh) padded to 16
  unsigned int ff;                      // An integer that holds 0xffffffff
  unsigned char *hash;                  // Where we calculate the hashes for each string.
  ifstream f;
  string s1;
  int findex;
  int i, j, p, sz;
  double Intersection;                  
  int kmerLen;                          // q-gram size
  int bsize;                            // number of buckets

  /* Read the command line arguments. */

  if (argc <= 1) usage("");

  k = atoi(argv[1]); // size of minhash 
                 
  if (k <= 0) usage("k must be a number > 0");
  
  /* read the text file that contains the name of fasta files, chromosome, start, end */
  ifstream infile; 
  string infile_name(argv[2]);
  infile.open(infile_name);
  string infile_line;
  vector<int> n0_kmers;
  int kk = atoi(argv[7]);

  vector<int> SP;
  vector<int> LP;
  //set<string> SET_KMERS;

  /* read the fasta files and store all k-mers */
  int num_files = 0;

  vector<string> SEQQ;
  while ( getline (infile, infile_line) ){
    std::stringstream ss(infile_line);
    std::string token;
    int no_token = 0;
    string infile_name;
    string chr_name;
    int start_pos, end_pos;
    while(std::getline(ss, token, ',')) {
      if(no_token == 0) infile_name = token;
      else if(no_token == 1) chr_name = token;
      else if(no_token == 2) { start_pos = atoi(token.c_str()); SP.push_back(start_pos); }
      else if(no_token == 3) { end_pos = atoi(token.c_str()); LP.push_back(end_pos); }
      no_token++;    
    }
    num_files++;
    gzFile fp;
    kseq_t *seq;
    int l;
    ofstream seqFile;
    string seqFileName = infile_name+"actual_seq.fa";
    seqFile.open (seqFileName);
    ofstream myfile3;
    fp = gzopen(infile_name.c_str(), "r"); //file name
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0) { 
      string seqname(seq->name.s);
      if(seqname == chr_name){
        string seqn(seq->seq.s);
        //SEQQ.push_back(seqn);
        string actual_seq = seqn.substr(start_pos, end_pos-start_pos+1);
        transform(actual_seq.begin(), actual_seq.end(), actual_seq.begin(), ::toupper);
        int len = actual_seq.length();
        cout << len << endl;
        SEQQ.push_back(actual_seq);
        seqFile << ">" << infile_name << len << "\n" << actual_seq << "\n";
        for(int i=0; i<len-kk+1; i++){
          string kkmer = actual_seq.substr(i, kk);
          seqns.push_back(kkmer);
        }
        n0_kmers.push_back(len-kk+1);
      }
      if(seqname == chr_name) break;
    }
    kseq_destroy(seq);
    gzclose(fp);
    seqFile.close();
  }

  cout << "n0_kmers.size: " << n0_kmers.size()<< endl;
  cout << "size: " << seqns.size() << endl;

  map<uint64_t, vector<string>> HASH_to_STRING;

  kmerLen = atoi(argv[3]); // q-gram size
  
  map<string, uint64_t> hashV;
   
  vector<vector<uint64_t>> MHASH;
  
  /* Read the data sets.  For each value, you're going to calculate the k hashes
     and then update the minimum hashes for the data set. */
 
  for (findex = 0; findex < seqns.size(); findex++) {
    uint64_t hash1 = 0;
    vector<uint64_t> mmhash;
    string ss = seqns.at(findex);
    ntHashIterator itr(ss, k, kmerLen);
    int t = 0;
    while (itr != itr.end()) {
      for(int u = 0; u < k; u++){
        hash1 = (*itr)[u];
        if(t == 0) mmhash.push_back(hash1);
        else if(hash1 < mmhash.at(u)) mmhash[u] = hash1;
      } 
      ++itr;
      ++t;
    }
    MHASH.push_back(mmhash);
  }

  string s = "";
  cout << "MHASH.size(): " << MHASH.size() << endl;
  
  bsize = atoi(argv[4]);
  s = "";
  cout << "minhash stage completed...bucketing started\n" ;
  
  vector<unordered_map<uint64_t, vector<int>>> MMM;
  
  int initial_count = 0;
  unordered_map<uint64_t, vector<int>> MM1; 
  unordered_map<uint64_t, int> MM_COUNT;
    
  vector<vector<uint64_t>> BHASH; 
  for(int m=0; m<n0_kmers.size(); m++){
    unordered_map<uint64_t, vector<int>> MM;
    cout << endl; 
    for (int i = 0; i < n0_kmers.at(m); i++) {
      vector<uint64_t> bh;
      if(MHASH[initial_count+i].size() == k){
        for (int p = 0; p < k; p += bsize) {
        //for (int p = 0; p < k; p ++) {
          s = "";
          uint64_t hash1 = 0;
          MetroHash64::Hash((uint8_t*)&(MHASH[initial_count+i][p]), sizeof(uint64_t)*bsize, (uint8_t *)&hash1, 0);
          bh.push_back(hash1);
          if(MM.find(hash1) != MM.end()) {
            if(find(MM[hash1].begin(), MM[hash1].end(), initial_count+i) == MM[hash1].end()) 
              MM[hash1].push_back(initial_count+i); 
          }
          else {
            vector<int> v; 
            v.push_back(initial_count+i); 
            MM.insert(make_pair(hash1, v));
            if(MM_COUNT.find(hash1) == MM_COUNT.end()) MM_COUNT.insert(make_pair(hash1,1));
            else MM_COUNT[hash1]++;  
          }
        }
      }
      if(m==0) BHASH.push_back(bh);
    }
    initial_count += n0_kmers.at(m);
    cout << "initial_count(" << m <<"): " << initial_count << endl;
    MMM.push_back(MM);
  }
  
  int cnt1 = 0, cnt2=0, cnt3=0;  
  for (auto& x: MM_COUNT) {
    if(x.second == 2) cnt2++;
    else if(x.second > 2) cnt3++;
    else cnt1++;
  }

  for(int m=0; m<n0_kmers.size(); m++){
    unordered_map<uint64_t, vector<int>> MM = MMM.at(m);
    for (auto& x: MM) {  
      std::vector<int>::iterator it1; 
      it1 = std::unique (x.second.begin(), x.second.end());
      x.second.resize( std::distance(x.second.begin(),it1) ); 
    }
  }
 
  cout << "bucketing completed...\n" ;

  vector<set<int>> CLST2;
  int ind = 0;
  
  ofstream myfile1;  
  myfile1.open ("similar1.txt");
  ofstream myfile2;  
  myfile2.open ("similar2.txt");
  ofstream myfile3;
  myfile3.open ("CNE.txt");
  ofstream myfile4;
  myfile4.open ("CNE.fa");
  ofstream myfile5;
  myfile5.open ("final_contigs.faa");
  ofstream myfile6;
  myfile6.open ("final_contigs_residue.faa");
 
  int cid = 1, cid1 = 1;
  s = "";
  set< set<int> > SSS;
  int no_CNE = 0;
  int *start_pos = new int[num_files];
  int *end_pos = new int[num_files];

  for(int i=0; i<num_files; i++){
    start_pos[i] = -2;
    end_pos[i] = -2;
  }

  for (i = 0; i < n0_kmers.at(0); i++) {
    int myfile1_flag = 0;
    int* arrayFlag = new int[seqns.size()] ();
    //cout << "\r" << i << " processing..." << flush;
    if(MHASH[i].size() == k){  // the k-mer does not contain non-ACTG character 
      vector<int> sseq;
      vector<double> sseq_score;
      set <int> extra_sseq;
      sseq.push_back(i);
      extra_sseq.insert(i);
      sseq_score.push_back(0);

      /* Check for each bucket signature associated with ith k-mer of first sequence */
      //for(int p=0; p<BHASH[i].size(); p++){
      for(int j=1; j<n0_kmers.size(); j++){
        int jjj = i;
        double last_per_sim = 0.0;
        for(int p=0; p<BHASH[i].size(); p++){
          uint64_t hash1 = BHASH[i][p];
          //int jjj = i;
          //double last_per_sim = 0.0;
          if(MM_COUNT[hash1] == atoi(argv[5])) {
            if(myfile1_flag == 0) { myfile1 << i <<" / " << (SP[0] + i) << ": <"; myfile1_flag = 1;}
            //for(int j=1; j<n0_kmers.size(); j++){
            last_per_sim = 0.0;
            if(MMM[j][hash1].size() <= 100){
              for(int k=0; k<MMM[j][hash1].size(); k++){
                int jj = MMM[j][hash1][k];
                if(arrayFlag[jj] == 0){
                  arrayFlag[jj] = 1;
                  EdlibAlignResult result;
                  double per_sim = 0.0;
                  int l1 = seqns.at(i).length(), l2 = seqns.at(jj).length(); 
                  if(l1 < l2)
                    result = edlibAlign(seqns.at(i).c_str(), l1, seqns.at(jj).c_str(), l2, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                  else
                    result = edlibAlign(seqns.at(jj).c_str(), l2, seqns.at(i).c_str(), l1, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                  per_sim = ((double)(result.alignmentLength-result.editDistance)/(double)result.alignmentLength);
                  edlibFreeAlignResult(result);
                  //if(per_sim > atof(argv[6])) {//extra_sseq.insert(jj);
                  if(per_sim > atof(argv[6]) && per_sim > last_per_sim) { jjj = jj; /*sseq.insert(jj);*/ last_per_sim = per_sim;}    
                }
              }
            }
          }
        }
        if(jjj != i) {sseq.push_back(jjj); sseq_score.push_back(last_per_sim);}
      }
      /* checking is done */
      
      if(sseq.size() >= atoi(argv[5])){
        int yyy = 0;
        myfile2 << (SP[0] + i) << ": <";
        int flag = 0;
        if(i == end_pos[0]+1 || i == end_pos[0]+2 || i == end_pos[0]+3) flag = 1; 
        else {
          if(start_pos[0] == -2 && end_pos[0] == -2) {
            start_pos[0] = end_pos[0] = i;
            for(int it=1; it<sseq.size(); it++){
              int sid =sseq.at(it);
              cout << it << " " << sid << endl;
              start_pos[it] = end_pos[it] = sid;
            }
          }
          else{
            cout << "0: " << start_pos[0] << " " << end_pos[0] << " " << (end_pos[0]-start_pos[0]+kk) << endl;
            string s1 = SEQQ[0].substr(start_pos[0], end_pos[0]-start_pos[0]+kk);
            myfile4 << ">" << (no_CNE+1) << "_0:\n" << s1 << "\n";
            myfile3 << (no_CNE+1) << ": " << (start_pos[0] + SP[0]) << "-" << (end_pos[0]+SP[0]+kk) << "," << (end_pos[0]-start_pos[0]+kk);
         
             
            start_pos[0] = i; end_pos[0] = i;
          
            for(int pp=1; pp<num_files; pp++) {
              int new_pos1 = 0, new_pos2 = 0;
              getPos(start_pos[pp], n0_kmers, &new_pos1);
              getPos(end_pos[pp], n0_kmers, &new_pos2);
              double per_sim = 0.0;
                cout << pp << " " << start_pos[pp] << " " << end_pos[pp] << " " << (end_pos[pp]-start_pos[pp]+kk) << endl;
                cout << pp << " " << new_pos1 << " " << new_pos2 << " " << (new_pos2-new_pos1+kk) << endl; //exit(0); 
                string ss1 = SEQQ[pp].substr(new_pos1, new_pos2-new_pos1+kk);
                //cout << "s1:\n" << s1 << endl << "ss1:\n" << ss1 << endl;
               
                EdlibAlignResult result;
                if(s1.length() < ss1.length())
                  result = edlibAlign(s1.c_str(), s1.length(), ss1.c_str(), ss1.length(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                else    
                  result = edlibAlign(ss1.c_str(), ss1.length(), s1.c_str(), s1.length(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
                per_sim = ((double)(result.alignmentLength-result.editDistance)/(double)result.alignmentLength);
                edlibFreeAlignResult(result);
             
                myfile4 << ">" << (no_CNE+1) << "_" << pp << ":\n" << s << "\n";
              myfile3 << "," << (new_pos1 + SP[pp]) << "-" << (new_pos2 + SP[pp] + kk) << "," << (new_pos2-new_pos1+1+kk) << "(" << to_string(per_sim) << ")";
              start_pos[pp] = sseq.at(pp); end_pos[pp] = sseq.at(pp);
              cout << "start and end" << pp << " " << start_pos[pp] << " " << sseq.at(pp) << " " << end_pos[pp] << " " <<  sseq.at(pp) << endl;
            }
          }
          myfile3 <<"\n";
          myfile4 <<"\n\n";
          no_CNE++;
        } 
        for(int it=0; it<sseq.size(); it++){
          int sid =sseq.at(it);
          //if(start_pos[yyy] == -2) {start_pos[yyy] = sid; end_pos[yyy] = sid;}
          if(flag == 1) end_pos[it] = sid; if(it == 2 && start_pos[it] == 1640629) cout << "end_pos[it]: " << end_pos[it] << endl;
          int new_pos = 0;
          int posn = getPos(sid, n0_kmers, &new_pos);
          if(sid != i) myfile2 << " " << (SP[posn] + new_pos) << "(" << to_string(sseq_score.at(it)) << ")";
          myfile5 << ">" << cid << "_" << yyy << ":" << posn << "-" << sid << "\n" << seqns.at(sid) << "\n";
          ++yyy;
        }
        myfile2 << ">\n";
        ++cid;
      }
      else{
        int yyy = 0;
        for(int it=0; it<sseq.size(); it++){
          int sid = sseq.at(it);
          int new_pos = 0;
          int posn = getPos(sid, n0_kmers, &new_pos);
          myfile6 << ">" << cid1 << "_" << yyy << ":" << posn << "-" << sid << "\n" << seqns.at(sid) << "\n";
        }
        ++cid1;
      } 

      std::set<int>::iterator it;

      if(extra_sseq.size() >= 2){
        myfile1 << i <<" / " << (SP[0] + i) << ": <";
        for(it =extra_sseq.begin(); it != extra_sseq.end(); ++it){
          int sid = *it;
          int new_pos;
          int posn = getPos(sid, n0_kmers, &new_pos);
          if(sid != i) myfile1 << " " << (SP[posn] + new_pos);
        }
        myfile1 << ">\n";
      }    
    }
    free(arrayFlag); 
  }

  myfile1.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();
  myfile5.close();
  myfile6.close();
  exit(0);

}


/***********************************/
/* Pre-processing of the reads for */
/* error correction based on a     */
/* consensus procedure with a      */
/* frequence threshold for SNPs    */
/***********************************/

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <set>

#include <seqan/sequence.h>
#include <seqan/file.h>

#define READ_LEN 64

using namespace seqan;
using namespace std;

int main(int argc, char* argv[]){
    if(argc < 3){
        std::cerr << "Usage: confronto <RNA-seq_file> --gene <Gene>" << std::endl;
        std::cerr << std::endl;
        return 1;
    }
  
    long file_dimension = 0;
    ifstream f(argv[1],::std::ios_base::in|::std::ios_base::binary|ios::ate);

    if(f.is_open() ){
        file_dimension = f.tellg();
        f.close();
    }
    else{
        std::cerr << "Unable to open files " << std::endl;
        return 1;	
    }
    long read_size = 0;
    int perc = 0;
    std::fstream data_file;
     data_file.open(argv[1], std::ios_base::in | std::ios_base::binary);
     if(data_file.is_open()){
         string gene(argv[3]);
         std::cerr << "Filtering " << gene << " gene" << std::endl;
         String<char> fasta_tag;
         String<Dna5> read_seq;
         
         //RNA-seq Reads
         clock_t tStart = clock();
         std::cerr << "Processing RNA-seq file..." << std::endl;
         while(!data_file.eof()){
             readMeta(data_file, fasta_tag, Fasta());
             read(data_file, read_seq, Fasta());
             read_size+=((length(fasta_tag)*sizeof(char))+(length(read_seq)*sizeof(char)));
             if ((read_size/(double)file_dimension)*100 - 1 >= perc) {
                 perc++;
                 if(perc%10 == 0){
                     std::cerr << "Processed: " << perc << "%" << std::endl;
                 }
             }
             
             string read_tag;
             assign(read_tag,fasta_tag);
             if(read_tag.find(gene) != string::npos){
                 std::cout << ">" << read_tag << std::endl;
                 std::cout << read_seq << std::endl;
             }
         }
         
         std::cerr << "Processing RNA-seq file...done!" << std::endl;
         std::cerr << "Processing took " << (double)(clock() - tStart)/CLOCKS_PER_SEC;
         std::cerr << " seconds." << std::endl << std::endl;
         data_file.close();
         
     }else{
         std::cerr << "Unable to open RNA-seq file" << std::endl;
         std::cerr << std::endl;
         return 1;
     }
     return 0;
}

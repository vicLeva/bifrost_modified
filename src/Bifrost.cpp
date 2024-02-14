#include "CompactedDBG.hpp"
#include "ColoredCDBG.hpp"

#include <cstdio>
#include <sys/stat.h>

using namespace std;

int main(int argc, char **argv){
    //First part for building bifost index over 1 dataset, writing it in memory and query positive and negative sequences.
    //At end of main() is a commented part for determining bifrost's index size
    std::cout << "Start expe" << std::endl;
    std::cout << "Start 6_1" << std::endl;

    
    string seq;
    CDBG_Build_opt opt;
    opt.k = 31;
    opt.nb_threads = 4;
    opt.verbose = false;

    auto ttot = std::chrono::high_resolution_clock::now();
    opt.filename_seq_in.push_back(string("~/data/AHX_ACXIOSF_6_1_C2FGHACXX.IND4_clean.fastq.gz"));
    CompactedDBG<> cdbg(opt.k);
    cdbg.build(opt);
    std::cout << std::to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (proprocessing+build+inserts)\n";

    cdbg.writeBinary("~/6_1_cdbg");

    std::ifstream infile1("~/data/6_1_reads.fasta");
    ttot = std::chrono::high_resolution_clock::now();
    while (!infile1.eof()) {
        std::getline(infile1, seq);  // skip header line
        std::getline(infile1, seq);

        cdbg.searchSequence(seq, true, false, false, false, false);
    }
    std::cout << std::to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (pos queries)\n";
    
    std::ifstream infile("~/data/neg_reads.fasta");
    ttot = std::chrono::high_resolution_clock::now();
    while (!infile.eof()) {
        std::getline(infile, seq);  // skip header line
        std::getline(infile, seq);

        cdbg.searchSequence(seq, true, false, false, false, false);
    }
    std::cout << std::to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms (neg queries)\n";


    // This part is for computing bifrost's index size after it's been built 
    // To do so : comment upper part, uncomment this part and user /usr/bin/time -v ./bifrost
    /*
    std::cout << "Start expe" << std::endl;
    CDBG_Build_opt opt;
    opt.k = 31;
    opt.nb_threads = 4;
    opt.verbose = false;
    CompactedDBG<> cdbg(opt.k);
    cdbg.readBinary("~/gut_cdbg");
    std::ifstream infile1("~/data/few_queries.fasta");
    while (infile1 >> a) {
        cdbg.searchSequence(a, true, false, false, false, false);
	}
    */
}

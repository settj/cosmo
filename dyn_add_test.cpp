#include <chrono>

#include <boost/date_time.hpp>

#include <sdsl/wavelet_trees.hpp>

// TCLAP
#include "tclap/CmdLine.h"

#include "kmer-counter.hpp"
#include "dyn_debruijn_graph.hpp"

struct parameters_t {
    std::string input_filename = "";
    std::string kmer_filename = "";


    //boost::posix_time::ptime local_time = boost::posix_time::second_clock::local_time();
    //boost::gregorian::date date = local_time.date();
    //std::string date_str = to_simple_string(date);
    //std::string outfile= "addResults" + date_str + ".tsv";
    std::string outfile= "addResults.tsv";
};


void parse_arguments(int argc, char **argv, parameters_t & params)
{
    TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
    //TCLAP::UnlabeledValueArg<std::string> matrix_filename_arg("matrix", "Matrix file. Currently only supports DSK's binary format (for k<=64).", true, "", "matrix_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> kmer_filename_arg("kmers", ".fasta file of kmers to add.", true, "", "kmer file", cmd);

    cmd.parse( argc, argv );

    params.input_filename  = input_filename_arg.getValue();
    params.kmer_filename   = kmer_filename_arg.getValue();
}

int main(int argc, char* argv[]) {
    parameters_t p;
    parse_arguments(argc, argv, p);

    //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    // Can add this to save a couple seconds off traversal - not really worth it.
    std::cerr << "loading dbg "<< p.input_filename  << std::endl;
    debruijn_graph_shifted<> dbg;
    load_from_file(dbg, p.input_filename);
    //input.close();
    std::cerr << "k             : " << dbg.k << endl;
    std::cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    std::cerr << "num_edges()   : " << dbg.num_edges() << endl;
    std::cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
    std::cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl <<endl;

    //debruijn_graph_shifted<> dbg2(dbg.k, dbg.m_node_flags, dbg.m_edges, dbg.m_symbol_ends, dbg.m_alphabet);
    dyn_debruijn_graph<> dyn_dbg(dbg.k, dbg.m_node_flags, dbg.m_edges, dbg.m_symbol_ends, dbg.m_alphabet);


    size_t absent_kmers = 0;
    size_t counter = 0;
    size_t num_kmers; 
    std::set<string> kmers;
    getKmers( num_kmers, dyn_dbg.k, kmers, p.kmer_filename );

    typedef std::chrono::high_resolution_clock Clock;
    //typedef std::chrono::milliseconds ms;
    //typedef std::chrono::duration<double> ddur;
    auto t_start = Clock::now();
    double time_elapsed = 0;

    std::ofstream ofs(p.outfile);
    for (auto it = kmers.begin(); it != kmers.end(); ++it) {
        if (counter % 1000 == 0) {
            std::cerr<<"\r           \r";
            std::cerr<<(double)counter*100/kmers.size()<<"% of kmers were processed";
        }
        std::string kmer = *it;

        t_start = Clock::now();
        dyn_dbg.add(kmer);
        //double query_time = double (Clock::now() - t_start);
        double query_time = (Clock::now() - t_start).count();
        time_elapsed += query_time;

        //ofs<<counter<<"\t"<<kmer<<"\t"<<present<<"\t"<<query_time/CLOCKS_PER_SEC<<endl;
        //if (present == 0)
            //absent_kmers+=1;
        counter+=1;
    }
    ofs.close();

    std::cerr<<"\n100% of kmers were processed\n";
    cout<<kmers.size()-absent_kmers<<" ("<<(double)(kmers.size()-absent_kmers)*100/kmers.size()<<"%) of the kmers were in the graph\n";
    cout<<kmers.size()<< " kmers were queried in "<< time_elapsed/CLOCKS_PER_SEC<<" (s)"<<endl;
    cout<<"Time per operation: "<<(time_elapsed/CLOCKS_PER_SEC)/kmers.size()<<" (s)"<<endl;
}

#include <iostream>
#include <fstream>
#include <string>
//#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "debruijn_graph_shifted.hpp"

struct parameters_t {
  std::string input_filename = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);

  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
}


void dump_edges(const debruijn_graph_shifted<> &dbg) {
  for (size_t i = 0; i < dbg.size(); i++) {
      cout <<  /* << "e:" <<*/ dbg.edge_label(i) << std::endl;
  }
}

int main(int argc, char* argv[]) {
    parameters_t p;
    parse_arguments(argc, argv, p);
    cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
    //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
    // Can add this to save a couple seconds off traversal - not really worth it.
    cerr << "loading dbg "<< p.input_filename  << std::endl;
    debruijn_graph_shifted<> dbg;
    load_from_file(dbg, p.input_filename);
    //input.close();

 
    cerr << "k             : " << dbg.k << endl;
    cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
    cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
    dump_edges(dbg);
 }

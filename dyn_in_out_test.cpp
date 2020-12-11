#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "tclap/CmdLine.h"

#include <boost/dynamic_bitset.hpp>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "debruijn_graph_shifted.hpp"
#include "dyn_debruijn_graph.hpp"

struct io_parameters_t {
  std::string input_filename = "";
};


void parse_arguments(int argc, char **argv, io_parameters_t & params)
{
    TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
    //TCLAP::UnlabeledValueArg<std::string> matrix_filename_arg("matrix", "Matrix file. Currently only supports DSK's binary format (for k<=64).", true, "", "matrix_file", cmd);

    cmd.parse( argc, argv );
    std::cout << "arg0" << argv[0] << std::endl;
    std::cout << "arg1" << argv[1] << std::endl;
    std::cout << "???parsing arg" << input_filename_arg.getValue() << std::endl;
    params.input_filename  = input_filename_arg.getValue();
    std::cout << "done parsing arg" << params.input_filename << std::endl;
}

int main(int argc, char* argv[]) {
    io_parameters_t p;
    parse_arguments(argc, argv, p);

    std::cout << "done parsing arg" << p.input_filename << std::endl;
    //Load in raw dbg
    //debruijn_graph_shifted<> dbg;
    //load_from_file(dbg, p.input_filename);
 
    //Make buffboss structure out of dbg
    dyn_debruijn_graph<> dyn_dbg;
    std::cout << "time to load" << std::endl;
    dyn_dbg.load_from_file(p.input_filename);
    std::cout << "loaded" << std::endl;

    std::cout << "time to insert" << std::endl;
    //add/delete stuff
    dyn_dbg.insert(4, 0);
    std::cout << std::endl;
    std::cout << "time to insert" << std::endl;
    dyn_dbg.insert(4, 1);
    std::cout << std::endl;
    std::cout << "time to insert" << std::endl;
    dyn_dbg.insert(4, 2);
    std::cout << std::endl;
    std::cout << "time to insert" << std::endl;
    dyn_dbg.insert(4, 3);
    std::cout << std::endl;
    
    std::cout << "time to remove" << std::endl;
    dyn_dbg.remove(2, 0);
    std::cout << std::endl;
    std::cout << "time to remove" << std::endl;
    dyn_dbg.remove(2, 1);
    std::cout << std::endl;
    std::cout << "time to remove" << std::endl;
    dyn_dbg.remove(2, 2);
    std::cout << std::endl;
    std::cout << "time to remove" << std::endl;
    dyn_dbg.remove(2, 3);
    std::cout << std::endl;

    std::cout << "STORE" << std::endl;
    //Write out to file
    dyn_dbg.store_to_file("dyn_io_test_initial");
    std::cout << std::endl;

    //Load from file 
    dyn_debruijn_graph<> dyn_dbg_from_files;
    std::cout << "reload from file after storing" << std::endl;
    dyn_dbg_from_files.load_from_file("dyn_io_test_initial");
    std::cout << std::endl;

    //Write back out to file to compare
    std::cout << "RE-STORE" << std::endl;
    dyn_dbg_from_files.store_to_file("dyn_io_test_after_load");
    std::cout << std::endl;
}

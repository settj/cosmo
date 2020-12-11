// TCLAP
#include "tclap/CmdLine.h"

// KMC
#include "kmc_api/kmc_file.h"

#include "io.hpp"

struct parameters_t {
    std::string input_filename = "";
};

parameters_t parse_arguments(int argc, char **argv) {
    parameters_t params;
    TCLAP::CmdLine cmd(banner, ' ', version);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", "Input file.", true, "", "input_file", cmd);
    cmd.parse( argc, argv );
    params.input_filename  = input_filename_arg.getValue();
    return params;
}

int main(int argc, char** argv) {
    auto params = parse_arguments(argc, argv);
    std::string file_name = params.input_filename;


    std::vector<CKMCFile*> kmer_data_bases;
    size_t num_colors;
    size_t min_union;
    size_t max_union;
    uint32_t kmc_k;
    kmc_read_header(file_name, kmc_k, min_union, max_union, num_colors, kmer_data_bases);


    std::ofstream kmer_count_file;
    kmer_count_file.open("kmers_per_color.txt");
    int color=0;
    for( CKMCFile* kmc_file : kmer_data_bases) {
        kmer_count_file << color << "," << kmc_file->KmerCount() << "\n";
        ++color;
    }
    kmer_count_file.close();

    return 0;
}

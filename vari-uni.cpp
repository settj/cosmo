//3rd Party
#include "tclap/CmdLine.h"

//Local Headers
#include "debruijn_graph_shifted.hpp"

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string output_filename = "";
  unsigned int unitig_threshold = 0;
  std::string gene_color_filename = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params) {
  TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color", ".sd_vector file.", true, "", "color_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> output_filename_arg("output", ".txt file.", true, "", "output_file", cmd);
  TCLAP::UnlabeledValueArg<unsigned int> unitig_threshold_arg("unitig_threshold", "value of unitig threshold", true, 0, "integer", cmd);
  TCLAP::UnlabeledValueArg<std::string> gene_color_filename_arg("gene_color", "file that holds colors of MEGARes genes present in the sample", true, "", "gene_color_file", cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.output_filename  = output_filename_arg.getValue();
  params.unitig_threshold  = unitig_threshold_arg.getValue();
  params.gene_color_filename  = gene_color_filename_arg.getValue();
}

void find_unitigs(debruijn_graph_shifted<>& dbg, sdsl::sd_vector<>& colors, std::vector<std::string> megares_sequences, unsigned int unitig_threshold) {
    //ith bit equals 1 when that node has been visited
    sdsl::bit_vector visited = sdsl::bit_vector(dbg.num_nodes(), 0);

    //One bit vector has 1s for start nodes, the other for end nodes
    sdsl::bit_vector start_nodes = sdsl::bit_vector(dbg.num_nodes(), 0);

    //We want significant unitigs, meaning ones that are long enough to matter. We will count the length
    //and compare against the threshold.
    unsigned int unitig_length = 1;
    unsigned int num_unitigs = 0;
    
    int num_colors = colors.size() / dbg.size();
    std::cout << num_colors << " is number of colors in graph." << std::endl;

    ofstream unitigs_fasta_during;
    unitigs_fasta_during.open("unitigs_during.fasta");

    //Loop over all the nodes in the graph looking for potential unitigs
    //for (size_t start_node = 0; start_node < dbg.num_nodes(); start_node++) {
    for (std::string kmer : megares_sequences) {
   
        auto kmer_iter = kmer.begin(); 

        while(kmer_iter != kmer.end()) {
            //Index returns a pair of size_ts that (I think?) correspond to first and last outgoing edges of a node
            boost::optional<std::pair<size_t, size_t>> start_node_opt = dbg.index(kmer_iter);
            
            if(!start_node_opt) {
                //If index does not return a valid start node, then we increment the iter and try again from the top
                kmer_iter++;
                continue;
            }

            std::pair<size_t, size_t> start_node_edge_range = *start_node_opt;
            //std::cerr << "Node's first edge: " << get<0>(start_node_edge_range) << std::endl;
            //std::cerr << "Node's second edge: " << get<1>(start_node_edge_range) << std::endl;
            //std::cerr << "First edge to node: " << dbg._edge_to_node(get<0>(start_node_edge_range)) << std::endl;
            //std::cerr << "Second edge to node: " << dbg._edge_to_node(get<1>(start_node_edge_range)) << std::endl;
            size_t start_node = dbg._edge_to_node(get<0>(start_node_edge_range));

            //If this node hasn't been visited AND it has outdegree 1 AND it has indegree >1, then it potentially starts a unitig
            //If the indegree is equal to 1, that means it's in the middle of a unitig
            if (!visited[start_node] && dbg.outdegree(start_node) == 1 && dbg.indegree(start_node) >1) {

                stringstream unitig_header;
                unitig_header << ">Unitig|" << start_node << "|";
                stringstream unitig;
                unitig << dbg.node_label(start_node).at(0);

                color_bv unitig_colors = 0;

                visited[start_node] = 1;
                unitig_length = 1;

                //Go through our alphabet and look for the edge to follow
                for (unsigned long label = 1; label < dbg.sigma + 1; label++) {
                    ssize_t edge = dbg.outgoing_edge(start_node, label);

                    //If that label is wrong, try again
                    if(edge == -1) {
                        continue;
                    }

                    for (int c = 0; c < num_colors; c++) {
                        unitig_colors |= colors[edge * num_colors + c] << c;
                    }
                    unitig << dbg.edge_label(edge);

                    //Walk along the unitig
                    ssize_t curr_node = dbg._edge_to_node(edge);
                    while(dbg.indegree(curr_node) == 1 && dbg.outdegree(curr_node) == 1) {
                        visited[curr_node] = 1;
                        unitig_length++;
                        
                        for (unsigned long label2 = 1; label2 < dbg.sigma + 1; label2++) { //Iterate over alphabet again
                            edge = dbg.outgoing_edge(curr_node, label2);
                            if(edge != -1) {
                                break;
                            }
                        }

                        for (int c = 0; c < num_colors; c++) {
                            unitig_colors |= colors[edge * num_colors + c] << c;
                        }
                        unitig << dbg.edge_label(edge).at(dbg.edge_label(edge).length() -1);

                        //Update current node and keep traversing
                        curr_node = dbg._edge_to_node(edge);
                    }

                    //When we get here, the unitig has ended. Check if it meets the threshold
                    if(unitig_length > unitig_threshold) {
                        start_nodes[start_node] = 1;
                        num_unitigs++;
                    
                        //Get color information
                        //TODO check if colors include genometrakr samples, o.w. this unitig doesn't matter
                        for (int c = 0; c < num_colors; c++) {
                            bool color_bit = unitig_colors[c];
                            if(color_bit) {
                                unitig_header << c << ", ";
                            }
                        }
                        unitigs_fasta_during << unitig_header.str() << "\n";
                        unitigs_fasta_during << unitig.str() << "\n";


                        std::cout << start_node << " is beginning of unitig with length " << unitig_length << std::endl;
                    }
                    break;
                }
            }

            //Going into the next iteration, we want to move into the next kmer of the MEGARes gene's sequence
            kmer_iter++;
        }
    }
    unitigs_fasta_during.close();
}

int main(int argc, char* argv[]) {
    parameters_t p;
    parse_arguments(argc, argv, p);

    std::cerr << "loading dbg" << std::endl;
    debruijn_graph_shifted<> dbg;
    load_from_file(dbg, p.input_filename);

    std::cerr << "loading colors" << std::endl;
    sdsl::sd_vector<> colors;
    load_from_file(colors, p.color_filename);

    cerr << "k             : " << dbg.k << endl;
    cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << dbg.num_edges() << endl;
    cerr << "colors        : " << colors.size() / dbg.size() << endl; 
    cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
    cerr << "Color size    : " << size_in_mega_bytes(colors) << " MB" << endl;


    std::cout << "NODES" << std::endl;
    for (size_t start_node = 0; start_node < dbg.num_nodes(); start_node++) {
        std::cout << start_node << ": " << dbg.node_label(start_node) << " with indegree " << dbg.indegree(start_node);
        for(int i=1; i<5; ++i) {
            ssize_t outgoing_node = dbg.outgoing(start_node, i);
            if(outgoing_node != -1) {
                std::cout << " " << i << "--" << outgoing_node;
            }
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "EDGES" << std::endl;
    for (size_t start_edge=0; start_edge < dbg.num_edges(); start_edge++) {
        std::cout << start_edge << ": " << dbg.edge_label(start_edge);
        std::cout << std::endl;
    }

    //std::cout << std::endl << "Number of preds" << std::endl;
    //for (size_t node=0; node < dbg.num_nodes(); ++node) {
        //std::cout << dbg.indegree(node) << std::endl;
    //}
    return 0; //REMOVE JUST USING TO LOOK AT NODES

    int max_megares_gene = 3824; //No color should be larger than this
    std::vector<int> megares_gene_colors;
    std::vector<int> gene_occurence_count(max_megares_gene, 0);

    //Read in colors from file
    std::ifstream gene_color_file;
    gene_color_file.open(p.gene_color_filename);
    if( gene_color_file.is_open() ) { //Make sure open was successful before going in
        std::string color_str;
        while( getline(gene_color_file, color_str) ) {
            int color_num = std::stoi(color_str);
            int curr_count = gene_occurence_count.at(color_num);
            gene_occurence_count.at(color_num) = curr_count+1;
            if(color_num == 0) {
                std::cout << curr_count << " is current count, which was updated to " << gene_occurence_count.at(color_num) << std::endl;
            }
            if (color_num < max_megares_gene && !(megares_gene_colors.end() != std::find(megares_gene_colors.begin(), megares_gene_colors.end(), color_num))) {
                megares_gene_colors.push_back(color_num);
            }
        }
    }
    gene_color_file.close();

    std::sort(megares_gene_colors.begin(), megares_gene_colors.end());
    //for(int num : megares_gene_colors) {
        //std::cerr << num << std::endl;
    //}
    for(int n=0; n<max_megares_gene;++n) {
        std::cerr << "Color " << n << " has " << gene_occurence_count.at(n) << " appearances in the color matrix." << std::endl;
    }

    //REMOVE TESTING
    return 0;


    std::cout << "Read in " << megares_gene_colors.size() << " colors.";

    std::vector<std::string> megares_filenames;
    std::ifstream megares_filenames_file;
    megares_filenames_file.open("megares_fasta_list.txt");
    if( megares_filenames_file.is_open()) {
        std::string filename;
        while( getline(megares_filenames_file, filename)) {
            megares_filenames.push_back(filename);
        }
    }
    megares_filenames_file.close();
    std::cout << "Read in " << megares_filenames.size() << " filenames.";

    //This is the goal of the whole process here: the list of kmers that we need to use to index
    //into the graph. The search for unitig will start from these nodes
    std::vector<std::string> megares_sequences;

    //For each color, look up MEGARes fasta file
    //Each color corresponds to one MEGARes gene
    std::ifstream current_gene_file;
    for(int color : megares_gene_colors) {
        std::string filename = megares_filenames.at(color);
        current_gene_file.open(filename);

        if(!current_gene_file.is_open()) {
            cerr << "Could not open file: " << filename;
            continue;
        }

        std::string sequence;
        //Get header info out of the way
        getline(current_gene_file, sequence);

        //Now get actual sequence data
        getline(current_gene_file, sequence);

        //Make kmers from the fasta file
        for(unsigned int i=0; i <= sequence.length()-dbg.k; ++i) {
            megares_sequences.push_back(sequence.substr(i, dbg.k));
        }

        current_gene_file.close();
    }

    //Pass kmers to find_unitigs
    find_unitigs(dbg, colors, megares_sequences, p.unitig_threshold);
}

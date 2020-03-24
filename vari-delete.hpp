#ifndef VARI_DELETE_H
#define VARI_DELETE_H

#include "SDIter.h"

#include "debruijn_graph_shifted.hpp"

struct parameters_t {
  std::string input_filename = "";
  std::string input_filename2 = "";    
  std::string matrix_filename = "";    
  int num_colors = 0;    
  int colors_start = 0;    
  int colors_end = 0;    
  std::string kmers_filename = "";
};

int getMilliCount();
int getMilliSpan(int nTimeStart);
void parse_arguments(int argc, char **argv, parameters_t & params);
void dump_nodes(debruijn_graph_shifted<> dbg, uint64_t * colors);
void dump_edges(debruijn_graph_shifted<> dbg, uint64_t * colors);
int dbg_delete(const debruijn_graph_shifted<> &g1, const debruijn_graph_shifted<> &g2, SDIter& color_iter, const int num_colors, const int colors_start, const int colors_end);


#endif

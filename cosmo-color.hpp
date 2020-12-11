#ifndef COSMO_COLOR_H
#define COSMO_COLOR_H

#include <vector>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "config.hpp"

template <size_t t_sigma,
          class  t_bit_vector_type,
          class  t_bv_rank_type,
          class  t_bv_select_type,
          class  t_edge_vector_type,
          class  t_symbol_type,
          class  t_label_type>
class debruijn_graph_shifted;

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string output_prefix = "";
  std::string color_mask1 = "";
  std::string color_mask2 = "";

};

int getMilliCount();
int getMilliSpan(int nTimeStart);
void parse_arguments(int argc, char **argv, parameters_t & params);
void test_symmetry(debruijn_graph_shifted<> dbg);
void dump_nodes(debruijn_graph_shifted<> dbg, uint64_t * colors);
void dump_edges(debruijn_graph_shifted<> dbg, uint64_t * colors);
void find_bubbles(debruijn_graph_shifted<> dbg, sdsl::rrr_vector<63> &colors, uint64_t color_mask1, uint64_t color_mask2);













#endif

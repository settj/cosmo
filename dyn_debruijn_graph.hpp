#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>

#include "debruijn_graph_shifted.hpp"
#include "sort.hpp"

template <size_t t_sigma            = 4, // default: DNA, TODO: change to 0 and make dynamic
          class  t_bit_vector_type  = sdsl::sd_vector<>,
          class  t_bv_rank_type     = typename t_bit_vector_type::rank_0_type, // We invert the bits so it is 0 instead
          class  t_bv_select_type   = typename t_bit_vector_type::select_0_type,
          class  t_edge_vector_type = sdsl::wt_huff<sdsl::rrr_vector<63>>,
          class  t_symbol_type      = typename t_edge_vector_type::value_type,
          class  t_label_type       = std::string> // can define basic_string<t_symbol_type>, but need to use correct print func
class dyn_debruijn_graph:
    public debruijn_graph_shifted<t_sigma, t_bit_vector_type, t_bv_rank_type, t_bv_select_type, t_edge_vector_type, t_symbol_type, t_label_type>
{
    const static size_t sigma = 4;
    typedef typename t_bit_vector_type::size_type size_type;
    typedef t_label_type  label_type;
    typedef t_symbol_type  symbol_type;
    typedef size_t edge_type;
    typedef std::pair<edge_type, edge_type> node_type;

    typedef debruijn_graph_shifted<> dbg_t;

private:
    //debruijn_graph_shifted<t_sigma, t_bit_vector_type, t_bv_rank_type, t_bv_select_type, t_edge_vector_type, t_symbol_type, t_label_type> dbg;
    boost::hash<std::string> bitstring_hash;

    boost::dynamic_bitset<> node_bits;
    boost::unordered_map<std::string, size_t> added_edges;
    boost::unordered_map<std::string, size_t> deleted_edges;
    boost::dynamic_bitset<> label;
    //TODO add color as satellite data?
    std::unordered_map<size_t, ssize_t> added_edges_node_rank;

    size_t static_graph_size;
    float resize_fraction = 0.2; //arbitrary


    //Just used for building
    struct parameters_t {
        std::string input_filename = "";
        std::string output_prefix  = "";
        std::string output_base    = "";
        size_t k = 0;
        size_t m = 0;
        bool swap = false;
        bool variable_order = false;
        bool shift_dummies  = false;
    };

    //Helper function to go from character string to string of "bits" 
    //Both strings are ascii encoded
    std::string string_to_bitstring(std::string str) {
        std::cout << str << " being converted to... ";
        std::string bitset_str(str.size()*2, '0');
        for(unsigned int i=0; i<str.size(); ++i) {
            switch(str[i]) {
                case 'A':
                    bitset_str[2*i]    = '0';
                    bitset_str[2*i +1] = '0';
                    break;
                case 'C':
                    bitset_str[2*i]    = '0';
                    bitset_str[2*i +1] = '1';
                    break;
                case 'G':
                    bitset_str[2*i]    = '1';
                    bitset_str[2*i +1] = '0';
                    break;
                case 'T':
                    bitset_str[2*i]    = '1';
                    bitset_str[2*i +1] = '1';
                    break;
                default:
                    std::cout << "ERROR CONVERTING STRING TO BITSTRING" << std::endl;
            }
        }
        std::cout << bitset_str << std::endl;
        return bitset_str;
    }

public:
    dyn_debruijn_graph() {}
    dyn_debruijn_graph(size_t in_k, const t_bit_vector_type & node_flags, const t_edge_vector_type & edges,
                           const std::array<size_t, 1+sigma>& symbol_ends, const label_type& alphabet)
    : debruijn_graph_shifted<t_sigma, t_bit_vector_type, t_bv_rank_type, t_bv_select_type, t_edge_vector_type, t_symbol_type, t_label_type>(in_k, node_flags, edges, symbol_ends, alphabet),
      label(in_k*2) {
        node_bits.resize(this->m_num_nodes, false);
        static_graph_size = this->num_edges() * 5; //Overestimate of Vari's bits/kmer
    }
 
    /**
     * L
     * @param[in] v     Node index
     *
     */
    std::map<symbol_type, size_t> list(size_t v) {
        std::map<symbol_type, size_t> result;

        //If node is in original graph, start listing the edges leaving the node from original graph
        if(node_bits.size() > v) {
            node_type n = this->get_node(v);
            size_t first_edge_idx = std::get<0>(n);
            size_t last_edge_idx  =  std::get<1>(n);
            for(size_t i=first_edge_idx; i<=last_edge_idx; ++i) {
                symbol_type out_edge_idx = _strip_edge_flag(this->m_edges[i]);

                //Don't include outgoing dummy edges? TODO make sure this is valid
                if(out_edge_idx > 0) {
                    result[out_edge_idx-1] = i;
                }
            }

            //If the node's outgiong edges have been modified -- check buffers!
            if(node_bits[v] == 1) {
                boost::dynamic_bitset<> edge_label = label << 2;

                for (int i=0; i < 4; ++i) { //Iterate over alphabet

                    //Form bitstring of the potential edge
                    std::string bitset_str;
                    to_string((edge_label | boost::dynamic_bitset<>(2, i)), bitset_str);

                    if ( deleted_edges.find(bitset_str) != deleted_edges.end()) {
                        result.erase(i);
                    }
                    if ( added_edges.find(bitset_str) != added_edges.end()) {
                        //add this character to result
                        result[i] = this->num_edges() + 1; //TODO make sure that this lines up with symbol type
                    }
                    //Clear edge label for next iteration
                    edge_label &= std::bitset<2>(0);
                }
            }
        }
        else { //If the node is not in the original graph, only have to check added edges

            /* Don't have to check if the node is modified, we know if we got here without the node being in
               the graph, then it had to be formed from added edges */
            boost::dynamic_bitset<> edge_label = label << 2;

            for (int i=0; i < 4; ++i) { //Iterate over the alphabet
                //Form bitstring of the potential edge
                std::string bitset_str;
                to_string((edge_label | boost::dynamic_bitset<>(2, i)), bitset_str);

                if ( added_edges.find(bitset_str) != added_edges.end()) {
                    //add this character to result
                    result[i] = node_bits.size() + 1; //TODO make sure that this lines up with symbol type
                }

                //Clear edge label for next iteration
                edge_label &= std::bitset<2>(0);
            }
        }
        return result;
    }

    /**
     * DESCRIBE, mention neg ret val
     * @param[in] v     Node index
     * @param[in] c     Outgoing symbol from alphabet
     *
     */
    ssize_t forward(size_t v, symbol_type c) {
        //If node is in original graph and its edges not modified
        if(node_bits.size() > v && node_bits[v] == 0) {
            label = label << 2;
            label |= boost::dynamic_bitset<>(2, c); //TODO make sure that constructor is ok
            return this->outgoing(v,c);
        } else { //Either node is not in original graph, or it is in the graph but its outgoing edges have been modified
            size_t out_edge = list(v).at(c); //TODO wrap in try-catch to see if it fails
           
            //If edge is in orginal graph, just use normal navigation
            if(this->num_edges() > out_edge) {
                label = label << 2;
                label |= boost::dynamic_bitset<>(2, c); //TODO make sure that constructor is ok
                return this->outgoing(v,c);
            } else { //Edge listed for c is not in original graph, get node rank from satellite data
                label = label << 2;
                label |= boost::dynamic_bitset<>(2, c); //TODO make sure that constructor is ok
                return added_edges_node_rank[out_edge];
            }
        }
    }

    /**
     * DESCRIBE, mention neg ret val
     * @param[in] v     Node index
     * @param[in] c     Outgoing symbol from alphabet
     *
     */
    //TODO make return size_t??
    ssize_t search(std::string target_kmer) {

        auto kmer_itr = target_kmer.begin();
        auto node_ptr = this->index(kmer_itr);
        node_type node;
        if(node_ptr) node = *node_ptr;
        std::cout << "INDEXING 0-30 OF KMER: " << std::get<0>(node) << ", " << std::get<1>(node) << std::endl;

        ++kmer_itr;
        node_ptr = this->index(kmer_itr);
        if(node_ptr) node = *node_ptr;
        std::cout << "INDEXING 1-31 OF KMER: " << std::get<0>(node) << ", " << std::get<1>(node) << std::endl;


        //Check to see if any edges leading to the target node have been added
        //for each sigma in Sigma {
            //if node_label |= std::bitst<2>(sigma) is in added_edges
                //update label
                //return the satellite data of the added_edges entry
        //}
        int target_idx = 0;
        auto c = target_kmer[target_idx];
        symbol_type first_symbol = this->_encode_symbol(c);
        // Range is from first edge of first, to last edge of last
        size_t start = this->_symbol_start(first_symbol);
        size_t end   = this->m_symbol_ends[first_symbol]-1;
        size_t first = 0, last = 0; //Initial values don't matter, but this removes warning

        bool is_in_graph = true;

        // find c-labeled pred edge
        // if outside of range, find c- labeled pred edge
        std::cout << "KMER: ";
        std::cout << c;
        for (size_t i = 0; i < this->k - 2; i++) {
            ++target_idx;
            c = target_kmer[target_idx];
            std::cout << c;
            symbol_type x = this->_encode_symbol(c);
            // update range; Within current range, find first and last occurence of c or c-
            // first -> succ(x, first)
            for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
                first = this->m_edges.select((this->m_edges.rank(start, y)) + 1, y);
                if (start <= first && first <= end) break;
            }
            if (!(start <= first && first <= end)) {
                 //Need to look in maps
                 is_in_graph = false;
                 break; 
            }
            // last -> pred(x, last)
            if (start == end) {
                last = first;
            } else {
                for (uint8_t y=x<<1; y<(x<<1)+1; y++) {
                    last = this->m_edges.select((this->m_edges.rank(end + 1, y)), y);
                    if (start <= last && last <= end) break;
                }
            }
            assert(start <= last && last <= end);
            // Follow each edge forward
            start = this->_forward(first, x);
            end   = this->_forward(last, x);
            end   = this->_last_edge_of_node(this->_edge_to_node(end));
        }
        ++target_idx;
        c = target_kmer[target_idx];
        std::cout << c;
        std::cout << std::endl;

        //If node was in graph, get kmer/edge
        if(is_in_graph) {
            std::cout << "NODE WAS IN GRAPH" << std::endl;
            symbol_type x = this->_encode_symbol(c);
            for(size_t i = start; i <= end; ++i) {
                std::cout << "CHECKING EDGE WITH LABEL: " << this->edge_label(i) << std::endl;
                std::cout << "c: " << c << ", x: " << x << ", stripped edge character: " << (this->m_edges[i] >> 1) << std::endl;
                if(x == (this->m_edges[i] >> 1)) {
                    return this->_edge_to_node(end);
                }
            }
            //If we get here, the kmer we want was not in the graph! Only the k-1mer
            is_in_graph = false;
        }

        if(!is_in_graph) {
            std::cout << "EDGE WAS NOT IN GRAPH, CHECK HASH" << std::endl;
            std::string target_kmer_bitstring = string_to_bitstring(target_kmer);
            if ( added_edges.find(target_kmer_bitstring) != added_edges.end()) {
                return added_edges[target_kmer_bitstring];
            }
            else if ( deleted_edges.find(target_kmer_bitstring) != deleted_edges.end()) {
                return deleted_edges[target_kmer_bitstring];
            }
        }
        std::cout << "EDGE NOT IN GRAPH OR HASH" << std::endl;
        return -1;
    }

    /**
     * DESCRIBE
     * @param[in] kmer  String
     *
     */
    void add(std::string kmer) {
        //Look for node (k-1mer) in graph
        boost::optional<node_type> node = this->index(kmer.begin());

        if(node) {
            //Can convert to unsigned becausee existence of node guarantees existence of edge
            size_t node_id = this->_node_to_edge((size_t)std::get<0>(*node));
            this->insert(node_id, this->_encode_symbol(kmer[kmer.size()-1]));
        }
        else {
            std::string bitset_str = this->string_to_bitstring(kmer);
            added_edges.emplace(bitset_str, bitstring_hash(bitset_str));
        }
    }

    /**
     * DESCRIBE, mention neg ret val
     * @param[in] v     Node index
     * @param[in] c     Outgoing symbol from alphabet
     *
     */
    void insert(size_t v, symbol_type c) {
        //Get bitstring to hash with
        boost::dynamic_bitset<> tmp(2, c);
        boost::dynamic_bitset<> tmp_label(label.size(), 0);
        tmp_label = label << 2;
        tmp_label = tmp_label | tmp;
        std::string bitset_str(label.size()+2, '0');
        boost::to_string(tmp_label, bitset_str);

        //Now, do the insertion
        
        //If edge not in graph
        if(this->outgoing_edge(v, c) == -1) {
            added_edges.emplace(bitset_str, bitstring_hash(bitset_str));

            if(too_much_space()) {
                build_new_dbg();
            }
        }
        else { //Edge in graph
            //If edge was deleted, un-delete it
            if ( deleted_edges.find(bitset_str) != deleted_edges.end()) {
                deleted_edges.erase(bitset_str);
            }
        }
    }

    /**
     * DESCRIBE, mention neg ret val
     * @param[in] v     Node index
     * @param[in] c     Outgoing symbol from alphabet
     *
     */
    void remove(size_t v, symbol_type c) {
        std::string bitset_str;
        to_string(((label << 2) | boost::dynamic_bitset<>(2, c)), bitset_str);

        //If edge not in graph
        if(this->outgoing_edge(v, c) == -1) {
            //If edge was added, un-add it
            if ( added_edges.find(bitset_str) != added_edges.end()) {
                added_edges.erase(bitset_str);
            }
        }
        else { //Edge in graph
            std::cout << bitset_str << std::endl;
            deleted_edges.emplace(bitset_str, bitstring_hash(bitset_str));
            std::cout << deleted_edges.size() << std::endl;

            if(too_much_space()) {
                build_new_dbg();
            }
        }

    }
    
    /**
     * DESCRIBE
     * @param[in] 
     * @param[in]
     *
     */
    bool too_much_space() {
        //TODO standardize UNITS OF SPACE (bytes?)
        //Space in buffers is number of total buffer entries times size of each key (bitstring) and each value (size_t)
        //Space in satellite data is number of entires times size of each key and each value (both size_t)
        size_t space = (added_edges.size() + deleted_edges.size())*(this->k*2 + sizeof(size_t))
                       + added_edges_node_rank.size()*2*sizeof(size_t);
        return space > (resize_fraction * static_graph_size);
    }

    /**
     * DESCRIBE
     * @param[in] 
     * @param[in]
     *
     */
    unsigned long long bitstring_to_int64(std::string bitstring) {
        //TODO use bitset
        return 0;
    }

    /**
     * DESCRIBE
     * @param[in] 
     * @param[in]
     *
     */
    void build_new_dbg() {
        typedef dbg_builder<dbg_t, kmer_t> builder_t;
        parameters_t params;
        builder_t builder(params);

        auto map_itr = deleted_edges.begin();
        while(map_itr != deleted_edges.end()) {
            //builder.push(bitstring_to_int64(*map_itr),0); //TODO color
        }



        map_itr = added_edges.begin();
        while(map_itr != added_edges.end()) {
            //builder.push(bitstring_to_int64(*map_itr),0); //TODO color
        }

        size_t num_colors = 0; //TODO pull this from satellite data
        size_t num_set = 0;
        size_t edge_idx = 0;
        auto added_dbg = builder.build([&](auto x) { // Pre merge
                // color_bv = bit_vector(x * num_colors);
            },[&](auto x) { // Merge visitor
                //auto kmer = x.edge;
                //cerr << kmer_to_string(kmer,k) << endl;
                //auto color = boost::get<0>(x.payload);
                //auto color = x.payload.get<0>();
                //serialize_color_bv(cfs, color);
                for (size_t color_idx = 0; color_idx < num_colors; color_idx++) {
                    // TODO: try complemented bits in color-major form in a *sd_vector* for large data set.
                    //color_bv[color_idx * num_colors + edge_idx] = !color[color_idx];
                    //color_bv[edge_idx * num_colors + color_idx] = color[color_idx];
                    //num_set += color[color_idx];
                }
                edge_idx++;
            }
        );
    }


    /**
     * DESCRIBE
     * @param[in] 
     * @param[in]
     *
     */
    void load_from_file(std::string basename) {
        //Load boss dbg
        sdsl::load_from_file((*this), basename+".dbg");
        cerr << "k             : " << this->k << endl;
        cerr << "num_nodes()   : " << this->num_nodes() << endl;
        cerr << "num_edges()   : " << this->num_edges() << endl;
        //cerr << "Total size    : " << size_in_mega_bytes(*this) << " MB" << endl;
        //cerr << "Bits per edge : " << bits_per_element(*this) << " Bits" << endl <<endl;
        label.resize(this->k);
        std::cout << label.size() << " is bitset size" << std::endl;

        //Load H+
        std::cout << "LOADING H+" << std::endl;
        std::ifstream hplus_ifs(basename + "_Hplus.txt");
        std::string line;
        typename decltype(added_edges)::iterator hplus_itr = added_edges.begin();
        while(std::getline(hplus_ifs, line)) {
            std::cout << line << std::endl;
            added_edges.emplace(line, bitstring_hash(line));
        }

        //Load H-
        std::cout << "LOADING H-" << std::endl;
        std::ifstream hminus_ifs(basename + "_Hminus.txt");
        typename decltype(added_edges)::iterator hminus_itr = deleted_edges.begin();
        while(std::getline(hminus_ifs, line)) {
            std::cout << line << std::endl;
            deleted_edges.emplace(line, bitstring_hash(line));
        }

    }


    /**
     * DESCRIBE
     * @param[in] 
     * @param[in]
     *
     */
    void store_to_file(std::string basename) {
        //Write out base graph
        sdsl::store_to_file(static_cast<debruijn_graph_shifted<>>(*this), basename + ".dbg");

        //Write out H+
        std::ofstream hplus_ofs(basename + "_Hplus.txt");
        typename decltype(added_edges)::iterator hplus_itr = added_edges.begin();
        while(hplus_itr != added_edges.end()) {
            std::pair<std::string, size_t> entry = *hplus_itr;
            hplus_ofs << std::get<0>(entry);
            hplus_ofs << "\n";
            hplus_itr++;
        }

        //Write out H-
        std::ofstream hminus_ofs(basename + "_Hminus.txt");
        typename decltype(added_edges)::iterator hminus_itr = deleted_edges.begin();
        while(hminus_itr != deleted_edges.end()) {
            std::pair<std::string, size_t> entry = *hminus_itr;
            hminus_ofs << std::get<0>(entry);
            hminus_ofs << "\n";
            hminus_itr++;
        }
    }

};

#include <iostream>
#include <string>
#include <sys/timeb.h>
#include <vector>

#include "tclap/CmdLine.h" //parse args

#include "sort.hpp"

struct test_vd_parameters_t {
    std::string target_filename = "";
    std::string input_filename = "";    
};


//shit i needed from vari-delete.cpp because the build process is so fucked i can't link vari-delete.

int getMilliCount(){
    timeb tb;
    ftime(&tb);
    int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
    return nCount;
}


int getMilliSpan(int nTimeStart){
    int nSpan = getMilliCount() - nTimeStart;
    if(nSpan < 0)
        nSpan += 0x100000 * 1000;
    return nSpan;
}


// Returns the number of 0's in a run beginning at set_start
uint64_t length(const uint64_t set_start, const std::vector<bool> &sets)
{
    uint64_t i = 0;
    for ( i= set_start; sets[i] != true; ++i)
        ;
    return i - set_start;
}

// Returns the next set_start given some set_start
uint64_t advance(const uint64_t set_start, const std::vector<bool> &sets)
{
    return length(set_start, sets) + 1;
}

// Returns true if the given set_start is in the last set
bool last(const uint64_t set_start, const std::vector<bool> &sets)
{

    return set_start + length(set_start, sets)  == sets.size() - 1;
}

// Returns the number of sets in positions greater or equal to start
uint64_t num_sets_ge(const uint64_t set_start, const std::vector<bool> &sets)
{
    assert (set_start <= sets.size() - 1);
    uint64_t count = 0;
    for (uint64_t i = set_start; !last(i, sets); i += advance(i, sets)) {
        count++;
    }
    return count + 1;

}

uint64_t num_sets(const std::vector<bool> &sets)
{
    return num_sets_ge(0, sets);
}


uint64_t validate(const std::vector<bool> &sets)
{
    uint64_t sum = 0;
    uint64_t num = num_sets(sets);
    uint64_t i = 0;
    uint64_t pos = 0;
        
    for (; i < num - 1; ++i) {
        sum += length(pos, sets);
        pos += advance(pos, sets);
    }
    assert(last(pos, sets));
    return sum + length(pos, sets);
}
uint64_t maxsize(const std::vector<bool> &sets)
{
    uint64_t sum = 0;
    uint64_t num = num_sets(sets);
    uint64_t i = 0;
    uint64_t pos = 0;
    uint ms = 0;
    for (; i < num - 1; ++i) {
        sum += length(pos, sets);
        if (length(pos, sets) > 1) {
            ms = length(pos, sets);
            std::cout << "max set sized bumpbed to " << ms << " at set " << i << " covering interval ending at element " << sum << std::endl;
            
        }
        pos += advance(pos, sets);
    }
    assert(last(pos, sets));
    if (length(pos, sets) > 1) {
        ms = length(pos, sets);
        std::cout << "max set sized bumpbed to " << ms << " at final set covering interval ending at element " << sum << std::endl;        
    }
    return ms;
//    return sum + length(pos, sets);
}

void subdivide2(const std::vector<unsigned char> &g1_col, const uint64_t g1_ptr, const uint64_t g1_num,
              const std::vector<unsigned char> &g2_col, const uint64_t g2_ptr, const uint64_t g2_num,
              /*const std::vector<unsigned char> &test_g2_col, const uint64_t test_g2_ptr, const uint64_t test_g2_num,*/
              /*const int colno,*/ std::vector<bool> &g1_out_set, std::vector<bool> &g2_out_set, /*std::vector<bool> &test_g2_out_set,*/ int &active_alpha_size)
{
    uint64_t g1_local_ptr = g1_ptr;
    uint64_t g2_local_ptr = g2_ptr;
    char chars_seen = 0;

    while (g1_local_ptr < g1_ptr + g1_num || g2_local_ptr < g2_ptr + g2_num) {
        
        // assumtion 1: we just finished a run or we're at the beginning

        // if only col1 has characters:
        //   consume characters (adding zeros) while the run continues, then add one to both out sets, then go back to assumption 1
        
        if (g2_local_ptr == g2_ptr + g2_num) {
            unsigned char run_char = g1_col[g1_local_ptr];
            while (g1_local_ptr <  g1_ptr + g1_num && run_char == g1_col[g1_local_ptr]) {
                g1_local_ptr++;
                g1_out_set.push_back(false);
            }
            chars_seen++;
        }

        // if only col2 has chracters:
        //   consume characters (adding zeros) while the run continues, then add one to both out sets, then go back to assumption 1
        else if (g1_local_ptr == g1_ptr + g1_num) {
            unsigned char run_char = g2_col[g2_local_ptr];
            while (g2_local_ptr < g2_ptr + g2_num && run_char == g2_col[g2_local_ptr]) {
                g2_local_ptr++;
                g2_out_set.push_back(false);
            }
        }


        // if col1 and col2 have characters:
        //   if they are the same character:
        //      process all of col1, process all of col2, add ones to both bit sets
        else  if (g1_col[g1_local_ptr] ==  g2_col[g2_local_ptr]) {
            
            unsigned char run_char = g1_col[g1_local_ptr];
            while (g1_local_ptr <  g1_ptr + g1_num && run_char == g1_col[g1_local_ptr]) {
                g1_local_ptr++;
                g1_out_set.push_back(false);
            }
            while (g2_local_ptr < g2_ptr + g2_num && run_char == g2_col[g2_local_ptr]) {
                g2_local_ptr++;
                g2_out_set.push_back(false);
            }
        }

        
        //   else:
        //      process all characters from lexicographically smaller one as if the other one was empty, add one to both out sets
        else if (g1_col[g1_local_ptr] < g2_col[g2_local_ptr]) {
            unsigned char run_char = g1_col[g1_local_ptr];
            while (g1_local_ptr <  g1_ptr + g1_num && run_char == g1_col[g1_local_ptr]) {
                g1_local_ptr++;
                g1_out_set.push_back(false);
            }
            chars_seen++;
        } else {
            unsigned char run_char = g2_col[g2_local_ptr];
            while (g2_local_ptr < g2_ptr + g2_num && run_char == g2_col[g2_local_ptr]) {
                g2_local_ptr++;
                g2_out_set.push_back(false);
            }
        }
        //chars_seen++;
        g1_out_set.push_back(true);
        g2_out_set.push_back(true);
     

    
        // look at head of each queue
        // if either contains cur_char, read from queue
    
    }

    active_alpha_size = chars_seen;

    // if c1 is empty
}

int refine_sets(const std::vector<unsigned char> &g1_col, const std::vector<unsigned char> &g2_col, /*const std::vector<unsigned char> &test_g2_col,*/
                const std::vector<bool>& g1_sets, const std::vector<bool> &g2_sets, /*const std::vector<bool> &test_g2_sets,*/
                const int colno,
                std::vector<bool> &g1_out_set, std::vector<bool> &g2_out_set, /*std::vector<bool> &test_g2_out_set,*/ std::vector<bool> &Lval)
{

    uint64_t g1_ptr = 0;
    uint64_t g2_ptr = 0;    
    //uint64_t test_g2_ptr = 0;    
    
    uint64_t g1_set_start = 0;
    uint64_t g2_set_start = 0;
    //uint64_t test_g2_set_start = 0;
    //std::cout << "subdividing " << num_sets(g1_sets) << " and " << num_sets(g2_sets) << " sets" << std::endl;
    do {
        uint64_t g1_num = length(g1_set_start, g1_sets);
        uint64_t g2_num = length(g2_set_start, g2_sets);
        //uint64_t test_g2_num = length(test_g2_set_start, test_g2_sets);
        int active_alpha_size = 0;
//        int active_alpha_size2 = 0;        
        //std::cout << "subdividing ranges (col " << colno <<" ) " << std::endl << "\t" << g1_ptr << ":+" << g1_num << " = " ;
        //dump_range(g1_ptr, g1_ptr+ g1_num, g1_col);
        //std::cout << std::endl << "\t" << g2_ptr << ":+" << g2_num;
        //dump_range(g2_ptr, g2_ptr+ g2_num, g2_col);
        //std::cout <<std::endl;
        //uint64_t g1_out_set_initsize = g1_out_set.size();
        //uint64_t g2_out_set_initsize = g2_out_set.size();
//         std::vector<bool> g1t1,g2t1,g1t2,g2t2;
//         subdivide(g1_col, g1_ptr, g1_num,
//                   g2_col, g2_ptr, g2_num,
// //                  colno,
//                   g1t1,//g1_out_set,
//                   g2t1,//g2_out_set,
//                   active_alpha_size2);


          subdivide2(g1_col, g1_ptr, g1_num,
                  g2_col, g2_ptr, g2_num,
                  /*test_g2_col, test_g2_ptr, test_g2_num,*/
//                  colno,
                     g1_out_set,
                     g2_out_set,
                     //test_g2_out_set,
                  active_alpha_size);

          // if (g1t1 != g1t2 || g2t1 != g2t2 || active_alpha_size != active_alpha_size2) {
          //     dumprange(g1_col, g1_ptr, g1_num);
          //     dumprange(g2_col, g2_ptr, g2_num);
          //     dumpboolrange(g1t2,0,g1t2.size());
          //     dumpboolrange(g2t2,0,g2t2.size());
          //     exit(1);
          // }

          // g1_out_set.insert(g1_out_set.end(), g1t2.begin(), g1t2.end());
          // g2_out_set.insert(g2_out_set.end(), g2t2.begin(), g2t2.end());
        //std::cout << "\talpha size: " << active_alpha_size << std::endl;
        //std::cout << "\tg1_out_set size: " << g1_out_set.size() << ", g2_out_set size: " << g2_out_set.size() << std::endl;
        //std::cout << "validations  " << validate(g1_out_set) << " " << validate(g2_out_set) << std::endl;
        //std::cout << "\tg1_out_set additional sets: " << ((g1_out_set_initsize > g1_out_set.size()) ? (num_sets_ge(g1_out_set_initsize + 1, g1_out_set)) : 0)
        //         << ", g2_out_set additional sets: " << ((g2_out_set_initsize > g2_out_set.size()) ? num_sets_ge(g2_out_set_initsize + 1, g2_out_set) : 0) << std::endl;        
        //std::cout << "alphabet size: " << active_alpha_size << std::endl;
        if (colno == 0 && active_alpha_size > 0) {
            Lval.push_back(0);
            Lval.insert(Lval.end(), active_alpha_size - 1, 1);
            std::cout << "LVal stuff" << std::endl;
            for(unsigned int i=0; i<Lval.size(); ++i) {
                std::cout << Lval.at(i) << std::endl;
            }
        }
        
        g1_ptr += g1_num;
        g2_ptr += g2_num;
        //test_g2_ptr += test_g2_num;
        assert(g1_ptr <= g1_col.size());
        assert(g2_ptr <= g2_col.size());
            
        if (last(g1_set_start, g1_sets) || last(g2_set_start, g2_sets)) break;
        //if (last(g1_set_start, g1_sets) || last(g2_set_start, g2_sets) || last(test_g2_set_start, test_g2_sets)) break;
        
        g1_set_start += advance(g1_set_start, g1_sets);
        g2_set_start += advance(g2_set_start, g2_sets);
        //test_g2_set_start += advance(test_g2_set_start, test_g2_sets);
    } while (true);

    return 0;
}

int get_column(const debruijn_graph_shifted<> &g, const int col_num, std::vector<unsigned char> &g_col)
{
    assert (g_col.size() == 0);
    for (uint64_t i = 0; i < g.num_edges(); ++i) {
        if (col_num == 0 ) {
            g_col.push_back(g._map_symbol(g._strip_edge_flag(g.m_edges[i])));
        } else {
            std::string label = g.edge_label(i); 
            const char *lab = label.c_str();
            if (g.k - col_num - 1>= strlen(lab)) {
                g_col.push_back('$');
            } else {
                char c = lab[g.k - col_num - 1];
                g_col.push_back(c ? c : '$');
            }
        }
    }
    return 0;
}

class Flags {

public:
    Flags(const std::vector<bool> &g1_flagsets, const std::vector<bool> &g2_flagsets);
    const std::vector<bool> &g1_flagsets;
    const std::vector<bool> &g2_flagsets;
    bool seen(char nt);
    void add(char nt);
    void adv_g1();
    void adv_g2();
private:
    std::set<char> newflags;
    void check_flagsets();
    uint64_t g1_ptr = 0;
    uint64_t g2_ptr = 0;
    uint64_t g1_base = 0;
    uint64_t g2_base = 0;
};

Flags::Flags(const std::vector<bool> &g1_flagsets_a, const std::vector<bool> &g2_flagsets_a) : g1_flagsets(g1_flagsets_a), g2_flagsets(g2_flagsets_a)
{
}

bool Flags::seen(char c)
{
    return newflags.count(c);
}


void Flags::add(char c)
{
    newflags.insert(c);
}


void Flags::check_flagsets()
{
    if (g1_flagsets[g1_base + g1_ptr] && g2_flagsets[g2_base + g2_ptr]) {
        g1_base += g1_ptr + 1;
        g2_base += g2_ptr + 1;
        g1_ptr = 0;
        g2_ptr = 0;
        newflags.clear();
        //std::cout << "g1 and g2 base for Flags are now: " << g1_base << " and " << g2_base << " out of " << g1_flagsets.size() << " and " << g2_flagsets.size() << std::endl;
    }
}
    
void Flags::adv_g1()
{
    g1_ptr += 1;
    check_flagsets();
}

void Flags::adv_g2()
{
    g2_ptr += 1;
    check_flagsets();
}



void fill_Lcol(const debruijn_graph_shifted<> &g, std::vector<bool> &Lcol)
{
    size_t p = 1;
    size_t Lindex = 0;
    do {
        Lindex = g.m_node_select(p);
        Lcol[Lindex] = true;
        p += 1;
    } while (Lindex < Lcol.size() - 16 /*FIXME: this is assumed to be realted to other truncation*/);
}

char combine(char symbol, bool flag)
{
    char encoded=0;
    switch(symbol) {
    case '$' : encoded = 0;
        break;
    case 'A' : encoded = 1;
        break;
    case 'C' : encoded = 2;
        break;
    case 'G': encoded = 3;
        break;
    case 'T': encoded = 4;
        break;
    }
    char ret = (encoded << 1 ) | flag;
    //std::cout << " (" << symbol << ":" << flag << ")=>" << (int)ret << " " << std::endl;
    return ret;
                    
}

bool curr_label_precedes_node_label_2(std::string& dummy_label, int dummy_label_idx, std::string node_label) {

    std::cout << "comparing for preceding with " << node_label << std::endl;
    int node_label_idx = node_label.length()-1;
    while(dummy_label_idx > 0) {
        if(dummy_label.at(dummy_label_idx) < node_label[node_label_idx]) {
            return true;
        }

        --dummy_label_idx;
        --node_label_idx;
    }

    return false;
}

bool curr_label_precedes_node_label(std::string& dummy_label, int dummy_label_idx, std::string node_label) {

    std::cout << "comparing for preceding with " << node_label << std::endl;
    if(dummy_label_idx < 2) { return true; }

    int node_label_idx = node_label.length()-1;
    while(dummy_label_idx > 0) {
        if(dummy_label.at(dummy_label_idx) < node_label[node_label_idx]) {
            return true;
        }

        --dummy_label_idx;
        --node_label_idx;
    }

    return false;
}

void insert_incoming_dummy_chain(std::vector<char>& edges, std::vector<bool>& L, std::map<char, uint64_t>& ntcounts, std::vector<ssize_t>& edge_mapping, size_t edge_idx, const debruijn_graph_shifted<>& g1, std::unordered_multimap<ssize_t, std::pair<std::string, char> >& dummies_to_insert) {
    int k = g1.k;
    char ebwt_char;
    //size_t edge_idx = g1._node_to_edge(i);
    std::cout << g1.edge_label(edge_idx) << " is full edge label. " << std::endl;
    std::string full_edge_label = g1.edge_label(edge_idx);

    ssize_t mapped_edge_idx;
    for(int i=1; i < k-1; ++i) {
        ebwt_char = g1._map_symbol(g1._symbol_access(edge_idx));
        edge_idx = g1._backward(edge_idx);
        ntcounts[g1._map_symbol(g1._symbol_access(edge_idx))] += 1; //Go ahead and update counts before edge_idx is modified in map loop
        std::cout << g1.edge_label(edge_idx) << " is " << i << "th edge label" << std::endl;

        switch(ebwt_char) {
            case '$' : ebwt_char = 0; //std::cout << "0" << std::endl;
                break;
            case 'A' : ebwt_char = 1; //std::cout << "1" << std::endl;
                break;
            case 'C' : ebwt_char = 2; //std::cout << "2" << std::endl;
                break;
            case 'G': ebwt_char = 3; //std::cout << "3" << std::endl;
                break;
            case 'T': ebwt_char = 4; //std::cout << "4" << std::endl;
                break;
        }

        //do the correct search in order to identify the edge in the original graph before which the dummy
        //should be inserted, then do mapping
        char sorted_col_char = g1._map_symbol(g1._symbol_access(edge_idx));
        std::cout << "sorted col char: " << sorted_col_char << std::endl;
        size_t lower_node_bound = g1._symbol_start(sorted_col_char);
        size_t upper_node_bound = g1._edge_to_node(edge_idx);
        size_t chosen_node=0;

        //TODO make sure chosen node will be found no matter what
        for(size_t j=lower_node_bound; j<=upper_node_bound; ++j) {
            if(curr_label_precedes_node_label(full_edge_label, i, g1.node_label(j))) {
                chosen_node = j;
                break;
            }
        }

        std::cout << "chosen node, chosen edge: " << chosen_node << ", " << g1._node_to_edge(chosen_node) << std::endl;
        size_t chosen_edge = g1._node_to_edge(chosen_node);
        mapped_edge_idx = edge_mapping.at(chosen_edge);
        std::cout << mapped_edge_idx << " is original mapped edge idx" << std::endl;
        while(mapped_edge_idx == -1) {
            ++chosen_edge;
            mapped_edge_idx = edge_mapping.at(chosen_edge);
            std::cout << mapped_edge_idx << " is potential actually mapped edge idx" << std::endl;
        }
        std::cout << mapped_edge_idx << " is mapped edge idx we are going to use" << std::endl;


        //Input right character from node label in edges, '0' in L (and modify next maybe?), add to ntcounts['$']
        std::string kmer = std::string(k-i-1, '$') + full_edge_label.substr(0, i);
        std::cout << kmer << " is " << (k-i-1) << "th dummy kmer" << std::endl;
        std::pair<std::string, char> dummy_kmer_char(kmer, ebwt_char);
        dummies_to_insert.insert(std::pair<ssize_t, std::pair<std::string, char> >(mapped_edge_idx, dummy_kmer_char));
        //std::cout << "incoming dummy edge at " << mapped_edge_idx << std::endl;
        //edges.insert(edges.begin() + mapped_edge_idx, (ebwt_char << 1));
        //L.insert(L.begin() + mapped_edge_idx, false);

    }

    std::cout << "edges before inserting very first dummy:" << std::endl;
    for(size_t i=0; i < edges.size(); ++i) {
        std::cout << (edges.at(i) >> 1) << std::endl;//((char)edges.at(i) & 1) << std::endl;
    }
    //First, see if node consisting of all dummies exists. If not, insert it.
    //std::cout << "inserting dummy node first!!! " << std::endl;
    ebwt_char = g1._map_symbol(g1._symbol_access(edge_idx));

    switch(ebwt_char) {
    case '$' : ebwt_char = 0; //std::cout << "0" << std::endl;
        break;
    case 'A' : ebwt_char = 1; //std::cout << "1" << std::endl;
        break;
    case 'C' : ebwt_char = 2; //std::cout << "2" << std::endl;
        break;
    case 'G': ebwt_char = 3; //std::cout << "3" << std::endl;
        break;
    case 'T': ebwt_char = 4; //std::cout << "4" << std::endl;
        break;
    }


    int offset = 0;
    bool insert = true;
    for(size_t i=0; i < ntcounts['$']; ++i) {
        if( (edges.at(i) >> 1) < ebwt_char ) {
            std::cout << (edges.at(i) >> 1) << " " << ebwt_char << " means add to offset " << std::endl;
            
            ++offset;
        }
        else if( (edges.at(i) >> 1) == ebwt_char ) {
            insert = false;
            break;
        }
    }

    if(insert) {
        //std::cout << "inserting into beginning of eges with offset: " << offset << std::endl;
        //edges.insert(edges.begin() + offset, (ebwt_char << 1));
        //L.insert(L.begin(), false);
        ntcounts['$'] += 1;
        std::string kmer = std::string(k-1, '$') + full_edge_label[k-1];
        std::pair<std::string, char> dummy_kmer_char(kmer, ebwt_char);
        dummies_to_insert.insert(std::pair<ssize_t, std::pair<std::string, char> >(offset, dummy_kmer_char));
    }

    //std::cout << "edges after inserting very first dummy:" << std::endl;
    //for(size_t i=0; i < edges.size(); ++i) {
        //std::cout << (edges.at(i) >> 1) << std::endl;//((char)edges.at(i) & 1) << std::endl;
    //}
}

void insert_outgoing_dummy(std::vector<char>& edges, std::vector<bool>& L, std::map<char, uint64_t>& ntcounts, std::vector<ssize_t>& edge_mapping, size_t edge_idx, const debruijn_graph_shifted<>& g1, std::unordered_multimap<ssize_t, std::pair<std::string, char> >& dummies_to_insert) { //, std::map<size_t) {
    //Determine where to put dummy edge in new graph based on mapping of indices of outgoing edges of node in g1
    //size_t edge_idx = g1._node_to_edge(i);
    ntcounts[g1._map_symbol(g1._symbol_access(edge_idx))] += 1; //go ahead and increment the count here before edge_idx is modified in loop

    std::string full_edge_label = g1.edge_label(edge_idx);
    //do the correct search in order to identify the edge in the original graph before which the dummy
    //should be inserted, then do mapping
    char sorted_col_char = g1._map_symbol(g1._symbol_access(edge_idx));
    std::cout << "sorted col char: " << sorted_col_char << std::endl;
    size_t lower_node_bound = g1._symbol_start(sorted_col_char);
    size_t upper_node_bound = g1._edge_to_node(edge_idx);
    size_t chosen_node=0;

    //TODO make sure chosen node will be found no matter what
    for(size_t j=lower_node_bound; j<=upper_node_bound; ++j) {
        if(curr_label_precedes_node_label_2(full_edge_label, g1.k-2, g1.node_label(j))) {
            chosen_node = j;
            break;
        }
    }

    std::cout << "chosen node, chosen edge: " << chosen_node << ", " << g1._node_to_edge(chosen_node) << std::endl;
    size_t chosen_edge = g1._node_to_edge(chosen_node);
    ssize_t mapped_edge_idx = edge_mapping.at(chosen_edge);
    std::cout << mapped_edge_idx << " is original mapped edge idx" << std::endl;
    while(mapped_edge_idx == -1) {
        ++chosen_edge;
        mapped_edge_idx = edge_mapping.at(chosen_edge);
        std::cout << mapped_edge_idx << " is potential actually mapped edge idx" << std::endl;
    }
    std::cout << mapped_edge_idx << " is mapped edge idx we are going to use" << std::endl;


    //Input '$' in edges, '0' in L (and modify next maybe?)
    //
    //TODO actually track where edges are inserted; adding one will only work for this current test 2/5/2020
    std::string kmer = g1.edge_label(edge_idx);
    kmer = kmer.substr(0, kmer.size()-1) + "$";
    std::pair<std::string, char> dummy_kmer_and_char(kmer, 0);
    dummies_to_insert.insert(std::pair<ssize_t, std::pair<std::string, char> >(mapped_edge_idx, dummy_kmer_and_char));
    //edges.insert(edges.begin() + mapped_edge_idx+1, 0);
    //L.insert(L.begin() + mapped_edge_idx+1, false);
}

/*
 * TODO 
 */
int dbg_delete(const debruijn_graph_shifted<> &g1, const debruijn_graph_shifted<> &g2/*, kmer_vector_t::bufreader_type& kmers_reader, SDIter& color_iter, const int num_colors, const int colors_start, const int colors_end*/)
{

    //
    // planning phase (same as planning done for a merge)
    //
    
    std::vector<bool> g1_sets(g1.num_edges() + 1, false);
    g1_sets[g1_sets.size() - 1] = true;
    std::cout << "created g1_sets with " << num_sets(g1_sets) << " sets of size " << length(0, g1_sets) << std::endl;

    std::vector<bool> g2_sets(g2.num_edges() + 1, false);
    //std::vector<bool> test_g2_sets(g2.num_edges() + 1, false);
    g2_sets[g2_sets.size() - 1] = true;
    //test_g2_sets[test_g2_sets.size() - 1] = true;
    std::cout << "created g2_sets with " << num_sets(g2_sets) << " sets of size " << length(0, g2_sets) << std::endl;
    std::cout << std::endl;

    std::vector<bool> L; // node flags

    //assert(g1.k == 10); // FIXME: remove after prototyping, just want to make sure behavior matches Debby's
    

    std::vector<unsigned int> cols;  // column iteration order, counting columns from the right
    for (unsigned int i = 1; i < g1.k ; ++i) {
        cols.push_back(i);
    }
    cols.push_back(0);

    std::vector<bool> g1_flagsets;//(/*g1_sets*/);
    std::vector<bool> g2_flagsets;//(/*g2_sets*/);

    // fill the cache

    // first for m_edges
    std::vector<unsigned char> g1_edges(g1.num_edges(),0); // cache WT data
    std::vector<unsigned char> g2_edges(g2.num_edges(),0); // cache WT data
    std::cerr << "caching wavelet tree edges " << std::endl << std::flush;
    int startgettime = getMilliCount();
    g1.get_edges(g1_edges);
    int delta1 = getMilliSpan(startgettime);
    std::cerr << "got abbreviated, flagged edges in " << delta1 << " milliseconds." << std::endl << std::flush;


    int startgettime2 = getMilliCount();
    g2.get_edges(g2_edges);
    int delta2 = getMilliSpan(startgettime2);
    std::cerr << "got abbreviated, flagged edges in " << delta2 << " milliseconds." << std::endl << std::flush;    

    std::vector<unsigned char> g1_col(g1.num_edges());
    g1.get_edge_column(g1_edges, g1_col);
    std::vector<unsigned char> g2_col(g2.num_edges());
    g2.get_edge_column(g2_edges, g2_col);
    
    // then for m_node_flags
    boost::dynamic_bitset<> g1_node_flags(g1.num_edges());
    boost::dynamic_bitset<> g2_node_flags(g2.num_edges());    

    std::cerr << "caching node flags for g1" << std::endl << std::flush;
    int startgettime3 = getMilliCount();
    size_t num_set1 = g1.get_node_flags(g1_node_flags);
    int delta3 = getMilliSpan(startgettime3);
    std::cerr << "got  " << num_set1 << ", reported " <<g1_node_flags.count() << " bits  in g1_node_flags in " << delta3 << " milliseconds." << std::endl << std::flush;


    std::cerr << "caching node flags for g2" << std::endl << std::flush;
    int startgettime4 = getMilliCount();
     size_t num_set2 =    g2.get_node_flags(g2_node_flags);
    int delta4 = getMilliSpan(startgettime4);
    std::cerr << "got " << num_set2 << ", reported " << g2_node_flags.count() << " bits in g2_node_flags in " << delta4 << " milliseconds." << std::endl << std::flush;


    std::vector<size_t> g1_node_rank_cache;
    g1.init_rank_cache(g1_node_rank_cache);

    std::vector<size_t> g2_node_rank_cache;
    g2.init_rank_cache(g2_node_rank_cache);

    std::cout << std::endl << "*** Starting delete planning phase. ***" << std::endl;

    
    int colno = 0;
    int planstarttime = getMilliCount();
    for (auto col: cols) {
        std::cout << "planning phase at col: " << col << std::endl;

        // get_column // FIXME: be more careful here, maybe use col^1 from g._symbol_starts
        if (col == 0) {
            g1.get_edge_column(g1_edges, g1_col);
            g2.get_edge_column(g2_edges, g2_col);
    
                
        } else  {
            std::vector<unsigned char> h1_col(g1_col.size(),0);
            g1.get_column(g1_node_rank_cache, g1_node_flags, g1_edges, g1_col, h1_col);
            g1_col.assign(h1_col.begin(), h1_col.end());
            // g1_col.clear();
            // g1_col.insert(g1_col.begin(), h1_col.begin(), h1_col.end());
            
            std::vector<unsigned char> h2_col(g2_col.size(),0);
            g2.get_column(g2_node_rank_cache, g2_node_flags, g2_edges, g2_col, h2_col);
            g2_col.assign( h2_col.begin(), h2_col.end());
            // g2_col.clear();
            // g2_col.insert(g2_col.begin(), h2_col.begin(), h2_col.end());
        }

        //std::vector<unsigned char> g1_col; // FIXME: change 'char' type to something less static
        //std::vector<unsigned char> g2_col;

        // g1_col.clear();
        // g2_col.clear();
        // get_column(g1, col, g1_col);
        // get_column(g2, col, g2_col);

            // std::stringstream fname;
            // fname << "ecolipre1.dbg";
            // fname << ".new";
            // fname << col;
            // std::string s = fname.str();
            // std::cout << "writing " << s << std::endl;
            // ofstream f(s.c_str());
            // for (auto c: g1_col) {
            //     f << c << std::endl;
            // }
            // f.close();

        
        if (col == 1) {
            std::vector<unsigned char> g1_col1(g1_col);
            std::vector<unsigned char> g2_col1(g2_col);
        }
        

        std::vector<bool> g1_new_sets;
        std::vector<bool> g2_new_sets;
        std::vector<bool> test_g2_new_sets;
        // https://github.com/facebook/folly/blob/master/folly/docs/FBVector.md
        g1_new_sets.reserve(g1_sets.size() * 1.5);
        g2_new_sets.reserve(g2_sets.size() * 1.5);
        test_g2_new_sets.reserve(g2_sets.size() * 1.5);
        
        std::cout << "going to refine sets in bool vectors of sizes " << g1_sets.size() << ", " << g2_sets.size() << ";   " << std::endl << "calling refine_sets() with column " << col << std::endl;
         int startgettime = getMilliCount();

        
        //std::vector<unsigned char> test_g2_col(kmers_reader.size());
        std::vector<unsigned char> test_g2_col();

    
        //for(int i=0; i < test_g2_col.size(); ++i){
            //test_g2_col.at(i) = kmer_to_string(*kmers_reader, g1.k, g1.k)[g1.k-col];
            //++kmers_reader;
        //}
        //kmers_reader.rewind();


        refine_sets(g1_col, g2_col, /*test_g2_col,*/ g1_sets, g2_sets, /*test_g2_sets,*/ col, g1_new_sets, g2_new_sets, /*test_g2_new_sets,*/ L);
        int delta = getMilliSpan(startgettime);
        std::cout << "refine_sets() completed in " << delta << " milliseconds." << std::endl << std::flush;    
        
        std::cout <<"    got back new vectors of bools of sizes " << g1_new_sets.size() << ", " << g2_new_sets.size() << ";   " << std::endl << std::endl;
        // assert(g1_col.size() == validate(g1_new_sets));
        // assert(g2_col.size() == validate(g2_new_sets));
        // std::cout << "validation passed g1: " << validate(g1_new_sets)
        //           << " g2: " << validate(g2_new_sets) << std::endl;
        //FIXME: avoid copying values in this next part
        g1_sets.assign(g1_new_sets.begin(), g1_new_sets.end()); 
        g2_sets.assign(g2_new_sets.begin(), g2_new_sets.end()); 

        if (col == g1.k - 2) {
            g1_flagsets.insert(g1_flagsets.end(), g1_sets.begin(), g1_sets.end());
            g2_flagsets.insert(g2_flagsets.end(), g2_sets.begin(), g2_sets.end());
        }

        float curplantime = getMilliSpan(planstarttime);
        colno++;
        std::cout << curplantime / 1000.0 / 60.0 << " minutes elapsed in planning phase." << std::endl;
        std::cout << ((curplantime / 1000.0 / 60.0) / (float)colno ) * (g1.k - colno) << " microsoft minutes remaining in planning phase." << std::endl;
    }
    //std::cout << "--------------------------" << std::endl;
    //
    //int count=0;
    //for(int i=0; i < g1_sets.size(); ++i) {
        //std::cout << g1_sets.at(i);
        //++count;
    //}
    //std::cout << ", " << count << std::endl;

    //for(int i=0; i< L.size(); ++i) {

        //std::cout << L.at(i);
    //}
    //std::cout << std::endl;

    //count=0;
    //for(int i=0; i < g2_sets.size(); ++i) {
        //std::cout << g2_sets.at(i);
        //++count;
    //}
    //std::cout << ", " << count << std::endl;

    //
    // Execution phase
    //
    
    //FIXME: find some other name than pointer, since really not
    uint64_t g1_ptr = 0; // tracks position in EBWT(g)_1
    uint64_t g2_ptr = 0; // tracks position in EBWT(g)_2

    uint64_t g1_set_ptr = 0; // tracks position in P_1
    uint64_t g2_set_ptr = 0; // tracks position in P_2
    
    uint64_t out_ptr = 0; // tracks position in EBWT(g)_M

    Flags flags(g1_flagsets, g2_flagsets);
    //std::cout << "flags: ";
    //dump_range(0, g1_flagsets.size(), g1_flagsets);
    std::cout << std::endl;
    //dump_range(0, g2_flagsets.size(), g2_flagsets);
    std::cout << std::endl;


    std::set<std::string> predecessor_dummy_nodes;
    std::set<std::string> successor_dummy_nodes;
        

    std::map<char, uint64_t> ntcounts;
    ntcounts['$'] = 0;
    ntcounts['A'] = 0;
    ntcounts['C'] = 0;
    ntcounts['G'] = 0;
    ntcounts['T'] = 0;




    // Lcol is just for debugging
    std::vector<bool> Lcol(g1.num_edges(), false); //DEBUG ONLY
    fill_Lcol(g1, Lcol);                           //DEBUG ONLY

    std::cout << "debug Lcol.size = " << Lcol.size() << " L.size = " << L.size() << std::endl;

    // BEGIN output stuff

    std::string outfile("deleted");
    string temp_edge_file = outfile + ".w.temp";
    std::cerr << "Creating file " << std::endl << std::flush;
    std::ofstream ef(temp_edge_file);
//    stxxl::syscall_file edge_file(temp_edge_file, stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT | stxxl::file::TRUNC);
//    edge_file.set_size(L.size()/*kmers_a.size() + outgoing_dummies_q.size() + outgoing_dummies_v.size() + incoming_dummies_v.size()*/);

    std::cerr << "Constructing vector bound to file" << std::endl << std::flush;
//    stxxl::vector<uint8_t> * output = new stxxl::vector<uint8_t>(&edge_file);
    std::cerr << "Creating writer bound to vector " << std::endl << std::flush;
    //  typename stxxl::vector<uint8_t>::bufwriter_type * edge_writer = new stxxl::vector<uint8_t>::bufwriter_type(*output);
    //char w_idx = 0;
    std::cerr << "writing to writer " << std::endl << std::flush;    

    // END output stuff
    std::ofstream pf("deleted.plan");

    //size_t g1_num_nodes = g1.num_nodes();
    size_t g1_num_edges = g1.num_edges();
    //std::vector<int> in_edges(g1_num_nodes);            //Tracks how many incoming edges remain on a node. Init to indegree
    //std::vector<int> out_edges(g1_num_nodes);           //Tracks how many outgoing edges remain on a node. Init to outdegree
    std::map<size_t, size_t> in_degrees;
    std::map<size_t, size_t> out_degrees;
    std::vector<ssize_t> edge_mapping(g1_num_edges, -1);    //Tracks where edge from g1 ends up in difference graph. Init to -1.
    std::vector<char> edges;                          //Representation of E-BWT that will be used to construct wavelet tree.

    //for(size_t i=0; i < g1.num_nodes(); ++i) {
        //in_edges.at(i) = g1.indegree(i);
        //out_edges.at(i) = g1.outdegree(i);
    //}

    std::cout << "g1 sets" << std::endl;
    for(size_t i=0; i < g1_sets.size(); ++i) {
        std::cout << g1_sets.at(i);
    }
    std::cout << std::endl;

    std::cout << "g2 sets" << std::endl;
    for(size_t i=0; i < g2_sets.size(); ++i) {
        std::cout << g2_sets.at(i);
    }
    std::cout << std::endl;


    do {
        // find the end of the next (equivalence class) set
        while (!g1_sets[g1_set_ptr]) ++g1_set_ptr;
        while (!g2_sets[g2_set_ptr]) ++g2_set_ptr;
        std::cout << g1_set_ptr;
        std::cout << g2_set_ptr;
        std::cout << "-------------" << std::endl;
        char symbol = 0; //This is what will be stored in E-BWT
        bool flag;
        // if the P_1 set is non-empty...
        if (g1_set_ptr > 0 && !g1_sets[g1_set_ptr - 1]) {
            
            symbol = g1._map_symbol(g1._strip_edge_flag(g1.m_edges[g1_ptr]));

            // ...and the P_2 set is empty, then "keep" the edge
            if (g2_set_ptr > 0 && g2_sets[g2_set_ptr - 1]) {
                ntcounts[g1._map_symbol(g1._symbol_access(g1_ptr))] += 1;
                flag = (flags.seen(symbol) > 0);
                //pf << (char)1; // include only g1 row
                //ef << combine(symbol, flag);
                edges.push_back(combine(symbol, flag));
                flags.add(symbol);
                std::cout << symbol << " added to edges with flag " << flag << std::endl;
                std::cout << "edge label: " << g1.edge_label((g1_ptr)) << std::endl;
                edge_mapping.at(g1_ptr) = out_ptr;
                ++out_ptr;

            }
            // ...and the P_2 set is non-empty, then we may want to "delete" this edge
            else /*(g2_set_ptr > 0 && !g2_sets[g2_set_ptr - 1])*/ {

                //TODO figure out how to update color matrix appropriately, but for now we can at least know whether or not to actually delete
                
                //First, advance iter to row of color matrix corresponding to this edge
                //size_t pos = color_iter.peek();
                //while(pos < (g1_ptr * num_colors) && pos != -1) {
                    //color_iter.advance();
                    //pos = color_iter.peek();
                //}
       
             
                //Now that we are in this edge's row, inspect the columns that belong to 
                //the minuend graph. If this edge is in some of those columns, do not delete 
                //because this edge is a component of other reads not being deleted!
                //bool should_del = true;
                //while(pos < ((g1_ptr+1) * num_colors) && pos != -1 && should_del) {
                    //if(pos < colors_start && pos > colors_end) {
                        //should_del = false;
                    //}
                    //color_iter.advance();
                    //pos = color_iter.peek();
                //}

                //if(should_del) {


                    std::cout << "edge label: " << g1.edge_label((g1_ptr)) << std::endl;

                    //std::cout << "succ label: " << g1.node_label(g1._edge_to_node(g1._forward(g1_ptr))) << std::endl;
                    //std::cout << in_edges.at(g1._edge_to_node(g1._forward(g1_ptr))) << std::endl;
                    //in_edges.at(g1._edge_to_node(g1._forward(g1_ptr))) -= 1;
                    //std::cout << in_edges.at(g1._edge_to_node(g1._forward(g1_ptr))) << std::endl;

                    //std::cout << "pred label: " << g1.node_label(g1._edge_to_node(g1_ptr)) << std::endl;
                    //std::cout << out_edges.at(g1._edge_to_node(g1_ptr)) << std::endl;
                    size_t g1_end_node = g1._edge_to_node(g1._forward(g1_ptr));
                    if(in_degrees.find(g1_end_node) == in_degrees.end()) {
                        in_degrees.insert( std::pair<size_t, size_t>(g1_end_node, g1.indegree(g1_end_node)-1) );
                        std::cout << "node label: " << g1.node_label(g1_end_node) << std::endl;
                        std::cout << "original indeg: " << in_degrees[g1_end_node] << std::endl;
                    }
                    else {
                        in_degrees[g1_end_node] -= 1;
                        std::cout << "node label: " << g1.node_label(g1_end_node) << std::endl;
                        std::cout << "new indeg: " << in_degrees[g1_end_node] << std::endl;
                    }

                    size_t g1_start_node = g1._edge_to_node(g1_ptr);
                    if(out_degrees.find(g1_start_node) == out_degrees.end()) {
                        out_degrees.insert( std::pair<size_t, size_t>(g1_start_node, g1.outdegree(g1_start_node)-1) );
                        std::cout << "node label: " << g1.node_label(g1_start_node) << std::endl;
                        std::cout << "original outdeg: " << out_degrees[g1_start_node] << std::endl;
                    }
                    else {
                        out_degrees[g1_start_node] -= 1;
                        std::cout << "node label: " << g1.node_label(g1_start_node) << std::endl;
                        std::cout << "new outdeg: " << out_degrees[g1_start_node] << std::endl;
                    }
                    //std::cout << out_edges.at(g1._edge_to_node(g1_ptr)) << std::endl;

                    g2_ptr += 1;
                    flags.adv_g2();
                    //pf << (char)3; // include g1 row and g2 row
                    //pf << (char)0;
                //}
                //else {
                    //std::cout << "Not deleting " << g1.edge_label(g1_ptr) << " because of color matrix considerations." << std::endl;
                    //ntcounts[g1._map_symbol(g1._symbol_access(g1_ptr))] += 1;
                    //flag = (flags.seen(symbol) > 0);
                    //edges.push_back(combine(symbol, flag));
                    //flags.add(symbol);
                    //edge_mapping.at(g1_ptr) = out_ptr;
                    //++out_ptr;
                //}
            } 
            g1_ptr += 1;
            flags.adv_g1();
        } else { // else, P_2 MUST be non-empty
            g2_ptr += 1;
            flags.adv_g2();
            //pf << (char)2; // include only g2 row
            //pf << (char)0;
        }
        ++g1_set_ptr;
        ++g2_set_ptr;
        //std::cout << "------------------------------------" << std::endl;
    } while (g1_set_ptr != g1_sets.size()  && g2_set_ptr != g2_sets.size() );
    pf.close();

    //std::cout << "edges before dummies:" << std::endl;
    //for(size_t i=0; i < edges.size(); ++i) {
        //std::cout << (edges.at(i) >> 1) << ((char)edges.at(i) & 1) << std::endl;
    //}
    //std::cout << ntcounts['$'] << " "
              //<< ntcounts['$'] + ntcounts['A'] << " "
              //<< ntcounts['$'] + ntcounts['A'] + ntcounts['C']   << " "
              //<< ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  << " "
              //<< ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  + ntcounts['T'] << " "
              //<< std::endl;


    ssize_t i = edge_mapping.size() - 1;
    while(edge_mapping.at(i) == -1) {
        --i;
    }

    std::cout << "I IS " << i << std::endl;
    size_t nonneg_edge_idx = edge_mapping.at(i);
    for(; i >= 0; --i) {
        std::cout << "I IS " << i << std::endl;
        if(edge_mapping.at(i) == -1) {
            edge_mapping.at(i) = nonneg_edge_idx;
        }
        else {
            nonneg_edge_idx = edge_mapping.at(i);
        }
    }
 
    std::cout << "L:" << std::endl;
    for(size_t i=0; i < L.size(); ++i) {
        std::cout << L.at(i) << std::endl;
    }

 
 
    //std::unordered_map<size_t, bool> edges_needing_incoming_dummies; 
    //std::unordered_map<size_t, bool> edges_needing_outgoing_dummies; 
    //size_t edge_idx;
    //for(size_t i=0; i < g1.num_nodes(); ++i) {
        //edge_idx = g1._node_to_edge(i);
        //std::cout << "node label: " << g1.node_label(i) << std::endl;
        //std::cout << in_degrees[i] << std::endl;
        //std::cout << out_degrees[i] << std::endl;
//
        //if(in_degrees[i] > 0 && out_degrees[i] == 0) {
            //std::cout << "node label: " << g1.node_label(i) << std::endl;
            //std::cout << "node needs outgoing dummy" << std::endl;

            //This edge has already been added
            //if( edges_needing_outgoing_dummies.find(edge_idx) != edges_needing_outgoing_dummies.end()) {
                //edges_needing_outgoing_dummies[edge_idx] = 1;
            //}
            //else {
                //edges_needing_outgoing_dummies[edge_idx] = 0;
            //}
        //}
        //else if(out_degrees[i] > 0 && in_degrees[i] == 0) {
            //std::cout << "node label: " << g1.node_label(i) << std::endl;
            //std::cout << "node needs incoming dummy chain" << std::endl;

            //This edge has already been added
            //if( edges_needing_incoming_dummies.find(edge_idx) != edges_needing_incoming_dummies.end()) {
                //edges_needing_outgoing_dummies[edge_idx] = 1;
            //}
            //else {
                //edges_needing_outgoing_dummies[edge_idx] = 0;
            //}
        //}
        //else {
            //std::cout << "node 'deleted'" << std::endl;
        //}

    //}

    for(size_t i=0; i < edges.size(); ++i) {
        std::cout << "EBWT at " << i << ": " << (edges.at(i) >> 1) << std::endl;
    }

    std::unordered_multimap<ssize_t, std::pair<std::string, char> > dummies_to_insert;
    for(size_t i=0; i < g1.num_nodes(); ++i) {
        if(in_degrees[i] > 0 && out_degrees[i] == 0) {
            std::cout << "node label: " << g1.node_label(i) << std::endl;
            std::cout << "node needs outgoing dummy" << std::endl;
            insert_outgoing_dummy(edges, L, ntcounts, edge_mapping, g1._node_to_edge(i), g1, dummies_to_insert);

            //std::cout << "edges after outgoing dummy:" << std::endl;
            //for(size_t i=0; i < edges.size(); ++i) {
                //std::cout << (edges.at(i) >> 1) << std::endl;//((char)edges.at(i) & 1) << std::endl;
            //}
        }
        else if(out_degrees[i] > 0 && in_degrees[i] == 0) {
            std::cout << "node label: " << g1.node_label(i) << std::endl;
            std::cout << "node needs incoming dummy chain" << std::endl;
            std::cout <<  "edge label: " << g1.edge_label(g1._node_to_edge(i)) << std::endl;
            insert_incoming_dummy_chain(edges, L, ntcounts, edge_mapping, g1._node_to_edge(i), g1, dummies_to_insert);

            //std::cout << "edges after incoming dummy:" << std::endl;
            //for(size_t i=0; i < edges.size(); ++i) {
                //std::cout << (edges.at(i) >> 1) << std::endl;//((char)edges.at(i) & 1) << std::endl;
            //}
        }
        else {
            //std::cout << "node 'deleted'" << std::endl;
        }
        std::cout << "==========================================================" << std::endl;
    }

    std::cout << "L:" << std::endl;
    for(size_t i=0; i < L.size(); ++i) {
        std::cout << L.at(i) << std::endl;
    }

    std::vector<char> edges_with_dummies;
    std::vector<bool> L_with_dummies; // node flags
    for(size_t i=0; i<edges.size(); ++i) {
        auto dummies_range = dummies_to_insert.equal_range(i);
    
        if(std::distance(dummies_range.first, dummies_range.second) == 1) {
            edges_with_dummies.push_back( ((*(dummies_range.first)).second).second << 1 );
            L_with_dummies.push_back(false);
        }
        else if(std::distance(dummies_range.first, dummies_range.second) > 1) {
            std::vector<std::pair<std::string, char> > sorted_dummies;
            for(auto itr = dummies_range.first; itr != dummies_range.second; ++itr) {
                std::cout << (*itr).second.first << std::endl;
                sorted_dummies.push_back((*itr).second);
            }
            std::sort(sorted_dummies.begin(), sorted_dummies.end(),
                [](const std::pair<std::string, char> a, const std::pair<std::string, char> b) -> bool {
                    return a.first < b.first;
                });

            for(auto dummy_info : sorted_dummies) {
                edges_with_dummies.push_back(dummy_info.second << 1);
                L_with_dummies.push_back(false);
            }
        }

        edges_with_dummies.push_back(edges.at(i));
        L_with_dummies.push_back(L.at(i));
    }

    std::cout << "L with dummies:" << std::endl;
    for(size_t i=0; i < L_with_dummies.size(); ++i) {
        std::cout << L_with_dummies.at(i) << std::endl;
    }

    std::cout << "edges with dummies:" << std::endl;
    for(size_t i=0; i < edges_with_dummies.size(); ++i) {
        std::cout << (edges_with_dummies.at(i) >> 1) << std::endl;//((char)edges.at(i) & 1) << std::endl;
    }

    std::cout << "set pointers " << g1_set_ptr << "/" <<g1_sets.size() << " " << g2_set_ptr << "/" << g2_sets.size() << std::endl;
    std::cout << "EBWT(G) ptrs  _1: " << g1_ptr << " _2: " << g2_ptr << " _M: " << out_ptr << std::endl;
    std::cout << "elements in family: " << validate(g1_sets) << " " << validate(g2_sets) << std::endl;
    std::cout << "max size in family g1_sets: " << maxsize(g1_sets) << std::endl;
    std::cout << "max size in family g2_sets: " << maxsize(g2_sets) << std::endl;
    assert (g1_set_ptr == g1_sets.size() );
    assert ( g2_set_ptr == g2_sets.size() );

    std::cout << "m_symbol_ends for merged graph: ";
    std::cout << ntcounts['$'] << " "
              << ntcounts['$'] + ntcounts['A'] << " "
              << ntcounts['$'] + ntcounts['A'] + ntcounts['C']   << " "
              << ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  << " "
              << ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  + ntcounts['T'] << " "
              << std::endl;


    for(size_t i=0; i < edges.size(); ++i) {
        ef << edges.at(i);
    }

    ef.close();

    std::cerr << "constructing wavelet tree " << std::endl << std::flush;
    typedef sdsl::wt_huff<sdsl::rrr_vector<63>> wt_t;
    wt_t e_bwt;
    construct(e_bwt, temp_edge_file, 1);
    std::cerr << "closing file " << std::endl << std::flush;    

    std::cerr << "creating succinct de bruijn graph " << std::endl << std::flush;
    size_t n = L.size(); // number of elements

    
    size_t m = 0; // number of 1's
    for (size_t i = 0; i < L.size(); i++)
        m += L[i];
    
    std::cout << "n: " << n << " m: " << m << std::endl;
    sdsl::sd_vector_builder *b_builder = new sdsl::sd_vector_builder(n, m);
    for (size_t i = 0; i < L.size(); i++)
        if (L[i]) b_builder->set(i);
        
    sd_vector<> node_bv(*b_builder);
    delete b_builder;
    array<size_t, 5> counts{
        ntcounts['$'] ,
        ntcounts['$'] + ntcounts['A'] ,
        ntcounts['$'] + ntcounts['A'] + ntcounts['C'] ,
        ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  ,
        ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  + ntcounts['T'] };
    sd_vector<> dum_pos_bv;
    vector<kmer_t> dummies;

    debruijn_graph_shifted<> dbgm = cosmo::make_dbg<debruijn_graph_shifted<>>()(g1.k, node_bv, e_bwt, counts, "$ACGT", dum_pos_bv, dummies);
    cerr << "k             : " << dbgm.k << endl;
    cerr << "num_nodes()   : " << dbgm.num_nodes() << endl;
    cerr << "num_edges()   : " << dbgm.num_edges() << endl;
    cerr << "Total size    : " << size_in_mega_bytes(dbgm) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(dbgm) << " Bits" << endl;

    std::cerr << "storing succinct de bruijn graph to file" << std::endl << std::flush;    
    sdsl::store_to_file(dbgm, "deleted.dbg");

    return 0;
}


void parse_arguments(int argc, char **argv, test_vd_parameters_t & params){
    TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
    TCLAP::UnlabeledValueArg<std::string> target_filename_arg("target_graph", ".dbg file.", true, "", "target_graph_file", cmd);
    TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input_graph", ".dbg file.", true, "", "input_graph_file", cmd);
    cmd.parse( argc, argv );

    params.target_filename  = target_filename_arg.getValue();
    params.input_filename  = input_filename_arg.getValue();
}
void find_kmers_to_delete(const debruijn_graph_shifted<> &target_dbg, const debruijn_graph_shifted<>& input_dbg/*, std::vector<std::string>& resultant_kmers*/)
{

    std::string target_kmers_file("target_kmers");
    std::cerr << "Creating target kmers file " << std::endl << std::flush;
    std::ofstream tk_file(target_kmers_file);

    for(size_t edge_idx=0; edge_idx < target_dbg.num_edges(); ++edge_idx) {
        tk_file << target_dbg.edge_label(edge_idx) << std::endl;
    }
    tk_file.close();


    std::string input_kmers_file("input_kmers");
    std::cerr << "Creating input kmers file " << std::endl << std::flush;
    std::ofstream ik_file(input_kmers_file);


    for(size_t edge_idx=0; edge_idx < input_dbg.num_edges(); ++edge_idx) {
        ik_file << input_dbg.edge_label(edge_idx) << std::endl;
    }
    ik_file.close();
}
void output_final_kmers(const debruijn_graph_shifted<>& resultant_dbg/*, const std::vector<string>& expected_kmers*/)
{
    std::vector<std::string> deleted_kmers;
    std::string deleted_kmers_file("deleted_kmers");
    std::cerr << "Creating deleted kmers file " << std::endl << std::flush;
    std::ofstream dk_file(deleted_kmers_file);
    for(size_t edge_idx=0; edge_idx < resultant_dbg.num_edges(); ++edge_idx) {
        dk_file << resultant_dbg.edge_label(edge_idx) << std::endl;
    }
    dk_file.close();
}
int main(int argc, char** argv) {
    //Load in two graphs

    test_vd_parameters_t p;
    parse_arguments(argc, argv, p);

    std::cout << "Loading target debruijn graph" << p.target_filename << std::endl;
    debruijn_graph_shifted<> target_dbg;
    load_from_file(target_dbg, p.target_filename);
    // dumpcolumns(target_dbg, p);
    // dump_edges(target_dbg);
    cerr << "k             : " << target_dbg.k << endl;
    cerr << "num_nodes()   : " << target_dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << target_dbg.num_edges() << endl;
    cerr << "Total size    : " << size_in_mega_bytes(target_dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(target_dbg) << " Bits" << endl <<endl;
 
    std::cout << "Loading input debruijn graph" << p.input_filename << std::endl;
    debruijn_graph_shifted<> input_dbg;
    load_from_file(input_dbg, p.input_filename);
    // dumpcolumns(input_dbg, p);
    // dump_edges(input_dbg);
    cerr << "k             : " << input_dbg.k << endl;
    cerr << "num_nodes()   : " << input_dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << input_dbg.num_edges() << endl;
    cerr << "Total size    : " << size_in_mega_bytes(input_dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(input_dbg) << " Bits" << endl;

    assert(target_dbg.k == input_dbg.k);

    //Enumerate all kmers, and figure out which ones will be in resultant.
    //std::vector<std::string> resultant_kmers;

    find_kmers_to_delete(target_dbg, input_dbg/*, resultant_kmers*/);

    //Delete input graph from target graph to get resultant graph
    dbg_delete(target_dbg, input_dbg);

    //Check what kmers are present in the resultant graph
    std::string resultant_filename("deleted.dbg");
    std::cout << "Loading resultant debruijn graph" << resultant_filename << std::endl;
    debruijn_graph_shifted<> resultant_dbg;
    load_from_file(resultant_dbg, resultant_filename);
    // dumpcolumns(resultant_dbg, p);
    // dump_edges(resultant_dbg);
    cerr << "k             : " << resultant_dbg.k << endl;
    cerr << "num_nodes()   : " << resultant_dbg.num_nodes() << endl;
    cerr << "num_edges()   : " << resultant_dbg.num_edges() << endl;
    cerr << "Total size    : " << size_in_mega_bytes(resultant_dbg) << " MB" << endl;
    cerr << "Bits per edge : " << bits_per_element(resultant_dbg) << " Bits" << endl;

    output_final_kmers(resultant_dbg/*, resultant_kmers*/);
}

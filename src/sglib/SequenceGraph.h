//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#ifndef SG_SEQUENCEGRAPH_HPP
#define SG_SEQUENCEGRAPH_HPP

#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <iostream>
#include <array>
#include <sglib/factories/KMerNeighbourFactory.h>
#include <sglib/readers/StringReader.h>
#include <sglib/readers/FileReader.h>
#include <sglib/SMR.h>
#include <deque>

typedef int64_t sgNodeID_t; //first node is 1; negatives are RC

enum sgNodeStatus_t {sgNodeActive,sgNodeDeleted};

class Node{
public:
    Node(std::string _seq) : sequence(_seq),status(sgNodeActive){};
    std::string sequence;
    sgNodeStatus_t status;
    bool is_canonical();
    void make_rc();
};

class Link{
public:
    Link( sgNodeID_t _src, sgNodeID_t _dst, int32_t _dist) : source(_src), dest(_dst), dist(_dist) {};
    sgNodeID_t source,dest;
    int32_t dist;

    bool operator==( const  Link);
    bool operator<(const Link)const;

    friend std::ostream &operator<<(std::ostream &os, const Link &link) {
        os << link.source << " -> " << link.dest;
        return os;
    }

};

class SequenceGraphPath;
class SequenceSubGraph;
class SequenceGraph {
public:
    SequenceGraph(){};

    // Builds a DBG from a vector of seqs
    SequenceGraph(const std::vector<std::string> &seqs, unsigned int k, unsigned int min_cov = 1){
        /*
         * Get context augmented k-mers from the sequences filtered by min_cov
         */
        SMR<KmerNeighbour,
        kmerNeighbourFactory<FastaRecord>,
        StringReader<FastaRecord>,
        FastaRecord, StringReaderParams, KMerNeighbourFactoryParams> kmerNeigh_SMR(StringReaderParams{seqs,1}, {k}, {4*GB, min_cov, 10000, "default"});

        auto kmers = kmerNeigh_SMR.process_from_memory();
        std::unordered_map<std::string, KmerNeighbour> kmerDict;
        std::transform(kmers.begin(),kmers.end(), std::inserter(kmerDict, kmerDict.begin()), [](KmerNeighbour &kmn){
            return std::pair<std::string, KmerNeighbour>(kmn.kmer, kmn);
        });
        /*
         * Cleanup the kmer neighbours if they don't exist
         */
        for (auto &km : kmers) {
            int context(km.pre_post_context);
            // Have fwd neighs
            if (km.canExtendFwd()) {
                std::string kmer(km.kmer.substr(1, km.kmer.size())+'N');
                for (unsigned succCode = 0; succCode < 4u; ++succCode) {
                    if (km.hasFwdNeighbour(succCode)) {
                        kmer.back() = "ACTG"[succCode];
                        if (kmerDict.find(kmer) == kmerDict.end())
                            context &= ~(1 << succCode);
                    }
                }
            }
            // Have bwd neighs
            if (km.canExtendBwd()) {
                std::string kmer('N'+km.kmer.substr(0,km.kmer.size()-1));
                for (unsigned succCode = 0; succCode < 4u; ++succCode) {
                    if (km.hasBwdNeighbour(succCode)) {
                        kmer.front() = "ACTG"[succCode];
                        if (kmerDict.find(kmer) == kmerDict.end())
                            context &= ~(1 << (succCode + 4u));
                    }
                }
            }
            km.pre_post_context = context;
        }

        /*
         * For each node, if node has just 1 neighbour extend until more than 1.
         * Remove nodes touched, and add a new node.
         * Restart from first non-touched.
         */

        // TODO: Finish the extend fwd, bwd functions for the nodes
        for (auto &km: kmers) {
            std::deque<char> node(km.kmer.begin(), km.kmer.end());
            // If can go fwd
            bool canContinue(false);
            if (km.canExtendFwd()) {
                // Can go bwd? leave
                if (km.canExtendBwd()) break;
                do {
                    canContinue = km.extendFwd(node, kmerDict);
                } while (canContinue);
            }
            else if (km.canExtendBwd()) {
                do {
                    canContinue = km.extendBwd(node, kmerDict);
                } while (canContinue);
            }
            // make singleton?
            add_node(Node(std::string(node.begin(), node.end())));
        }
    }
    //=== I/O functions ===
    void load_from_gfa(std::string filename);
    void write_to_gfa(std::string filename);


    //=== graph operations ===
    sgNodeID_t add_node(Node n);
    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);

    std::vector<Link> get_fw_links( sgNodeID_t n);
    std::vector<Link> get_bw_links( sgNodeID_t n);

    /*
     * Connected components, (TODO) optionally breaking up in repeats, nodes that class as repeats will be returned on their own
     */
    std::vector<std::vector<sgNodeID_t>> connected_components (int max_nr_totalinks=0, int max_nr_dirlinks=0, int min_rsize=0); //TODO: --> enable extra breaks in repeats

    // find bubbles in component of graph
    std::vector<std::vector<sgNodeID_t >> find_bubbles(std::vector<sgNodeID_t>);

    std::vector<sgNodeID_t> depth_first_search(sgNodeID_t node, unsigned int size_limit, unsigned int edge_limit, std::set<sgNodeID_t> tabu={});

    std::vector<sgNodeID_t> breath_first_search(std::vector<sgNodeID_t> &nodes, unsigned int size_limit);

    // remove_node
    void remove_node(sgNodeID_t);
    // remove_link
    void remove_link(sgNodeID_t source, sgNodeID_t dest);
    //These two need to mark expanded edges, and transfer read maps and unique kmers for non-expanded, but just read map for expanded.

    void join_path(const SequenceGraphPath p, bool consume_nodes=true);
    // expand_path --> creates an edge with the consensus of a path, eliminates old nodes if only in path and unused edges
    void join_all_unitigs();
    std::vector<SequenceGraphPath> get_all_unitigs(uint16_t min_nodes);
    // simplify --> executes expand_path on every multi-sequence unitig


    // tip_clip -> eliminates tips.


    //void explode_node( sgNodeID_t node, uint16_t k);
    //void explode_all_nodes ();
    //void collapse_identical_nodes ();

    //later
    //project spectra and use for flow


    //=== internal variables ===
    std::vector<sgNodeID_t> oldnames_to_nodes(std::string _oldnames);
    std::vector<Node> nodes;
    std::vector<std::vector<Link>> links;
    std::unordered_map<std::string,sgNodeID_t> oldnames_to_ids;
    std::vector<std::string> oldnames;
    std::string filename,fasta_filename;

    void consume_nodes(const SequenceGraphPath &p, const std::set<sgNodeID_t> &pnodes);

    std::vector<sgNodeID_t > find_canonical_repeats();
    bool is_loop(std::array<sgNodeID_t, 4> nodes) {
        for (auto j = 0; j < 3; ++j)
            for (auto i = j + 1; i < 4; ++i)
                if (nodes[i] == nodes[j] or nodes[i] == -nodes[j]) return true; //looping node
        return false;
    }
};


class SequenceGraphPath {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceGraphPath(SequenceGraph & _sg, const std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    std::string get_fasta_header();
    std::string get_sequence();
    void reverse();
    bool is_canonical();

private:
    SequenceGraph& sg;
};

class SequenceSubGraph {
public:
    std::vector<sgNodeID_t> nodes;
    explicit SequenceSubGraph(SequenceGraph & _sg, std::vector<sgNodeID_t> _nodes={})  : sg(_sg) ,nodes(_nodes) {};
    SequenceGraphPath make_path(); //returns empty path if not linear

private:
    SequenceGraph& sg;
};

#endif //SG_SEQUENCEGRAPH_HPP

//
// Created by Bernardo Clavijo (EI) on 25/05/2018.
//

#ifndef BSG_LINKAGEDIGRAPH_HPP
#define BSG_LINKAGEDIGRAPH_HPP


#include <sglib/graph/SequenceGraph.hpp>

/**
 * Contains a graph and links generated between nodes on that graph
 * Each link is represented on each node it appears as -A,B on A and -B,A on B
 */
class LinkageDiGraph {
public:
    LinkageDiGraph(SequenceGraph & _sg): sg(_sg){};

    void add_link( sgNodeID_t source, sgNodeID_t dest, int32_t d);
    void add_links(const LinkageDiGraph &other);

    void remove_link(sgNodeID_t source, sgNodeID_t dest);

    std::vector<Link> get_fw_links( sgNodeID_t n) const;
    std::vector<Link> get_bw_links( sgNodeID_t n) const;
    std::set<sgNodeID_t> fw_reached_nodes(sgNodeID_t n, int radius) const;


    void remove_transitive_links(int radius);
    void report_connectivity();
    //void solve();

    std::vector<std::vector<sgNodeID_t>> get_all_lines(uint16_t min_nodes) const;

    void dump_to_text(std::string filename);
    void load_from_text(std::string filename);

    SequenceGraph & sg;
    std::vector<std::vector<Link>> links;

};
#endif //BSG_LINKAGEDIGRAPH_HPP

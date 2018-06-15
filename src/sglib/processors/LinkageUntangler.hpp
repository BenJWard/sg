//
// Created by Bernardo Clavijo (EI) on 28/05/2018.
//

#ifndef BSG_LINKAGEUNTANGLER_HPP
#define BSG_LINKAGEUNTANGLER_HPP


#include <sglib/WorkSpace.hpp>
#include <sglib/LinkageDiGraph.hpp>

class LinkageUntangler {
    std::unordered_map<unsigned int, unsigned int> getValidTags();
    bool share_tags(const std::vector<uint64_t> &node_hash, const sgNodeID_t &n, const sgNodeID_t &m);
    void create_tag_buckets(std::vector<std::vector<bsg10xTag>> &node_tags, std::vector<uint64_t> &node_hash,
                                std::unordered_map<unsigned int, unsigned int> &valid_tags);
    void assign_tags_to_nodes(std::vector<std::vector<bsg10xTag>> &node_tags, std::unordered_map<unsigned int, unsigned int> &valid_tags);
public:

    explicit LinkageUntangler(WorkSpace & _ws): ws(_ws) { clear_node_selection();};

    //Node selection methods
    void clear_node_selection();
    void report_node_selection();
    void select_nodes_by_size_and_ci( uint64_t min_size, float min_ci, float max_ci);
    std::set<std::pair<sgNodeID_t, sgNodeID_t >> get_HSPNPs(uint64_t min_size, float min_ci, float max_ci);
    void select_nodes_by_HSPNPs(uint64_t min_size, float min_ci, float max_ci);
    void select_frontiers_by_size_and_ci();

    //Linkage creation methods (work on selected nodes)
    LinkageDiGraph make_topology_linkage(int radius);
    LinkageDiGraph make_paired_linkage(int min_reads);
    LinkageDiGraph make_tag_linkage(int min_reads,float end_perc=.3);


    //Problem localisation methods


    WorkSpace &ws;
    std::vector<bool> selected_nodes;
    std::vector<bool> frontier_nodes;
};


#endif //BSG_LINKAGEUNTANGLER_HPP

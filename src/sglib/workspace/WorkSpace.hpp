//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#ifndef BSG_WORKSPACE_HPP
#define BSG_WORKSPACE_HPP


#include <sglib/datastores/LinkedReadsDatastore.hpp>
#include <sglib/mappers/LinkedReadMapper.hpp>
#include <sglib/datastores/PairedReadsDatastore.hpp>
#include <sglib/mappers/PairedReadMapper.hpp>
#include <sglib/datastores/PathsDatastore.hpp>
#include <sglib/mappers/LongReadMapper.hpp>
#include "sglib/graph/SequenceGraph.hpp"
#include "sglib/KmerCompressionIndex.hpp"

class LogEntry{
public:
    LogEntry(std::time_t t, std::string v, std::string tx):timestamp(t),bsg_version(std::move(v)),log_text(std::move(tx)){};
    std::time_t timestamp;
    std::string bsg_version;
    std::string log_text;
};
class WorkSpace {
    //All status classes are public, treat them with care anyway ;)
    SequenceGraph sg;
    std::vector<LinkedReadsDatastore> linked_read_datastores;
    std::vector<LinkedReadMapper> linked_read_mappers;
    std::vector<PairedReadsDatastore> paired_read_datastores;
    std::vector<PairedReadMapper> paired_read_mappers;
    std::vector<LongReadMapper> long_read_mappers;
    std::vector<LongReadsDatastore> long_read_datastores;
    std::vector<PathsDatastore> path_datastores;
    KmerCompressionIndex kci;
    std::string verbose_log="";
    std::vector<LogEntry> log;

public:
    WorkSpace():kci(sg){};
    WorkSpace(const WorkSpace& that) = delete; //we definitely do not want copy constructors here, thank you
    void print_log();

    void add_log_entry(std::string text);

    void dump_to_disk(std::string filename);

    void load_from_disk(std::string filename,bool log_only=false);

    //general operations

    void remap_all();

    //Projected operations with info from the graph

    std::vector<sgNodeID_t>
    select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint32_t min_tags, uint32_t max_tags, float min_ci, float max_ci);

    KmerCompressionIndex& getKCI() {return kci;}
    SequenceGraph& getGraph() {return sg;}

    std::vector<LinkedReadMapper>& getLinkedReadMappers() {return linked_read_mappers;}
    std::vector<LinkedReadsDatastore>& getLinkedReadDatastores() {return linked_read_datastores;}
    std::vector<PairedReadMapper>& getPairedReadMappers() {return paired_read_mappers;}
    std::vector<PairedReadsDatastore>& getPairedReadDatastores() {return paired_read_datastores;}
    std::vector<LongReadMapper>& getLongReadMappers() {return long_read_mappers;}
    std::vector<LongReadsDatastore>& getLongReadDatastores() {return long_read_datastores;}
    std::vector<PathsDatastore>& getPathsDatastore() {return path_datastores;}
    std::vector<std::string> read_counts_header;
};


#endif //BSG_WORKSPACE_HPP

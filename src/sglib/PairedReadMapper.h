//
// Created by Bernardo Clavijo (EI) on 03/11/2017.
//

#ifndef SG_PAIREDREADMAPPER_HPP
#define SG_PAIREDREADMAPPER_HPP

#include <map>
#include <unordered_map>
#include <strings.h>
#include <vector>
#include <unordered_set>
#include <ostream>
#include <sglib/Types.h>
class SequenceGraph;

struct graphPosition{
    sgNodeID_t node;
    uint32_t pos;
};
class ReadMapping {
public:
    ReadMapping(){
        //just clean the structure, OSX doesn't give you clean memory
        bzero(this, sizeof(ReadMapping));
    }
    bool operator==(const ReadMapping &other){
        return this==&other;
    };
    bool operator<(const ReadMapping &other) const {
        return std::tie(node, read_id) < std::tie(other.node, other.read_id);
    };
    void merge(const ReadMapping &other){};
    friend std::ostream& operator<<(std::ostream& os, const ReadMapping& rm) {
        os << rm.node << "\t" << rm.unique_matches;
        return os;
    }

    sgNodeID_t node = 0;
    uint64_t read_id = 0;
    int32_t first_pos = 0;
    int32_t last_pos = 0;
    int32_t unique_matches = 0;
    bool rev=false;

};

/**
 * @brief A mapper for different kinds of paired short reads (Illumina).
 *
 * The mapper supports so far PE and 10x reads.
 * 10x reads must be trimmed and have the tags at the end of the read name.
 */
class PairedReadMapper {
public:

    enum prmReadType {prmPE, prmLMP, prm10x, prmLR};
    const std::vector<std::string> prmReadTypeDesc = {"Paired End", "Long Mate Pair", "10x Linked Reads", "Long Reads"};

    PairedReadMapper(SequenceGraph &_sg);;
    void map_reads(std::string , std::string , PairedReadMapper::prmReadType , uint64_t );
    void map_reads(std::string, uint64_t);
    void remove_obsolete_mappings();
    void remap_reads(std::unordered_set<uint64_t> const &  reads_to_remap={});
    uint64_t process_reads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t, bool tags=false, std::unordered_set<uint64_t> const & reads_to_remap={});
    uint64_t process_longreads_from_file(uint8_t, uint16_t, std::unordered_map<uint64_t , graphPosition> &, std::string , uint64_t );
    void save_to_disk(std::string filename);
    void load_from_disk(std::string filename);
    void print_stats();

    SequenceGraph & sg;
    std::string read1filename,read2filename;
    prmReadType readType;
    uint64_t memlimit;
    std::vector<std::vector<ReadMapping>> reads_in_node;
    std::vector<sgNodeID_t> read_to_node;//id of the main node if mapped, set to 0 to remap on next process
    std::vector<prm10xTag_t> read_to_tag;
};

#endif //SG_PAIREDREADMAPPER_HPP

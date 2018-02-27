//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/mappers/LongReadMapper.h>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::string gfa_filename, bubble_contigs_filename, output_prefix, long_reads;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(4);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
//@formatter:off
    cxxopts::Options options("map-lr", "LongRead Mapper");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>(),"")
            ("g,gfa", "input gfa file", cxxopts::value<std::string>(gfa_filename), "filepath")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("log_level", "output log level", cxxopts::value<unsigned int>(log_level), "uint");
    options.add_options("Long Read Options")
            ("r,long_reads", "input long reads", cxxopts::value<std::string>(long_reads), "filepath")
            ("d,dump_to","dump mapped reads to file",cxxopts::value<std::string>(dump_mapped), "filepath")
            ("l,load_from", "load mapped reads from file", cxxopts::value<std::string>(load_mapped), "filepath")
            ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb)->default_value("4"), "GB");
//@formatter:on
    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("g") != 1 or result.count("o") != 1 ) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
        }
        if (long_reads.empty()) {
            throw cxxopts::OptionException(" please specify a long reads file");

        }
    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        std::cout << options.help({""}) << std::endl;
        exit(1);
    }
    if (!sglib::check_or_create_directory(output_prefix)) {
        exit(1);
    }

    std::vector<std::string> seqs;
    SequenceGraph dbg(seqs, 100, 5);
}

//
// Created by Luis Yanes (EI) on 26/01/2018.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <sglib/factories/ContigBlockFactory.h>
#include <sglib/SequenceGraph.h>
#include "cxxopts.hpp"


int main(int argc, char * argv[]) {
    std::string in_file, output_prefix;
    std::string dump_mapped, load_mapped;
    unsigned int log_level(4);
    uint64_t max_mem_gb(4);
    bool stats_only(false);
//@formatter:off
    cxxopts::Options options("map-lr", "LongRead Mapper");
    options.add_options()
            ("help", "Print help", cxxopts::value<std::string>(),"")
            ("i,input", "input sequences file", cxxopts::value<std::string>(in_file), "filepath")
            ("o,output", "output file prefix", cxxopts::value<std::string>(output_prefix), "path")
            ("log_level", "output log level", cxxopts::value<unsigned int>(log_level), "uint");
            ("max_mem", "maximum_memory when mapping (GB, default: 4)", cxxopts::value<uint64_t>(max_mem_gb)->default_value("4"), "GB");
//@formatter:on
    try {
        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help({""}) << std::endl;
            exit(0);
        }

        if (result.count("i") != 1 or result.count("o") != 1 ) {
            throw cxxopts::OptionException(" please specify input files and output prefix");
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
    std::ifstream inputSeqs(in_file);
    std::string line;
    while (std::getline(inputSeqs, line)) {
        seqs.push_back(line);
    }
    SequenceGraph dbg(seqs, 8, 1);

    dbg.write_to_gfa(output_prefix+"_dbg.gfa");
}

//
// Created by Luis Yanes (EI) on 26/02/2018.
//

#ifndef BSG_STRINGREADER_H
#define BSG_STRINGREADER_H

#include <string>
#include <vector>
#include "Common.h"

struct StringReaderParams{
    explicit StringReaderParams(const std::vector<std::string> &ings, unsigned int minl=0) : strings(ings), min_length(minl) {}
    const std::vector<std::string> &strings;
    unsigned int min_length = 0;

};

template<typename FileRecord>
class StringReader {
public:
    explicit StringReader(StringReaderParams params) : params(params), numRecords(1), strings(params.strings) {
        min_length=params.min_length;//TODO: use this
    }

    bool next_record(FileRecord& rec) {
        if (numRecords<strings.size()) {
            rec.id = numRecords;
            rec.seq = strings[numRecords];
            rec.name = std::to_string(numRecords);
            ++numRecords;
            stats.totalLength += rec.seq.size();
            return true;
        } else return false;
    }
    ReaderStats getSummaryStatistics() {
        stats.totalRecords = numRecords-1;
        return stats;
    }
private:
    const std::vector<std::string> &strings;
    StringReaderParams params;
    uint32_t numRecords;
    ReaderStats stats;
    uint32_t min_length;

};


#endif //BSG_STRINGREADER_H

//
// Created by Luis Yanes (EI) on 15/09/2017.
//

#ifndef BSG_KMERNEIGHBOUR_H
#define BSG_KMERNEIGHBOUR_H

#include <vector>
#include <limits>
#include <tuple>
#include <deque>
#include <unordered_map>
#include "sglib/factories/KMerFactory.h"
struct KMerNeighbourFactoryParams {
    unsigned int k;
};

struct KmerNeighbour {
    KmerNeighbour() = default;
    KmerNeighbour(std::string kmer, int ctx) : kmer(kmer), pre_post_context(ctx), count(1) {}

    const bool operator<(const KmerNeighbour& other) const {
        return kmer<other.kmer;
    }

    const bool operator>(const KmerNeighbour &other) const {
        return kmer>other.kmer;
    }

    const bool operator==(const KmerNeighbour &other) const {
        return kmer==other.kmer;
    }
    void merge(const KmerNeighbour &other) {
        count += other.count;
        pre_post_context |= other.pre_post_context;
    }

    KmerNeighbour rc(){

    }

    friend std::ostream& operator<<(std::ostream& os, const KmerNeighbour& kmer) {
        os << kmer.kmer << "\t" << kmer.count << "\t" << kmer.pre_post_context;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, const KmerNeighbour& kmer) {
        is.read((char*)&kmer, sizeof(kmer));
        return is;
    }
    unsigned int bitCount(unsigned int x)
    {
        x = (((x >> 1) & 0b01010101010101010101010101010101)
             + x       & 0b01010101010101010101010101010101);
        x = (((x >> 2) & 0b00110011001100110011001100110011)
             + x       & 0b00110011001100110011001100110011);
        x = (((x >> 4) & 0b00001111000011110000111100001111)
             + x       & 0b00001111000011110000111100001111);
        x = (((x >> 8) & 0b00000000111111110000000011111111)
             + x       & 0b00000000111111110000000011111111);
        x = (((x >> 16)& 0b00000000000000001111111111111111)
             + x       & 0b00000000000000001111111111111111);
        return x;
    }


    bool canExtendFwd() const {
        auto bitcount(pre_post_context >> 4);
        return bitcount == 1 and active;
    }

    bool canExtendBwd() const {
        auto bitcount(pre_post_context & 0xf);
        return bitcount == 1 and active;
    }

    int getFwd() const {
        int bits(pre_post_context >> 4);
        int succCode=0;
        while (succCode<=4) {
            if (bits & 1 << succCode) return succCode;
            succCode++;
        }
    }

    int getBwd() const {
        int bits(pre_post_context & 0xf);
        int succCode=0;
        while (succCode<=4) {
            if (bits & 1 << succCode) return succCode;
            succCode++;
        }
    }
    bool hasFwdNeighbour(int neighbour) const {
        return ( (pre_post_context & 1 << neighbour) != 0);
    }
    bool hasBwdNeighbour(int neighbour) const {
        return ( (pre_post_context >> 4 & 1 << neighbour) != 0);
    }

    KmerNeighbour extendFwd(unsigned int k, std::deque<char> &node, const std::unordered_map<std::string, KmerNeighbour> &dict) {
        node.push_back("ACTG"[getFwd()]);
        this->active=false;
        std::string m(node.end()-k, node.end());
        auto km(dict.find(m));
        return km->second;
    }

    KmerNeighbour extendBwd(unsigned int k, std::deque<char> &node, const std::unordered_map<std::string, KmerNeighbour> &dict) {
        node.push_front("ACTG"[getBwd()]);
        this->active=false;
        std::string m(node.begin(), node.begin() + k);
        auto km(dict.find(m));
        return km->second;
    }

    std::string kmer = "";
    int pre_post_context = 0;
    uint8_t count = 1;
    bool active = true;
};

namespace std {
    template <>
    struct hash<KmerNeighbour> {
        size_t operator()(const KmerNeighbour& k) const {
            return std::hash<std::string>()(k.kmer);
        }
    };
}

template<typename FileRecord>
class kmerNeighbourFactory {

public:
    explicit kmerNeighbourFactory(KMerNeighbourFactoryParams params) : k(params.k) {
        fkmer.resize(k);
        rkmer.resize(k);
    }

    ~kmerNeighbourFactory() = default;

    void setFileRecord(FileRecord &rec) {
        currentRecord = rec;
        fkmer="";
        rkmer="";
        last_unknown=0;
    }

    // TODO: Adjust for when K is larger than what fits in uint64_t!
    const bool next_element(std::vector<KmerNeighbour> &mers) {
        uint64_t p(0);
        std::vector<std::string> v;
        v.emplace_back(
                currentRecord.seq.substr(0, k),
                buildCtx('N', currentRecord.seq[0+k+1])
        );
        for(int i = 1; i <= currentRecord.seq.size() - k - 1; ++i)
            v.emplace_back(
                    currentRecord.seq.substr(i, k),
                    buildCtx(currentRecord.seq[i-1], currentRecord.seq[i+k+1])
            );
        v.emplace_back(currentRecord.seq.substr(currentRecord.seq.size() - k, k),
                       buildCtx(currentRecord.seq[currentRecord.seq.size() - k - 1], 'N'));
        std::copy_if(v.begin(),v.end(), v.begin(), [](const std::string &s) { return s.find('N') == std::string::npos;});
        std::transform(v.begin(),v.end(), v.begin(), [&](std::string s) { return s < reverse_complement(s) ? s : reverse_complement(s);});
        return false;
    }

private:
    int buildCtx(char b, char f) const {
        int res=0;
        int bit = 0;
        if (b == 'A') res |= 1<<0;
        if (b == 'C') res |= 1<<1;
        if (b == 'T') res |= 1<<2;
        if (b == 'G') res |= 1<<3;
        if (f == 'A') res |= 1<<4;
        if (f == 'C') res |= 1<<5;
        if (f == 'T') res |= 1<<6;
        if (f == 'G') res |= 1<<7;
        return res;
    }
    std::string reverse_complement(const std::string &s) {
        std::string res;
        std::transform(s.rbegin(), s.rend(), res.begin(), [](const char &c) {
            if (c == 'A') return 'T';
            if (c == 'C') return 'G';
            if (c == 'T') return 'A';
            if (c == 'G') return 'C';
        });
        return res;
    }

    FileRecord currentRecord;

    unsigned int k = 31;
    std::string fkmer="";
    std::string rkmer="";
    unsigned int last_unknown=0;
    unsigned int pre_post=0;

};

#endif //BSG_KMERNEIGHBOUR_H
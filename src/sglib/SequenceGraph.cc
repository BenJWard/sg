//
// Created by Bernardo Clavijo (EI) on 18/10/2017.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cmath>
#include "SequenceGraph.h"
#include <sglib/mappers/LinkedReadMapper.hpp>
#include <list>
#include <queue>
#include <stack>
#include <tuple>
#include "sglib/readers/FileReader.h"

bool Node::is_canonical() {
    for (size_t i=0,j=sequence.size()-1;i<j;++i,--j){
        char f=sequence[i];
        char r=sequence[j];
        switch(r){
            case 'A':
                r='T';
                break;
            case 'C':
                r='G';
                break;
            case 'G':
                r='C';
                break;
            case 'T':
                r='A';
                break;
            case 'N':
                break;
            default:
            std::cout<<"unexpected character in fasta file: '"<<r<<"'"<<std::endl;
        }
        if (f<r) return true;
        if (r<f) return false;
    }
    return true;
};

void Node::make_rc() {
    std::string rseq;
    rseq.resize(sequence.size());
    for (size_t i=0,j=sequence.size()-1;i<sequence.size();++i,--j){
        switch(sequence[j]){
            case 'A':
                rseq[i]='T';
                break;
            case 'C':
                rseq[i]='G';
                break;
            case 'G':
                rseq[i]='C';
                break;
            case 'T':
                rseq[i]='A';
                break;
        }
    }
    std::swap(sequence,rseq);
};

sgNodeID_t SequenceGraph::add_node(Node n) {
    nodes.emplace_back(n);
    links.emplace_back();
    return (sgNodeID_t) nodes.size()-1;
}

void SequenceGraph::remove_node(sgNodeID_t n) {
    sgNodeID_t node=(n>0? n:-n);
    auto oldlinks=links[node];//this creates a copy to allow the iteration
    for (auto &l:oldlinks) remove_link(l.source,l.dest);
    nodes[node].status=sgNodeDeleted;
    //TODO: this is a lazy solution
    nodes[node].sequence.clear();
    //TODO: remove read mappings
}

void SequenceGraph::add_link(sgNodeID_t source, sgNodeID_t dest, int32_t d) {
    Link l(source,dest,d);
    links[(source > 0 ? source : -source)].emplace_back(l);
    std::swap(l.source,l.dest);
    links[(dest > 0 ? dest : -dest)].emplace_back(l);
}

void SequenceGraph::remove_link(sgNodeID_t source, sgNodeID_t dest) {
    auto & slinks = links[(source > 0 ? source : -source)];
    slinks.erase(std::remove(slinks.begin(), slinks.end(), Link(source,dest,0)), slinks.end());
    auto & dlinks = links[(dest > 0 ? dest : -dest)];
    dlinks.erase(std::remove(dlinks.begin(), dlinks.end(), Link(dest,source,0)), dlinks.end());

}

std::vector<Link> SequenceGraph::get_fw_links( sgNodeID_t n){
    std::vector<Link> r;
    for (auto &l:links[(n>0 ? n : -n)]) if (l.source==-n) r.emplace_back(l);
    return r;
}

std::vector<Link> SequenceGraph::get_bw_links(sgNodeID_t n) {
    return get_fw_links (-n);
}


bool Link::operator==(const Link a){
    if (a.source == this->source && a.dest == this->dest){
        return true;
    }
    return false;
}

bool Link::operator<(const Link a) const {
    if (a.source < this->source){
        return true;
    } if (a.source == this->source && a.dest < this->dest) {
        return  true;
    }
    return false;
}

std::vector<std::vector<sgNodeID_t >> SequenceGraph::find_bubbles(std::vector<sgNodeID_t> component){
    std::vector<std::vector<sgNodeID_t >> bubbles;
    std::vector<sgNodeID_t > checked;
    // loop over all links in component, if 2 (or x) nodes have same source and dest, they are bubbles
    for (auto n: component){
        //std::cout << "n " << "name: " << oldnames[n] << std::endl;
        std::vector<sgNodeID_t > bubble;
        bubble.push_back(n);
        // if we haven't checked this node
        if (std::find(checked.begin(), checked.end(), n) == checked.end()){
            auto links_n = links[n];
            std::set<Link> links_set;
            std::vector<Link> links_uniq;

            for (auto link:links_n){
                auto s = link.source > 0 ? link.source:-link.source;
                auto d = link.dest > 0 ? link.dest:-link.dest;

                //std::cout << "source: " << oldnames[s] << " dest: " << oldnames[d] << std::endl;
                if (links_set.find(link) == links_set.end()){
                    //std::cout << "not duplicate " <<std::endl;
                    links_uniq.push_back(link);
                    links_set.insert(link);
                }
            }
            /*for (auto u:  links_uniq){
                //std::cout << oldnames[u.source] << " " << oldnames[u.dest] << " \n";
                std::cout << u.source << " " << u.dest << " \n";
            }
            std::cout << "\nlinks unique: " << links_uniq.size() << " links set " << links_set.size() << std::endl;*/
            // if l has exactly 2 links, one source and one dest, it may be a bubble
            if (links_uniq.size() == 2) {
                std::vector<sgNodeID_t> linked_to;
                std::map<sgNodeID_t, int> linked_2nd_degree;
                for (auto l:links_uniq){
                    if (l.dest == n or l.dest == -n){
                        auto s = l.source > 0 ? l.source:-l.source;
                        linked_to.push_back(s);
                        std::set<sgNodeID_t > seen_s2;
                        std::set<sgNodeID_t > seen_d;
                        for (auto l_2:links[s]){
                            auto s2 = l_2.source > 0 ? l_2.source:-l_2.source;
                            if (seen_s2.find(s2) == seen_s2.end()) {
                                if (s2 != s && s2 != n) {
                                    linked_2nd_degree[s2] += 1;
                                }
                                seen_s2.insert(s2);
                            }

                            auto d = l_2.dest > 0 ? l_2.dest:-l_2.dest;
                            if (seen_d.find(d) == seen_d.end()) {
                                if (d != s && d != n) {

                                    linked_2nd_degree[d] += 1;
                                }
                                seen_d.insert(d);
                            }
                        }
                    } else if (l.source == n or l.source == -n) {
                        auto s = l.dest > 0 ? l.dest:-l.dest;
                        linked_to.push_back(s);

                        std::set<sgNodeID_t > seen_s2;
                        std::set<sgNodeID_t > seen_d;
                        for (auto l_2:links[s]){
                            auto s2 = l_2.source > 0 ? l_2.source:-l_2.source;
                            if (seen_s2.find(s2) == seen_s2.end()) {

                                if (s2 !=s && s2 !=n){
                                linked_2nd_degree[s2] += 1;

                            }
                                seen_s2.insert(s2);
                            }

                             auto d = l_2.dest > 0 ? l_2.dest:-l_2.dest;
                            if (seen_d.find(d) == seen_d.end()) {
                                if (d != s && d != n) {
                                linked_2nd_degree[d] += 1;
                                //checked.push_back(d);
                            }
                                seen_d.insert(d);
                            }
                        }
                    }

                }

/*for (auto link2:linked_2nd_degree){

    auto s = link2.first > 0 ? link2.first:-link2.first;
    std::cout << "2nd degree l: " << oldnames[s] << " " << link2.second << " l first: " << link2.first << std::endl;
}*/
                    // if n is a bubble, other bubble contigs will share source and dest
                    // nope... but each bubble contig should occur twice...


                    // if an element is 2nd degree joined to n twice, and linked to the same nodes as n, its in a bubble
                    for (auto j:linked_2nd_degree){
                        if (j.second == 2){
                            checked.push_back(j.first);
                            auto links_j = links[j.first];
                            std::set<sgNodeID_t > joined_j;
                            for (auto l_j:links_j){
                                if (l_j.source == j.first or l_j.source == -j.first){
                                    auto l_j_abs = l_j.dest > 0? l_j.dest: -l_j.dest;
                                    joined_j.insert(l_j_abs);
                                } else if (l_j.dest == j.first or l_j.dest == -j.first){
                                    auto l_j_abs = l_j.source > 0? l_j.source: -l_j.source;

                                    joined_j.insert(l_j_abs);

                                }
                            }

                            auto s = j.first > 0 ? j.first:-j.first;
                            /*std::cout << "j: " << oldnames[s] << " id: " << j.first;
                            for (auto jhg:joined_j){
                                auto d = jhg > 0 ? jhg:-jhg;

                                std::cout << " "<< oldnames[d] << " ";
                            }
                            std::cout << std::endl;*/
                            if (joined_j.size() == linked_to.size()){
                                std::vector<sgNodeID_t > joined_vec;
                                for (auto joined:joined_j){
                                    joined_vec.push_back(joined);
                                    //std::cout << joined << " ";
                                }
                                //std::cout <<std::endl;
                                for (auto joined: linked_to){
                                    joined_vec.push_back(joined);
                                    //std::cout << joined << " ";
                                }
                                //std::cout << std::endl;

                                if (std::is_permutation(linked_to.begin(), linked_to.end(), joined_vec.begin()) ) {
                                    bubble.push_back(j.first);
                                }
                            }
                        }

                    }


            }

        }
        if (bubble.size() > 1){
            bubbles.push_back(bubble);
        }
    }
    std::cout << "Found " << bubbles.size() << " bubbles in component of " << component.size() << " nodes " <<std::endl;
    return bubbles;
}


std::vector<std::vector<sgNodeID_t>> SequenceGraph::connected_components(int max_nr_totalinks, int max_nr_dirlinks,
                                                                         int min_rsize) {
    std::vector<bool> used(nodes.size());
    std::vector<std::vector<sgNodeID_t>> components;
    //TODO: first find all repeats, add them as independent components and mark them as used.
    size_t max_component = 0;
    for (sgNodeID_t start_node=1;start_node<nodes.size();++start_node){
        if (false==used[start_node]){
            used[start_node]=true;
            //if start node is repeat, just add a single-node component
            std::set<sgNodeID_t> in_component;
            in_component.insert(start_node);
            std::set<sgNodeID_t> to_explore;
            to_explore.insert(start_node);
            while(to_explore.size()){
                auto n=*to_explore.begin();
                to_explore.erase(n);
                //find all connections that are not already used, include both in component and in to_explore
                for (auto l:links[n]) {
                    sgNodeID_t next=(l.dest>0 ? l.dest : -l.dest);
                    if (in_component.count(next)==0) {
                        in_component.insert(next);
                        to_explore.insert(next);
                    }
                }
            }
            components.emplace_back();
            for (sgNodeID_t n:in_component){
                components.back().push_back(n);
                used[n]=true;
            }
            if (in_component.size() > max_component){
                max_component = in_component.size();
            }
        }
    }
    std::cout << "Max component size: " << max_component << std::endl;
    return components;
}

void SequenceGraph::write(std::ofstream & output_file) {
    uint64_t count;
    count=nodes.size();
    output_file.write((char *) &count,sizeof(count));
    for (auto &n:nodes){
        output_file.write((char *) &n.status,sizeof(n.status));
        count=n.sequence.size();
        output_file.write((char *) &count,sizeof(count));
        output_file.write((char *) n.sequence.c_str(),count);
    }
    count=links.size();
    output_file.write((char *) &count,sizeof(count));
    for (auto nl:links) {
        uint64_t lcount=nl.size();
        output_file.write((char *) &lcount,sizeof(lcount));
        output_file.write((char *) nl.data(), sizeof(Link) * lcount);
    }
}

void SequenceGraph::read(std::ifstream & input_file) {
    uint64_t count;
    input_file.read((char *) &count,sizeof(count));
    nodes.clear();
    nodes.reserve(count);
    uint64_t active=0;
    for (auto i=0;i<count;++i){
        uint64_t seqsize;
        std::string seq;
        sgNodeStatus_t status;
        input_file.read((char *) &status,sizeof(status));
        input_file.read((char *) &seqsize,sizeof(seqsize));
        seq.resize(seqsize);
        input_file.read((char *) seq.c_str(),seqsize);
        nodes.emplace_back(seq);
        nodes.back().status=status;
        if (nodes.back().status==sgNodeStatus_t::sgNodeActive) ++active;
    }
    input_file.read((char *) &count,sizeof(count));
    links.resize(count);
    for (auto &nl:links) {
        uint64_t lcount;
        input_file.read((char *) &lcount,sizeof(lcount));
        nl.resize(lcount,{0,0,0});
        input_file.read((char *) nl.data(), sizeof(Link) * lcount);
    }
}

void SequenceGraph::load_from_gfa(std::string filename) {
    std::string line;
    this->filename=filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        this->fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ifstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't read gfa file");
    std::getline(gfaf, line);
    if (line!="H\tVN:Z:1.0") std::cout<<"WARNING, first line of gfa doesn't correspond to GFA1"<<std::endl;

    std::ifstream fastaf(fasta_filename);
    std::cout << "fasta filesname: " << fasta_filename << std::endl;
    if (!fastaf) throw std::invalid_argument("Can't read fasta file");



    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Loading sequences from "<<fasta_filename<<std::endl;

    std::string name,seq="";
    seq.reserve(10000000); //stupid hack but probably useful to reserve
    oldnames_to_ids.clear();
    oldnames.push_back("");
    nodes.clear();
    links.clear();
    add_node(Node("",sgNodeStatus_t::sgNodeDeleted)); //an empty deleted node on 0, just to skip the space
    sgNodeID_t nextid=1;
    uint64_t rcnodes=0;
    while(!fastaf.eof()){
        std::getline(fastaf,line);
        if (fastaf.eof() or line[0]=='>'){
            if (!name.empty()) {
                //rough ansi C and C++ mix but it works
                if (oldnames_to_ids.find(name) != oldnames_to_ids.end())
                    throw std::logic_error("sequence " + name + " is already defined");
                oldnames_to_ids[name] = add_node(Node(seq));
                oldnames.push_back(name);
                //reverse seq if not canonical
                if (!nodes.back().is_canonical()) {
                    nodes.back().make_rc();
                    oldnames_to_ids[name] = -oldnames_to_ids[name];
                    ++rcnodes;
                }
            }
            name.clear();
            for (auto i = 1; i < line.size() and line[i] != ' '; ++i) name += line[i];
            seq = "";
        } else {
            seq+=line;
        }
    }
    std::cout<<nodes.size()-1<<" nodes loaded! "<<rcnodes<<" canonised"<<std::endl;

    //load store all conections.

    std::string gfa_rtype,gfa_source,gfa_sourcedir,gfa_dest,gfa_destdir,gfa_cigar;
    sgNodeID_t src_id,dest_id;
    int32_t dist;
    uint64_t lcount=0;
    uint64_t dist_egt0(0);
    while(std::getline(gfaf, line) and !gfaf.eof()) {
        std::istringstream iss(line);
        iss >> gfa_rtype;
        if (gfa_rtype == "L"){
            iss >> gfa_source;
            iss >> gfa_sourcedir;
            iss >> gfa_dest;
            iss >> gfa_destdir;
            iss >> gfa_cigar;
            //std::cout<<"'"<<source<<"' '"<<gfa_sourcedir<<"' '"<<dest<<"' '"<<destdir<<"'"<<std::endl;
            if (oldnames_to_ids.find(gfa_source)==oldnames_to_ids.end()){
                oldnames_to_ids[gfa_source] = add_node(Node(""));
                //std::cout<<"added source!" <<source<<std::endl;
            }
            if (oldnames_to_ids.find(gfa_dest)==oldnames_to_ids.end()){
                oldnames_to_ids[gfa_dest] = add_node(Node(""));
                //std::cout<<"added dest! "<<dest<<std::endl;
            }
            src_id=oldnames_to_ids[gfa_source];
            dest_id=oldnames_to_ids[gfa_dest];

            //Go from GFA's "Links as paths" to a normal "nodes have sinks (-/start/left) and sources (+/end/right)
            if (gfa_sourcedir=="+") src_id=-src_id;
            if (gfa_destdir=="-") dest_id=-dest_id;
            dist=0;
            if (gfa_cigar.size()>1 and gfa_cigar[gfa_cigar.size()-1]=='M') {
                //TODO: better checks for M
                dist=-atoi(gfa_cigar.c_str());
                if (dist>=0) {
                    dist_egt0++;
                }
            }
            add_link(src_id,dest_id,dist);
            ++lcount;
        }
    }
    if (dist_egt0 > lcount*0.5f) {
        std::cout << "Warning: The loaded graph contains " << dist_egt0 << " non-overlapping links out of " << lcount << std::endl;
    }
    std::cout<<nodes.size()-1<<" nodes after connecting with "<<lcount<<" links"<<std::endl;
}

void SequenceGraph::write_to_gfa(std::string filename, const std::unordered_set<sgNodeID_t> & mark_red, const std::vector<double> & depths, const std::unordered_set<sgNodeID_t> & selected_nodes){
    std::string fasta_filename;
    //check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ofstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't write to gfa file");
    gfaf<<"H\tVN:Z:1.0"<<std::endl;

    std::ofstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't write to fasta file");


    //load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Writing sequences to "<<fasta_filename<<std::endl;

    for (sgNodeID_t i=1;i<nodes.size();++i){
        if (nodes[i].status==sgNodeDeleted) continue;
        if (!selected_nodes.empty() and selected_nodes.count(i)==0 and selected_nodes.count(-i)==0) continue;
        fastaf<<">seq"<<i<<std::endl<<nodes[i].sequence<<std::endl;
        gfaf<<"S\tseq"<<i<<"\t*\tLN:i:"<<nodes[i].sequence.size()<<"\tUR:Z:"<<fasta_filename
                <<(mark_red.count(i)?"\tCL:Z:red":"")<<(depths.empty() or isnan(depths[i])?"":"\tDP:f:"+std::to_string(depths[i]))<<std::endl;
    }

    for (auto &ls:links){
        for (auto &l:ls)
            if (l.source<=l.dest and (selected_nodes.empty() or
                    selected_nodes.count(l.source)>0 or selected_nodes.count(-l.source)>0 or
                            selected_nodes.count(l.dest)>0 or selected_nodes.count(-l.dest)>0)) {
                gfaf<<"L\t";
                if (l.source>0) gfaf<<"seq"<<l.source<<"\t-\t";
                else gfaf<<"seq"<<-l.source<<"\t+\t";
                if (l.dest>0) gfaf<<"seq"<<l.dest<<"\t+\t";
                else gfaf<<"seq"<<-l.dest<<"\t-\t";
                gfaf<<(l.dist<0 ? -l.dist : 0)<<"M"<<std::endl;
            }
    }

}

std::vector<sgNodeID_t> SequenceGraph::oldnames_to_nodes(std::string _oldnames) {
    std::vector<sgNodeID_t> nv;
    const char * s = _oldnames.c_str();
    const char * sign;
    std::string oldname;
    sgNodeID_t node;
    while (*s!= NULL){
        for (sign=s;*sign!=NULL and *sign !='+' and *sign!='-';++sign);
        if (sign ==s or *sign==NULL ) throw std::invalid_argument("invalid path specification");
        oldname=s;
        oldname.resize(sign-s);
        node=oldnames_to_ids[oldname];
        if (0==node) throw std::invalid_argument("node "+oldname+" doesn't exist in graph");
        if (*sign=='-') node=-node;
        nv.push_back(node);
        s=sign+1;
        if (*s==',') ++s;
        else if (*s!=NULL) throw std::invalid_argument("invalid path specification");

    }
    return nv;
}

std::string SequenceGraphPath::get_fasta_header() {
    std::string h=">sgPath_";
    for (auto &n:nodes) {
        h += std::to_string(n)+",";
    }
    h.resize(h.size()-1);
    return h;
}

std::string SequenceGraphPath::get_sequence() {
    std::string s="";
    sgNodeID_t pnode=0;
    // just iterate over every node in path - contig names are converted to ids at construction
    for (auto &n:nodes) {
        std::string nseq;
        if (n>0){
            nseq=sg.nodes[n].sequence;
        } else {
            auto rcn=sg.nodes[-n];
            rcn.make_rc();
            nseq=rcn.sequence;
        }
        if (pnode !=0){
            //find link between pnode' output (+pnode) and n's sink (-n)
            auto l=sg.links[(pnode>0 ? pnode:-pnode)].begin();
            for (;l!=sg.links[(pnode>0 ? pnode:-pnode)].end();++l)
                if (l->source==pnode and l->dest==n) break;
            if (l==sg.links[(pnode>0 ? pnode:-pnode)].end()) {
                std::cout<<"can't find a link between "<<pnode<<" and "<<n<<std::endl;
                throw std::runtime_error("path has no link");
            } else {
                if (l->dist>0){
                    for (auto c=l->dist;c>0;--c) s+="N";
                }
                else {
                    auto ovl=-l->dist;
                    for (auto s1=s.c_str()+s.size()-ovl,s2=nseq.c_str();*s1!=NULL;++s1,++s2)
                        if (*s1!=*s2)
                            throw std::runtime_error("path overlap is invalid!");
                    nseq.erase(0,ovl);
                }
            }
        }
        s+=nseq;
        pnode=-n;
    }
    return s;
}

std::vector<SequenceGraphPath> SequenceGraph::get_all_unitigs(uint16_t min_nodes) {
    std::vector<SequenceGraphPath> unitigs;
    std::vector<bool> used(nodes.size(),false);

    for (auto n=1;n<nodes.size();++n){
        if (used[n] or nodes[n].status==sgNodeDeleted) continue;
        used[n]=true;
        SequenceGraphPath path(*this,{n});

        //two passes: 0->fw, 1->bw, path is inverted twice, so still n is +
        for (auto pass=0; pass<2; ++pass) {
            //walk til a "non-unitig" junction
            for (auto fn = get_fw_links(path.nodes.back()); fn.size() == 1; fn = get_fw_links(path.nodes.back())) {
                if (fn[0].dest != n and fn[0].dest != -n and get_bw_links(fn[0].dest).size() == 1) {
                    path.nodes.emplace_back(fn[0].dest);
                    used[fn[0].dest > 0 ? fn[0].dest : -fn[0].dest] = true;
                } else break;
            }
            path.reverse();
        }
        if (path.nodes.size()>=min_nodes) unitigs.push_back(path);
    }
    return unitigs;
}

void SequenceGraph::join_all_unitigs() {
    for (auto p:get_all_unitigs(2)){
        std::cout << "Unitig found: ";
        for (auto n:p.nodes) std::cout<< n<< " ";
        std::cout<<std::endl;
        join_path(p);
    }
}

void SequenceGraph::join_path(SequenceGraphPath p, bool consume_nodes) {
    std::set<sgNodeID_t> pnodes;
    for (auto n:p.nodes) {
        pnodes.insert( n );
        pnodes.insert( -n );
    }

    if (!p.is_canonical()) p.reverse();
    sgNodeID_t new_node=add_node(Node(p.get_sequence()));
    //TODO:check, this may have a problem with a circle
    for (auto l:get_bw_links(p.nodes.front())) add_link(new_node,l.dest,l.dist);
    for (auto l:get_fw_links(p.nodes.back())) add_link(-new_node,l.dest,l.dist);

    //TODO: update read mappings
    if (consume_nodes) {
        for (auto n:p.nodes) {
            //check if the node has neighbours not included in the path.
            bool ext_neigh=false;
            if (n!=p.nodes.back()) for (auto l:get_fw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (n!=p.nodes.front()) for (auto l:get_bw_links(n)) if (pnodes.count(l.dest)==0) ext_neigh=true;
            if (ext_neigh) continue;
            remove_node(n);
        }
    }
}

void SequenceGraphPath::reverse(){
    std::vector<sgNodeID_t> newn;
    for (auto n=nodes.rbegin();n<nodes.rend();++n) newn.emplace_back(-*n);
    //std::swap(nodes,newn);
    nodes=newn;
}

bool SequenceGraphPath::is_canonical() {
    auto rp=*this;
    rp.reverse();
    return this->get_sequence()<rp.get_sequence();
}

const bool SequenceGraphPath::operator<(const SequenceGraphPath &other) const {
    for (auto i=0;i<nodes.size();++i){
        if (other.nodes.size()<i) return true;
        if (nodes[i]<other.nodes[i]) return true;
        if (nodes[i]>other.nodes[i]) return false;
    }
    return false;
}

const bool SequenceGraphPath::operator==(const SequenceGraphPath &other) const {
    if (other.nodes.size()!=nodes.size()) return false;
    for (auto i=0;i<nodes.size();++i){
        if (nodes[i]!=other.nodes[i]) return false;
    }
    return true;
}

bool SequenceGraphPath::extend_if_coherent(SequenceGraphPath s) {
//    int offset=-1;
//    for (auto i=0;i<nodes.size();++i) {
//        if (nodes[i] == s.nodes[0]) {
//            offset = i;
//            break;
//        }
//    }
//    if (offset<=0) {
//    if (offset<0) {
//        offset=1;
//        for (auto i=0;i<nodes.size();++i) {
//            if (nodes[i] == s.nodes[0]) {
//                offset = -i;
//                break;
//            }
//        }
//    }
}

std::vector<SequenceSubGraph> SequenceGraph::get_all_tribbles() {


    for (sgNodeID_t n=1;n<nodes.size();++n) {
        //Heuristic to find "tribbles"
        // A --- B -- C -- H
        //  \     \      /
        //   \     E    /
        //    \     \  /
        //      F -- G
        // A->[B-F]
        // B->[C-E]
        // C->H
        // E->G
        // F->G
        // F->H
        auto a_fw=get_fw_links(n);
        if (a_fw.size()!=2) continue;
        auto b_fw=get_fw_links(a_fw[0].dist);

        sgNodeID_t A, B, C, D, E, F, G, H;
    }


}

std::vector<sgNodeID_t > SequenceGraph::find_canonical_repeats() {
    std::vector<sgNodeID_t > repeaty_nodes;

    uint64_t count=0, l700=0,l2000=0,l4000=0,l10000=0,big=0,checked=0,solvable=0;

    for (sgNodeID_t n=1; n < nodes.size(); ++n) {
        auto nfw_links = get_fw_links(n).size();
        auto nbw_links = get_bw_links(n).size();
        if ( nfw_links == nbw_links and nfw_links==2){
            ++count;

            if (nodes[n].sequence.size() < 700) ++l700;
            else if (nodes[n].sequence.size() < 2000) ++l2000;
            else if (nodes[n].sequence.size() < 4000) ++l4000;
            else if (nodes[n].sequence.size() < 10000) ++l10000;
            else ++big;

            if (nodes[n].sequence.size() > 1000) {
                // std::cout << "evaluating trivial repeat at " << n << "(" << sg.nodes[n].sequence.size() << "bp)" << std::endl;
                repeaty_nodes.push_back(n);
            }
        }
    }

    std::cout << "Candidates for canonical repeat expansion:                    " << count << std::endl;
    std::cout << "Candidates for canonical repeat expansion <700bp:             " << l700 << std::endl;
    std::cout << "Candidates for canonical repeat expansion >700bp & <2000bp:   " << l2000 << std::endl;
    std::cout << "Candidates for canonical repeat expansion >2000bp & <4000bp:  " << l4000 << std::endl;
    std::cout << "Candidates for canonical repeat expansion >4000bp & <10000bp: " << l10000 << std::endl;
    std::cout << "Candidates for canonical repeat expansion >10000bp:           " << big << std::endl;
    std::cout << "Trivially solvable canonical repeats:                         " << solvable << "/" << checked << std::endl;

    return repeaty_nodes;
}

std::vector<sgNodeID_t>
SequenceGraph::breath_first_search(std::vector<sgNodeID_t> &nodes, unsigned int size_limit) {
    std::queue<sgNodeID_t> to_visit(std::deque<sgNodeID_t>(nodes.begin(),nodes.end()));
    std::set<sgNodeID_t> visited;
    std::unordered_map<sgNodeID_t, sgNodeID_t > meta;

    while (!to_visit.empty() and visited.size() < size_limit) {
        const auto activeNode(to_visit.front());
        to_visit.pop();
        for (const auto &neighboor: get_fw_links(activeNode)) {
            if (visited.find(neighboor.dest) == visited.end()) {
                to_visit.push(neighboor.dest);
                visited.insert(neighboor.dest);
            }
        }
    }
    return std::vector<sgNodeID_t>(visited.begin(), visited.end());
}

std::vector<nodeVisitor>
SequenceGraph::depth_first_search(const nodeVisitor seed, unsigned int size_limit, unsigned int edge_limit, std::set<nodeVisitor> tabu) {
    // Create a stack with the nodes and the path length

    std::stack<nodeVisitor> to_visit;
    to_visit.push(seed);
    std::set<nodeVisitor> visited_set;
    for (const auto &n:tabu) {
        visited_set.insert(n);
    }

    while (!to_visit.empty()) {
        const auto activeNode(to_visit.top());
        to_visit.pop();
        auto looker = nodeVisitor(activeNode.node,0,0);
        auto visitedNode(visited_set.find( looker ) );
        if (visitedNode == visited_set.end() and
                (activeNode.path_length < edge_limit or edge_limit==0) and  // Is within path limits
                (activeNode.dist < size_limit or size_limit==0) ) // Is within sequence limits
        {
            visited_set.insert(nodeVisitor(activeNode.node, activeNode.dist, activeNode.path_length));
            for (const auto &l: get_fw_links(activeNode.node)) {
                to_visit.push(nodeVisitor(l.dest,
                                 static_cast<uint>(activeNode.dist + nodes[l.dest > 0 ? l.dest : -l.dest].sequence.length()),
                                 activeNode.path_length+1));
            }
        } else if (visitedNode != visited_set.end() and // If the node has been visited
                (visitedNode->path_length > activeNode.path_length or visitedNode->dist > activeNode.dist)) // but is closer than the last time I saw it
        {
            visited_set.erase(visitedNode); // Remove it from the visited nodes and
            to_visit.push(nodeVisitor(activeNode.node, activeNode.dist, activeNode.path_length)); // add it to the list of nodes to visit
        }
    }

    std::vector<nodeVisitor > result;
    for (const auto &v:visited_set) {
        if (v.node == seed.node) continue;
        result.push_back(v);
    }
    return result;
}

void SequenceSubGraph::write_to_gfa(std::string filename) {
    std::string fasta_filename;
//check the filename ends in .gfa
    if (filename.size()>4 and filename.substr(filename.size()-4,4)==".gfa"){
        fasta_filename=filename.substr(0,filename.size()-4)+".fasta";
    }
    else throw std::invalid_argument("filename of the gfa input does not end in gfa, it ends in '"+filename.substr(filename.size()-4,4)+"'");

    std::ofstream gfaf(filename);
    if (!gfaf) throw std::invalid_argument("Can't write to gfa file");
    gfaf<<"H\tVN:Z:1.0"<<std::endl;

    std::ofstream fastaf(fasta_filename);
    if (!fastaf) throw std::invalid_argument("Can't write to fasta file");


//load all sequences from fasta file if they're not canonical, flip and remember they're flipped
    std::cout<<"Writing sequences to "<<fasta_filename<<std::endl;
    std::unordered_set<sgNodeID_t > nodes_in_links;
    for (const auto &n:nodes) {
        for (auto &l:sg.links[std::abs(n)]) {
            nodes_in_links.insert(std::abs(l.source));
            nodes_in_links.insert(std::abs(l.dest));
        }
    }

    for (const auto &n:nodes_in_links){
        fastaf<<">"<<sg.oldnames[n]<<std::endl<<sg.nodes[std::abs(n)].sequence<<std::endl;
        gfaf<<"S\t"<<sg.oldnames[n]<<"\t*\tLN:i:"<<sg.nodes[std::abs(n)].sequence.length()<<"\tUR:Z:"<<fasta_filename<<std::endl;
    }

    for (const auto &n:nodes) {
        for (auto &l:sg.links[std::abs(n)]) {
            gfaf << "L\t";
            if (l.source > 0) gfaf << "" << sg.oldnames[std::abs(l.source)] << "\t-\t";
            else gfaf << "" << sg.oldnames[std::abs(l.source)] << "\t+\t";
            if (l.dest > 0) gfaf << "" << sg.oldnames[std::abs(l.dest)] << "\t+\t";
            else gfaf << "" << sg.oldnames[std::abs(l.dest)] << "\t-\t";
            gfaf << (l.dist < 0 ? -l.dist : 0) << "M" << std::endl;
        }
    }
}

    return std::vector<sgNodeID_t>(visited.begin(), visited.end());
}

SequenceGraph::SequenceGraph(const std::vector<std::string> &seqs, unsigned int k, unsigned int min_cov) {
    /*
     * Get context augmented k-mers from the sequences filtered by min_cov
     */

    uint countK(k-2);
    // K = k-2 as the kmers are augmented with a bwd and fwd context -> K+bwd+fwd = inputK
    SMR<KmerNeighbour,
            kmerNeighbourFactory<FastaRecord>,
            StringReader<FastaRecord>,
            FastaRecord, StringReaderParams, KMerNeighbourFactoryParams> kmerNeigh_SMR({seqs,1}, {countK}, {4*GB, min_cov, 10000, "default"});

    std::vector<KmerNeighbour> kmers = kmerNeigh_SMR.process_from_memory();
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

    std::vector<KmerNeighbour> node_neighs;
    for (auto &km: kmers) {
        KmerNeighbour &current = km;
        std::deque<char> node(km.kmer.begin(), km.kmer.end());
        // If can go fwd
        KmerNeighbour neighbour;
        if (current.canExtendFwd()) {
            // Can go bwd? leave
            if (current.canExtendBwd()) break;
            do {
                neighbour = current.extendFwd(countK, node, kmerDict);
                current = neighbour;
            } while (current.canExtendFwd());
        }
        else if (current.canExtendBwd()) {
            do {
                neighbour = current.extendBwd(countK, node, kmerDict);
                current = neighbour;
            } while (current.canExtendBwd());
        }
        // Build node
        Node t(std::string(node.cbegin(), node.cend()));
        t.make_rc();
        bool isCanonical(t.is_canonical());
        node_neighs.emplace_back(isCanonical ?
                                 KmerNeighbour(t.sequence,current.pre_post_context) :
                                 KmerNeighbour(std::string(node.cbegin(),node.cend()), current.pre_post_context)
        );
    }

    std::unordered_map<std::string, KmerNeighbour>().swap(kmerDict);
    // Create links for the nodes
    std::unordered_map<KmerNeighbour, std::vector<sgNodeID_t > > tip_to_node;
    for (std::vector<KmerNeighbour>::iterator node = node_neighs.begin(); node!=node_neighs.end(); ++node){

    }
}
#ifndef GENERES_H
#define GENERES_H

#include <string>

#include "helper.h"

class GeneRes
{
public:
    GeneRes();
    GeneRes(std::string _id) { id = _id; }

    GeneRes(std::string _id, std::string _chr, int _start, int _end, std::string _name, std::string _desc);
    GeneRes(std::string _id, std::string _gene_loc, std::string _gene_desc);

    void set_data(std::string _id, int _chr, int _start, int _end, std::string _name, std::string _desc);
    void set_data(std::string _id, std::string _gene_loc, std::string _gene_desc);

    std::string get_id() { return id; }
    std::string get_name() { return name; }
    std::string get_loc() { return chr + ":" + int2str(start) + ".." + int2str(end); }
    std::string get_desc() { return desc; }


    bool operator==(const GeneRes & rhs) const {
        if( id.compare(rhs.id)==0 ) {
            return true;
        }
        else {
            return false;
        }
    }

    std::string get_gene_details() {
        std::string details = "Gene ID : " + id + "\nGene Symbol : " + name + "\nDescription : " + desc + "\nLocus : " + get_loc();
        return details;
    }

    std::string id;

private:
    std::string chr;
    int start;
    int end;
    std::string name;
    std::string desc;
};

#endif // GENERES_H

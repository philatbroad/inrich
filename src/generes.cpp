#include "generes.h"
#include "helper.h"

#include <string>
#include <iostream>

GeneRes::GeneRes()
{
}


GeneRes::GeneRes(std::string _id, std::string _chr, int _start, int _end, std::string _name, std::string _desc) {
    id = _id;
    chr = _chr;
    start = _start;
    end = _end;
    name = _name;
    desc = _desc;
}


GeneRes::GeneRes(std::string _id, std::string _gene_loc, std::string _gene_desc) {
    set_data(_id, _gene_loc, _gene_desc);
}


void GeneRes::set_data(std::string _id, int _chr, int _start, int _end, std::string _name, std::string _desc) {
    id = _id;
    chr = _chr;
    start = _start;
    end = _end;
    name = _name;
    desc = _desc;
}

void GeneRes::set_data(std::string _id, std::string _gene_loc, std::string _gene_desc) {
    id = _id;

    // ID:CHR:START..END
    std::vector<std::string> fields = tokenize_string(_gene_loc, ':');
    if(fields.size()==3) {
        chr = fields.at(1);
        std::vector<std::string> loc = tokenize_string(fields.at(2), '.');
        str2int(loc.at(0), start);
        str2int(loc.at(loc.size()-1), end);
    }
    else if(fields.size()==2) {
        chr = fields.at(0);
        std::vector<std::string> loc = tokenize_string(fields.at(1), '.');
        str2int(loc.at(0), start);
        str2int(loc.at(loc.size()-1), end);
    }
    else {
        std::cerr << "ERROR: geneResult " << _gene_loc << "\n";
    }

    // NAME DESC
    fields = tokenize_string(_gene_desc, ' ');
    if(fields.size()>1) {
        name = fields.at(0);
        desc = "";
        for(unsigned int i=1; i<fields.size(); i++) {
            desc = desc + fields.at(i) + " ";
        }
    }
}

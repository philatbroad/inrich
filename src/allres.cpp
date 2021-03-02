#include "allres.h"
#include "helper.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

AllRes::AllRes()
{
}


AllRes::AllRes(std::string target_desc, std::string p, std::string gene_id, std::string gene_loc, std::string gene_desc, std::string interval_desc) {
    set_data(target_desc, p, gene_id, gene_loc, gene_desc, interval_desc);
}

AllRes::~AllRes()
{
}

PathwayRes* AllRes::get_pathway_result() {
    return new PathwayRes(pathway_id, pathway_name, pathway_p);
}

GeneRes* AllRes::get_gene_result() {
    return new GeneRes(gene_id, gene_chr, gene_start, gene_end, gene_symbol, gene_desc);
}

IntervalRes* AllRes::get_interval_result() {
    return new IntervalRes(interval_id, interval_loc);
}



void AllRes::set_data(std::string _target_desc, std::string _p, std::string _gene_id, std::string _gene_loc,
                         std::string _gene_desc, std::string _interval_desc) {


    // ==================================
    // Pathway Information
    // ==================================
    std::vector<std::string> fields = tokenize_string(_target_desc, ' ');
    pathway_id = fields.at(0);
    pathway_name = "";
    for(unsigned int i=1; i<fields.size(); i++) {
    	pathway_name = pathway_name + fields.at(i) + " ";
    }

    if(pathway_name.size()>0) {
        pathway_name = first2upper_str(pathway_name);
    }
    else {
        std::cerr << "ERROR: AllRes 45 " << _target_desc << "\n";
    }
    str2double(_p, pathway_p);


    // ==================================
    // Gene Information
    // ==================================
    gene_id = _gene_id;

    // ID:CHR:START..END

    fields = tokenize_string(_gene_loc, ':');
    if( fields.size()==3 ) {
        gene_chr = fields.at(1);
        std::vector<std::string> loc = tokenize_string(fields.at(2), '.');
        str2int(loc.at(0), gene_start);
        str2int(loc.at(loc.size()-1), gene_end);
    }
    else if( fields.size()==2 ) {
        gene_chr = fields.at(0);
        std::vector<std::string> loc = tokenize_string(fields.at(1), '.');
    }
    else {
        std::cerr << "ERROR: AllRes65" << _gene_loc << "\n";
    }

    // NAME DESC
    fields = tokenize_string(_gene_desc, ' ');
    gene_symbol = fields.at(0);
    gene_desc = "";
    for(unsigned int i=1; i<fields.size(); i++) {
    	gene_desc = gene_desc + fields.at(i) + " ";
    }
    if(gene_desc.size()>0) {
        gene_desc = first2upper_str(gene_desc);
    }


    // ==================================
    // Interval Information
    // ==================================
    // ID:CHR:START..END or ID:CHR:BP
    fields = tokenize_string(_interval_desc, ':');

    interval_id = fields.at(0);
    interval_loc = "";
    for(unsigned int i=1; i<fields.size(); i++) {
    	interval_loc = interval_loc + fields.at(i) + ":";
    }
    interval_loc = interval_loc.substr(0, interval_loc.size()-1);
}

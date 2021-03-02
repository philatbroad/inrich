#include "pathwayres.h"
#include "helper.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>


PathwayRes::PathwayRes(std::string _id) {
    id = _id;

    p = -1;
    pcorr = -1;
    exp_hit = -1;
    t_target = "na";
    name = "na";
    n_int = -1;
    low_warn = "na";
    unsh = "na";
    exp_gene_no = -1;
}

PathwayRes::PathwayRes(std::string _id, std::string _name, double _p) {
    id = _id;
    name = _name;
    p = _p;

    pcorr = -1;
    exp_hit = -1;
    t_target = "na";
    n_int = -1;
    low_warn = "na";
    unsh = "na";
    exp_gene_no = -1;
}


PathwayRes::PathwayRes(std::string _t_target, std::string _n_int, double _p, double _pcorr,
                             std::string _low_warning, std::string _unsh,
                             std::string _target_id, std::string _target_name,
                             double _exp_hit,
                             double _exp_gene_no) {
    t_target = _t_target;
    n_int = _n_int;
    p = _p;
    pcorr = _pcorr;

    if(_low_warning.compare("1")==0) {
        low_warn = "yes";
    }
    else {
        low_warn = "no";
    }
    unsh = _unsh;
    id = _target_id;

    if(_target_name.size()>0) {
        _target_name = first2upper_str(_target_name);
        name = _target_name;
    }
    exp_hit = _exp_hit;
    exp_gene_no = _exp_gene_no;
}


PathwayRes::PathwayRes(int _target_size, int _n_intervals, double _p, double _pcorr,
                             double _low_warning, double _unsh, std::string _id_name,
                             double _exp_hit,
                             double _exp_gene_no) {

    t_target = int2str(_target_size);
    n_int=int2str(_n_intervals);
    p = _p;
    pcorr = _pcorr;
    low_warn=double2str(_low_warning);
    unsh=double2str(_unsh);

    std::vector<std::string> tmp = tokenize_string(_id_name, ' ');
    id = tmp.at(0);
    name = "";
    for(unsigned int i=1; i<tmp.size(); i++) {
        name = name + tmp.at(i) + " ";
    }
    name = first2upper_str(name);
    exp_hit = _exp_hit;
    exp_gene_no = _exp_gene_no;
}

PathwayRes::PathwayRes(std::string _target_size, std::string _n_intervals, double _p, double _pcorr,
                             std::string _low_warning, std::string _unsh, std::string _id_name,
                             double _exp_hit,
                             double _exp_gene_no) {

    t_target = _target_size;
    n_int = _n_intervals;
    p = _p;
    pcorr = _pcorr;
    low_warn = _low_warning;
    unsh = _unsh;

    std::vector<std::string> tmp = tokenize_string(_id_name, ' ');
    id = tmp.at(0);
    name = "";
    for(unsigned int i=1; i<tmp.size(); i++) {
        name = name + tmp.at(i) + " ";
    }
    name = first2upper_str(name);
    exp_hit = _exp_hit;
    exp_gene_no = _exp_gene_no;
}


void PathwayRes::set_data(std::string _t_target, std::string _n_int,
                             double P, double _pcorr,
                             std::string _low_warning,
                             std::string _unsh,
                             std::string _target_id,
                             std::string _target_name,
                             double _exp_hit,
                             double _exp_gene_no) {
    t_target = _t_target;
    n_int = _n_int;
    p = P;
    pcorr = _pcorr;


    low_warn = _low_warning;

    unsh = _unsh;
    id = _target_id;
    if(_target_name.size()>0) {
        _target_name = first2upper_str(_target_name);
        name = _target_name;
    }
    else {
        name = "ERROR";
    }    

    exp_hit = _exp_hit;
    exp_gene_no = _exp_gene_no;
}

std::string PathwayRes::get_pathway_details() {
    std::string str_p;
    std::string str_pcorr;

    std::string details = "Target ID: " + id + "\n" +
                      "Target Term: " + name + "\n" +
                      "Total Genes in Target: " + t_target + "\n" +
                      "No of Significant Intervals/Genes in Target: " +  n_int  + "\n";

    if(exp_gene_no != -1 ) {
        details = details + "Expected No of Genes in Target: " + double2str(exp_gene_no) + "\n";
    }

    details = details + "Empirical P: " + double2str(p) + " ( " + double2str(pcorr) + " after correction )\n";

    if(exp_hit != -1) {
        details = details + "Expected Hits per Study: " + double2str(exp_hit) + "\n";
    }

    if(unsh.compare("")!=0 && unsh.compare("0")!=0 && unsh.compare("-1")!=0) {
        details = details + "** UnShuffled Proportion " + unsh + "!\n";
    }

    if(low_warn.compare("1")==0) {
        details = details + "** Low Permutation Warning!\n";
    }


    return details;
}

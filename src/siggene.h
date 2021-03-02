#ifndef SIGGENE_H
#define SIGGENE_H

#include <string>
#include "helper.h"

class SigGene {
public:
    SigGene() {; }
    SigGene(std::string _gene, std::string _loc, std::string _desc, std::string _snp, std::string _snp_loc, double _p )
        : gene(_gene), loc(_loc), desc(_desc), snp(_snp), snp_loc(_snp_loc), p(_p) { ; }

    std::string gene;
    std::string loc;
    std::string desc;
    std::string snp;
    std::string snp_loc;
    double p;

    std::string get_summary() {
        return gene + "\t" + loc + "\t" + desc + "\t" + snp + "\t" + snp_loc + "\t" + double2str(p) + "\n";
    }
};

#endif // SIGGENE_H

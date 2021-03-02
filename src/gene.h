#ifndef GENE_H
#define GENE_H

#include <string>

#include "interval.h"
#include "helper.h"

class Gene
{

public:
    Gene() {; }

    Gene(Interval gene, std::string desc)
        : gene(gene) , desc(desc) { }

    Interval gene;

    std::string desc;

    std::string get_chr_loc() { return "chr"  + int2str(gene.chr) + ":" + int2str(gene.bp1) + ".." + int2str(gene.bp2); }
    std::string get_loc() { return gene.get_desc(); }
    std::string get_desc() { return desc; }

};

#endif // GENE_H

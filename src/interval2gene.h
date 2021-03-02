#ifndef INTERVAL2GENE_H
#define INTERVAL2GENE_H

#include <string>

class Interval2Gene {
public:
    Interval2Gene()
        : interval("na"), n_snps(-1), n_altpos(-1), n_genes(-1), gene_id("."), gene_desc(".") {;}
    Interval2Gene(std::string _int, int _n_pos, int _n_altpos, int _n_genes)
        : interval(_int), n_snps(_n_pos), n_altpos(_n_altpos), n_genes(_n_genes), gene_id("."), gene_desc(".") {;}
    Interval2Gene(std::string _int, int _n_pos, int _n_altpos, int _n_genes, std::string _gene_id, std::string _gene_desc)
        : interval(_int), n_snps(_n_pos), n_altpos(_n_altpos), n_genes(_n_genes), gene_id(_gene_id), gene_desc(_gene_desc) {;}

    std::string interval;
    int n_snps;
    int n_altpos;
    int n_genes;
    std::string gene_id;
    std::string gene_desc;

    std::string get_inrich_summary() {
        return interval + "\t"
               + int2str(n_snps) + "\t"
               + int2str(n_altpos) + "\t"
               + int2str(n_genes) + "\t"
               + gene_id + "\t"
               + gene_desc + "\n";
    }
};

#endif // INTERVAL2GENE_H

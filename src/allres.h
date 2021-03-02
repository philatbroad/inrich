#ifndef ALLRES_H
#define ALLRES_H

#include "pathwayres.h"
#include "generes.h"
#include "intervalres.h"

class AllRes
{
public:
    AllRes();
    AllRes(std::string target_desc, std::string p, std::string gene_id, std::string gene_loc, std::string gene_desc, std::string interval_desc);
    ~AllRes();

    void set_data(std::string target_desc, std::string p, std::string gene_id, std::string gene_loc, std::string gene_desc, std::string interval_desc);
    //void set_data(std::string _id, std::string _all_loc, std::string _all_desc);


    bool operator==(const AllRes & rhs) const {
        if( pathway_id.compare(rhs.pathway_id)==0 &&
            gene_id.compare(rhs.gene_id)==0 &&
            interval_id.compare(rhs.interval_id)==0) {
            return true;
        }
        else {
            return false;
        }
    }

    // sort-order
    bool operator<(const AllRes & rhs) const
        {
            if( pathway_p < rhs.pathway_p ) return true;

            if( pathway_name.compare(rhs.pathway_name)==0 ) {
                if ( gene_id.compare(rhs.gene_id)<0 ) return true;
                else return false;
            }
            else if( pathway_name.compare(rhs.pathway_name)<0 ) {
                return true;
            }
            else {
                return false;
            }
        }


    std::string get_pathway_id() { return pathway_id; }
    std::string get_pathway_name() { return pathway_name; }
    double get_pathway_p() { return pathway_p; }
    std::string get_gene_id() { return gene_id; }
    std::string get_gene_name() { return gene_symbol; }
    std::string get_gene_desc() { return gene_desc; }
    std::string get_gene_loc() { return gene_chr + ":" + int2str(gene_start) + ".." + int2str(gene_end); }
    std::string get_interval_id() { return interval_id; }
    std::string get_interval_loc() { return interval_loc; }
    std::string get_interval_desc() { return interval_id + ":" + interval_loc; }

    PathwayRes* get_pathway_result();
    GeneRes* get_gene_result();
    IntervalRes* get_interval_result();

    std::string get_inrich_summary() {
        std::string sum =  get_interval_desc() + "\t"
                    + get_gene_loc() + "\t"
                    + gene_id + "\t"
                    + gene_symbol + " " + gene_desc + "\t"
                   + pathway_id + " " + pathway_name + "\t" + double2str(pathway_p) + "\n";
        return sum;
    }

    std::string get_aligator_summary() {
        std::string tmp;
        std::string sum = pathway_id + " " + pathway_name + "\t"
                     +  double2str(pathway_p) + "\t"
                    + gene_id + "\t"
                    + get_gene_loc() + "\t"
                    + gene_symbol + " " + gene_desc + "\t"
                    + get_interval_desc() + "\n";
        return sum;
    }


private:

    std::string pathway_id;
    std::string pathway_name;
    double pathway_p;

    std::string gene_id;
    std::string gene_symbol;
    std::string gene_desc;
    std::string gene_chr;
    int gene_start;
    int gene_end;

    std::string interval_id;
    std::string interval_loc;
};

#endif // RESULT_H

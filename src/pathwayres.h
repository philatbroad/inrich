#ifndef PATHWAYRES_H
#define PATHWAYRES_H

#include <string>
#include "helper.h"

class PathwayRes
{
public:
    PathwayRes() { }
    PathwayRes(std::string _id);
    PathwayRes(std::string _id, std::string _name, double _p);
    PathwayRes(std::string  _target_size, std::string _n_intervals, double _p, double _pcorr, std::string _low_warning, std::string  _unsh, std::string _target_id, std::string _target_name, double _exp_hit = -1, double _exp_gene_no = -1);
    PathwayRes(int _target_size, int _n_intervals, double _p, double _pcorr, double _low_warning, double _unsh, std::string _id_name, double _exp_hit = -1, double _exp_gene_no = -1) ;
    PathwayRes(std::string _target_size, std::string _n_intervals, double _p, double _pcorr, std::string _low_warning, std::string _unsh, std::string _id_name, double _exp_hit = -1, double _exp_gene_no = -1);

    void set_data(std::string _target_size, std::string _n_intervals, double _p, double _pcorr, std::string _low_warning, std::string  _unsh, std::string _target_id, std::string _target_name, double _exp_hit = -1, double _exp_gene_no = -1);

    double p;
    double pcorr;
    double exp_hit;
    double exp_gene_no;
    std::string t_target;
    std::string n_int;
    std::string low_warn;
    std::string unsh;
    std::string id;
    std::string name;

    bool operator==(const PathwayRes & rhs) const {
        if( id.compare(rhs.id)==0 ) {
            return true;
        }
        else {
            return false;
        }
    }

    // sort-order
    bool operator<(const PathwayRes & rhs ) const  {
            if ( p < rhs.p ) return true;

            if ( name.compare(rhs.name) < 0 ) return true;
            else return false;
    }

    std::string get_pathway_details();

    std::string get_aligator_summary() {
        return  t_target + "\t" +
         n_int + "\t" +
        double2str(exp_gene_no) + "\t" +
        double2str(p) + "\t" +
        double2str(pcorr) + "\t" +
        double2str(exp_hit) + "\t" +
        //low_warn + "\t" +
        id + " " + name + "\n";
    }

    std::string get_inrich_summary() {
        return t_target + "\t" +
                n_int + "\t" +
        double2str(p) + "\t" +
        double2str(pcorr) + "\t" +
        low_warn + "\t" +
        unsh + "\t" +
        id + " " + name + "\t" +
        double2str(exp_hit) + "\n";
    }

};

#endif // PATHWAYRES_H

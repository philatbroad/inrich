#ifndef ALIGATORMAIN_H
#define ALIGATORMAIN_H

#include <map>
#include <string>
#include <vector>
#include <cstdlib>

#include "loaders.h"
#include "crandom.h"
#include "gene.h"
#include "snp.h"
#include "siggene.h"
#include "allres.h"
#include "pathwayres.h"

class AligatorMain
{
public:
    explicit AligatorMain();
    ~AligatorMain();

    InputData globals;


    std::vector<AllRes> main_data;
    std::vector<PathwayRes> path_data;
    std::vector<SigGene> snp_data;

    void aligator_test();

    std::string get_sig_hit() { return sig_hit; }
    void set_globals(InputData & _globals);

    std::string summary_filename;

protected:
    std::set<SNP> sig_snps;
    std::set<SNP> all_snps;
    std::map<std::string, std::set<std::string> > targets;
    std::set<std::string> target_genes;
    std::map<int, std::set<Interval> > chr_genes;
    std::set<std::string> background_genes;
    std::map<std::string,Gene> all_genes;
    std::map<std::string, std::vector<Gene> > snp2genes;
    std::vector<std::string> snps_in_genes;
    std::map<std::string, std::set<SNP> > sig_genes;

    unsigned int n_targets;
    unsigned int n_sig_genes;

    std::vector<unsigned int> org_count;
    std::vector< std::map<int, int> > rep_dist;
    std::vector<double> emp_p;
    std::vector<double> corr_p;
    std::vector<double> exp_hit_count;
    std::vector<unsigned int> org_study_hit;
    std::vector<double> p_sig_study_hit;


    void init_data();

    bool load_files();
    bool load_file(int _modeNum, int _mode);
    int load_snps( std::string _snp_type, std::string _file, std::set<SNP> & snps );
    int load_targets( );
    int load_genes_with_chr_genes(  );

    bool gen_snp2gene();
    bool gen_sig_genes(  );
    bool update_fake_name( );
    bool cal_org_count();
    bool first_permutation();
    bool second_permutation();
    void make_cumul();


    int  count_overlap( std::map<std::string, std::set<SNP> > & genes, std::set<std::string> & pathway_genes );
    std::map<std::string, std::set<SNP> > gen_random_set( );
    std::vector<int>  gen_random_sets( int _observed );

    bool write_output( );
    bool log_target_path_data( );
    bool log_snp_data( );

    void log_summary_file();
    void log_main_file(std::ofstream & OUT);
    void log_gene_file(std::ofstream & OUT);
    void log_target_file(std::ofstream & OUT);
    void log_header();
    void report_sig_study_hit();

    void log(std::string _msg);
    std::string sig_hit;

    bool fake_name;

};

#endif // ALIGATORMAIN_H

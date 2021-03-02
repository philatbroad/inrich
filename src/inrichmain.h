#ifndef INRICHMAIN_H
#define INRICHMAIN_H


#include <string>
#include <vector>
#include <map>
#include <ostream>

#include "inputdata.h"
#include "common.h"
#include "gene.h"
#include "snp.h"
#include "interval.h"
#include "pairintervaldist.h"

#include "allres.h"
#include "pathwayres.h"

#include "interval2gene.h"
#include "log.h"

class InRichMain
{
public:
    explicit InRichMain();
    InRichMain(InputData & _globals);
    ~InRichMain();


    InputData globals;

    std::vector<AllRes> main_data;
    std::vector<PathwayRes> path_data;
    std::vector<Interval2Gene> interval_data;

    int load_background_data();
    int load_range_data();
    int load_assoc_data();
    int load_assoc_snps_data();
    int filter_non_map_data();
    bool bi_search( int min_ix, int max_ix, int bp1, int bp2 );

    void inrich_test();
    void clustering_test();
    //bool snp2interval();
    //void update_chr_bp_orig_snp( const std::string dbname );

    bool is_data_ok() {
        if ( globals.convert_mode ) {
            if( orig_snp.size() == 0 || snps.size() == 0 ) {
              return false; //"No Input Data: interval/genes/targets";
            }
        }
        else if ( globals.multi_set ) {
            if( all_genes.size() == 0 || targets.size() == 0 ) {
              return false; //"No Input Data: interval/genes/targets";
            }
        }
        else if ( globals.exam_bp_dist ) {
            if( all_genes.size() == 0 || orig.size() == 0 ) {
              return false; //"No Input Data: interval/genes";
            }
        }
        else {
            if( orig.size() == 0 || all_genes.size() == 0 || targets.size() == 0 ) {
              return false; //"No Input Data: interval/genes/targets";
            }

        }
        return true;
    }

    std::string my_log;

    std::string get_log() { return my_log; }
    std::string get_sig_hit() {return sig_hit;}
    std::string get_summary_file() { return summary_filename; }

protected:
    static const int filename_size = 20;
    std::string err_msg;

    std::set<Interval> test;
    std::set<Interval>  get_all_targets();
    std::vector<Interval> snps;

    std::set<std::string> background_genes;
    std::map<std::string,Gene> all_genes;
    std::set<Interval> all_genes_set;
    std::map<std::string, std::set<Interval> > targets;
    std::map<std::string, std::string> gene_results;
    std::set<Interval> orig;
    std::set<SNP> orig_snp;
    std::map<Interval, int> match_tracker;
    std::set<Interval> all_target_genes;

    std::set<PairIntervalDist> test_dist_interval;
    int val_dist_p_no;
    std::vector<int> test_dist_p;

    int n_interval;
    int n_set;


    std::vector< std::map<int,int> >  null_tracker;
    std::vector<int> count;
    std::vector<int> actual_count;

    void init_first_permutation();
    void init_second_permutation();
    std::vector<int> pcount;
    std::vector<int> p_stat_sum;
    std::vector<int> pcorr;


    std::vector<unsigned int> org_study_hit;
    std::vector<unsigned int> p_sig_study_hit;
    std::vector<unsigned int> org_study_gene_hit;
    std::vector<unsigned int> p_sig_study_gene_hit;

    void cal_study_sig_hit();
    void get_overlapping_genes( const std::set<Interval> & test ,
                                const std::set<Interval> & target ,
                                std::set<Interval> & otarget );

    std::map< std::string , Interval> of_interest;
    std::map< std::string , double> of_interest_bestp;
    std::map< std::string , std::string> of_interest_bestset;


    std::string main_filename;
    std::string interval_filename;
    std::string target_filename;
    std::string gene_filename;
    std::string log_filename;
    std::string summary_filename;
    void log_main_file();
    void log_interval_file();
    void log_target_file();
    void log_gene_file();

    void log_main_file(std::ofstream & OUT);
    void log_interval_file(std::ofstream & OUT);
    void log_target_file(std::ofstream & OUT);
    void log_gene_file(std::ofstream & OUT);

    void log_summary_file(bool _append=true);
    void log_summary_file(std::string _filename, bool _append=true);

    //bool make_lddb( const std::string & name, const std::string & dbname );
    //std::string ld_clumping( const std::string & lddb_name, std::set<SNP> & orig_snp);
    //std::string ld_clumping_2( const std::string & lddb_name, std::set<SNP> & orig_snp );
                                          

    std::string sig_hit;

    std::string register_interval ( std::vector<Interval> & int_list );

    // ------------------------------
    // moved to MainInputDialog
    // 2011-Feb
    // ------------------------------
    //QProgressDialog *progress_dialog;
    //bool cancel_progress;
    //void init_progress_dialog();
    //void remove_progress();


    // ===========================================
    // File Loading Functions : interval-func.cpp
    // ===========================================
    std::set<Interval> get_set( std::map< std::string, Gene > & );

    int num_in_interval( const Interval & interval ,
                         const std::set<Interval> & regions );

    std::map<Interval,std::set<Interval> >
    get_overlap( const std::set<Interval> & test ,
                 const std::set<Interval> & target ) ;

    std::set<Interval>
    match( const std::set<Interval> & t ,
           const std::vector<Interval> & pos ,
           std::map<Interval,int> & match_tracker ,
           std::set<Interval> * genes = NULL );

    void
    match_interval( std::set<Interval> & m ,
           std::set<Interval>::iterator & i,
           const std::vector<Interval> & pos ,
           const int s,
           std::map<Interval,int> & match_tracker ,
           std::set<Interval> * genes = NULL );

    bool overlap( const std::set<Interval> & s );

    bool
    overlap_interval( const std::set<Interval> & s,
                      std::set<Interval>::iterator & jtr );


    std::set<Interval> merge( const std::set<Interval> & intervals );

    int evaluate( const std::set<Interval> & test ,
                  const std::set<Interval> & target );

    // =========================================
    // PL 2010.07.18 ADDED
    // =========================================
    std::set<Interval> remove_subinterval( const std::set<Interval> & regions);


    // =========================================
    // PL 2010.07.29 ADDED
    // =========================================
    std::map< std::string, std::set<Interval> > remove_nonrelevant_targets(
         std::set<Interval> & _test,
         std::map< std::string, std::set<Interval> > & _org_targets,
         std::map< Interval, std::set<Interval> > & _overlap);

    // =========================================
    // PL 2010.08.18 ADDED
    // =========================================
    std::vector< std::vector<int> > init_study_hit( int nrep );
    void count_sig_study_hit_all( const int rep,
                              const int nrep,
                              std::vector< std::vector<int> > & sig_study_hit,
                              const std::vector<int> & pcount );
    void count_sig_study_hit( const int rep,
                              const int nrep,
                              std::vector< std::vector<int> > & sig_study_hit,
                              const int pcount );

    // =========================================
    // PL 2010.09.07 UPDATED
    // =========================================
    void report_sig_study_hit();

    std::string header();

    // ===========================================
    // PairIntervalDist : interval-func.cpp
    // ===========================================
    std::set<PairIntervalDist> compute_dist_interval ( std::set<Interval> & regions );
    std::set<PairIntervalDist> compute_dist_interval_adjacent ( std::set<Interval> & regions );
    std::set<PairIntervalDist> compute_dist_interval_all ( std::set<Interval> & regions );
    std::vector<int> init_test_p ( int nsize );
    int update_dist_p( int no_dist_p, std::vector<int> & dist_p, std::set<PairIntervalDist> & dist1, std::set<PairIntervalDist> &dist2 );
    void report_dist_p( int val_p_no,
                       std::vector<int> & dist_p,
                       std::set<PairIntervalDist> & test_dist,
                       std::set<Interval> & test,
                       std::string filename );



    std::set<Interval> make_compact( int cnt, std::set<Interval> & orig, std::set<Interval> & genes);

    bool log_intervals(std::string _file);
    // generate smaller region -- interval within tag SNP
    int get_closest_ref_snp_pos( int _snp_no, int _chr, int _pos, bool _left, int & start_pos );

    // generate larger region -- interval beyond tag SNP
    int get_closest_ref_snp_pos_2( int _snp_no, int _chr, int _pos, bool _left, int & start_pos );

    void init_data();

    bool load_files();
    bool load_file(int _mode);
    bool merge_data();
    void update_all_genes_map();
    void gen_test_set(std::map< Interval , std::set<Interval> > & overlap1, std::map< Interval , std::set<Interval> > & overlap2);
    std::set<Interval> gen_match_set(std::map< Interval , std::set<Interval> > & overlap1,
                                   std::map< Interval , std::set<Interval> > & overlap2,
				   std::set<Interval> & all_target_genes,
				                                  std::set<Interval> & match );

    bool init_optional_load();


    bool assign_snp_count(  std::set<Interval> & t ,
                                            const std::vector<Interval> & pos );
    bool pre_compute( const std::set<Interval> & t ,
                                          const std::vector<Interval> & pos ,
                                          std::set<Interval> * genes );
    void log(std::string _str1);
    bool run_clustering_permutation();
    bool run_first_permutation ();
    bool run_second_permutation ();
    bool write_output();
    bool report_clustering();
    void cal_unplaced();
    bool log_interval_data();
    bool log_target_path_data();

    std::string get_gene_symbol( Interval & _id );

};

#endif // INRICHMAIN_H

#ifndef INRICHMULTIMAIN_H
#define INRICHMULTIMAIN_H

#include <map>
#include <set>
#include <vector>
#include <string>

#include "interval.h"
#include "snp.h"
#include "inputdata.h"
#include "inrichmain.h"

class InRichMultiMain : public InRichMain
{
public:
    explicit InRichMultiMain();
    InRichMultiMain(InputData & _globals);
    ~InRichMultiMain();

    std::map<std::string, std::set<Interval> > multi_orig;
    std::map<std::string, std::set<SNP> > multi_orig_snps;
    std::map<std::string, std::set<Interval> > multi_test;

    std::map<std::string, std::map<unsigned int,std::vector<unsigned int> >  > multi_null_tracker;

    std::map<std::string, std::vector<int> > multi_count;

    std::map<std::string, std::vector<Interval> > multi_snps;
    std::map<std::string, std::map<int,Interval> > multi_scaffold;

    void inrich_multi_test();


protected:
    int n_multi_set_no;
    std::vector<std::string> path_summary;

    bool load_file(int _mode);
    int load_multi_assoc();
    int load_multi_assoc_snps();
    int load_multi_snps( std::string _list_file );
    int load_snp_map( std::string _set, std::string _mapfile );

    bool merge_data();

    bool init_optional_load();

    void init_first_permutation();
    int evaluate( const std::set<Interval> & test ,
                              const std::set<Interval> & target );
    bool run_first_permutation ();

    void init_second_permutation();
    bool run_second_permutation();

    bool write_output();
    bool log_path_data();
    std::string header(bool _line=true);
    void log_summary_file(std::string _filename, bool _append=true);
    void log_target_file(std::ofstream & OUT);
    void log_main_file(std::ofstream & OUT);
    void log_main_data(std::ofstream & OUT , std::string _set);
    void log_interval_file(std::ofstream & OUT);
    void log_interval_data(std::ofstream & OUT, std::string _set) ;

    bool load_files();

    bool is_data_ok();

    int get_random_match_no( std::string _set_name, int _pway );
    std::map<std::string, int> get_random_match_no_list( int _pway );

    void cal_study_sig_hit();
    //void report_sig_study_hit();

};

#endif // INRICHMULTIMAIN_H

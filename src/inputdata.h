#ifndef INPUTDATA_H
#define INPUTDATA_H

#include <map>
#include <set>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>


#include "common.h"
#include "interval.h"

class InputData
{
public:
     InputData();
     std::string get_cmd();
     std::string get_qt_cmd();

     bool verbose;
     bool precompute;
     bool is_snps;
     int topn;
     unsigned int tsize_min;
     unsigned int tsize_max;
     test_t test;
     bool use_map;
     bool background_list;
     int largest_gene;
     bool match_genes;
     bool use_range;
     int window;
     int nrep;
     int seed;
     double tol;
     double pthresh;
     unsigned int min_cnt;
     bool multi_set;
     bool multi_mapfile;
     //int  multi_testmode;

     bool is_human;

     long unsigned int total_seq;
     std::map<int,Interval> scaffold;
     std::vector<int> chr_list;
     std::map<Interval,std::vector<int> > acceptable_positions;

     //bool included( const Interval & i ) const;
     Interval random() const;


     bool exam_bp_dist;
     int exam_bp_dist_top;

     std::string mapfile;
     std::string testfile;
     std::string genefile;
     std::string targetfile;
     std::string rangefile;
     std::string backgroundfile;
     std::string create_db_name;
     std::string use_db_name;
     std::string outroot;

     bool compact;
     bool aligator;
     int  nboot;

     bool printPermutations;

     //bool create_lddb;
     //bool use_lddb;
     bool convert_mode;
     int  input_tab_index;

     void set_default_param(std::string _project_title,
                                     bool _aligator,
                                     std::string _testFile,
                                     std::string _geneFile,
                                     std::string _mapFile,
                                     std::string _targetFile );

     void set_advanced_param(std::string _project_title,
                                      bool _aligator,
                                      std::string _testFile,
                                      std::string _geneFile,
                                      std::string _mapFile,
                                      std::string _targetFile,
                                      std::string _rangeFile,
                                      std::string _backFile,
                                      int _minPathSize,
                                      int _maxPathSize,
                                      int _minCnt,
                                      int _window,
                                      double _p,
                                      int _nRep,
                                      int _nBoot,
                                      int _seed,
                                      bool _intervalMode,
                                      bool _preCompute,
                                      bool _matchGene,
                                      double _snpMappingDensity,
                                      int _topN,
                                      bool _physicalClustering,
                                      int _clusteringTop);

     void set_convert_param(
                    std::string _outroot,
                     std::string _testFile,
                     std::string _mapFile,
                     bool _create_db,
                     std::string _create_dbName,
                     std::string _useDBName ) ;

     void set_aligator_param(std::string _project_title,
                                      std::string _testFile,
                                      std::string _mapfile,
                                      std::string _geneFile,
                                      std::string _targetFile,
                                      int _minPathSize,
                                      int _maxPathSize,
                                      int _minCnt,
                                      int _window,
                                      int _nRep,
                                      int _nboot,
                                      int _seed);



     std::string check_input_data();
     std::string  check_file(std::string _file, std::string _filetype, bool _should = true );

     void set_data( int argc, char* argv[] );
     int parse_arg( int argc, char* argv[] );
     int check_arg();
     void set_arg();
     void log_arg();
     bool is_multi_mapfile( );

protected:
     std::string get_aligator_cmd();
     std::string get_inrich_cmd();
     std::string get_aligator_qt_cmd();
     std::string get_inrich_qt_cmd();

 };

#endif // INPUTDATA_H

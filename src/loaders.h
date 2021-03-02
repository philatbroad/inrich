#ifndef __INRICH_LOADER_H__
#define __INRICH_LOADER_H__

#include "interval.h"
#include "gene.h"
#include "snp.h"
#include "inputdata.h"

#include <map>
#include <vector>
#include <set>
#include <string>

enum file_type {
    BACKGROUND = 0 ,
    RANGE = 1 ,
    GENE = 2,
    SNP_MAP = 3 ,
    TARGET_SET = 4,
    ASSOC = 5,
    ASSOC_SNP = 6,
    ALI_GENE = 7,
    ALI_SNP_MAP = 8,
    ALI_TARGET_SET = 9,
    ALI_ASSOC = 10
} ;

int get_hum_chr(std::string _chrcode);
std::string get_hum_chrcode(int _chr);
int get_chr(std::string _chrcode, bool _human);

bool is_valid_chr(int _chr, bool _human);
bool is_valid_hum_chr(int _chr);

int 
load_ranges( const std::string &, InputData & globals, std::string & err_msg  );

std::vector<Interval> 
load_snps( const std::string & mapfile, InputData & globals , std::string & err_msg );

std::map<std::string,Gene> 
loader_genes( const std::string & genefile ,
            std::set<std::string> & bg, InputData & globals  , std::string & err_msg );

std::map<std::string, std::set<Interval> > 
loader_targets( const std::string & targetfile ,
              std::map<std::string,Gene> &, InputData & globals , std::string & err_msg, std::string & res_msg   );

std::set<Interval> 
load_assoc( const std::string & testfile, InputData & globals  , std::string & err_msg  );

std::set<std::string> 
load_background_genes( const std::string & , std::string & err_msg  );

std::set<SNP> 
load_assoc_snps( const std::string & testfile, InputData & globals , std::string & err_msg  );

std::map<std::string, SNP>
load_ld_table( const std::string & ldfile, InputData & globals );

#endif

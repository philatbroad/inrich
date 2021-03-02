#ifndef __INRICH_HELPER_H__
#define __INRICH_HELPER_H__

#include <string>
#include <vector>
#include <sstream>
#include <set>
#include <map>

#include "interval.h"

std::string get_chr(int _chr); 
bool is_all_numeric(std::string _str);
int read_line_num(std::string _file_name);
std::vector<std::string> tokenize_string(std::string _str, char _sep);
std::string extract_symbol(std::string _desc);
std::string extract_desc(std::string _desc);
std::string cur_date();

std::string int2str(int n);
std::string double2str(double n);
bool str2int( const  std::string & , int & );
bool str2double( const  std::string & , double & );
std::vector<std::string> tokenizeLine(std::ifstream & F1);

template <class T>
bool from_string(T& t,
		 const std::string& s,
		 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

bool fileExists( std::string f);

std::string first2upper_str(std::string s);
void stoupper(std::string & s);

int get_step(int num);

//std::set<Interval> get_set( std::map< std::string,Gene> & g );

#endif



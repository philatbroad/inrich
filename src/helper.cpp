
#include "helper.h"

#include <sstream>
#include <istream>
#include <iostream>
#include <string>
#include <fstream>
#include <ctime>


std::string get_chr(int _chr) {
	if(_chr>=1 && _chr<=26) {
		return int2str(_chr);
	}
	else if(_chr==123) {
		return "X";
	}
	else if(_chr==124) {
		return "Y";
	}
	else if(_chr==125) {
		return "XY";
	}
	else if(_chr==126) {
		return "M";
	}

	return "NA";
}


bool is_all_numeric(std::string _str) {
     for(unsigned int i = 0; i < _str.length(); i++) {
       		if(!std::isdigit(_str[i])) {
			return false;
		}
      }

      return true;
}


int get_step(int _num)
{
    int step=_num/100;

    return step <= 0 ? 1 : step;
}


int read_line_num(std::string _file_name) {
    std::ifstream OUT( _file_name.c_str(), std::ios::in );
    int num=0;
    std::string sline;
    while( std::getline(OUT,sline)) {
        num++;
    }
    //TARGET.seekg(0, std::ios::beg);  // This doesn work, here. So, head to close and reopen
    OUT.close();
    return num;
}

std::vector<std::string> tokenize_string(std::string _str, char _sep) {
    std::istringstream iss(_str);
    std::string token;
    std::vector<std::string> tokens;
    while ( getline(iss, token, _sep) )
    {
        tokens.push_back(token);
    }
    return tokens;
}

std::string extract_symbol(std::string _desc) {
    std::istringstream iss(_desc);
    std::string symbol;
    getline(iss, symbol, ' ');
    return symbol;
}

std::string extract_desc(std::string _desc) {
    std::istringstream iss(_desc);
    std::string token;
    std::string desc;
    getline(iss, token, ' ');
    while ( getline(iss, token,  ' ') )
    {
        desc = desc + token + " ";
    }
    return desc;
}


std::string cur_date() {
    time_t rawtime;

    time ( &rawtime );
    char buf[50];
    sprintf ( buf, "%s", ctime (&rawtime) );
    std::string cur_date = buf;
    return cur_date;

    /*
    struct tm *ptr;
       time_t sec;
       time(&sec);
       ptr = localtime(&sec);

       int month = (int) ptr->tm_mon + 1;
       std::string month_str;
       switch(month) {
       case 1:
           month_str = "Jan";
           break;
       case 2:
           month_str = "Feb";
           break;
       case 3:
           month_str = "Mar";
           break;
       case 4:
           month_str = "Apr";
           break;
       case 5:
           month_str = "May";
           break;
       case 6:
           month_str = "June";
           break;
       case 7:
           month_str = "July";
           break;
       case 8:
           month_str = "Aug";
           break;
       case 9:
           month_str = "Sep";
           break;
       case 10:
           month_str = "Oct";
           break;
       case 11:
           month_str = "Nov";
           break;
       case 12:
           month_str = "Dec";
           break;

       }

       int day   = (int) ptr->tm_mday;
       int year  = (int) ptr->tm_year + 1900;

       return int2str(day) + "-" + month_str + "-" + int2str(year);*/
}

std::string int2str(int n)
{
    std::ostringstream s2( std::stringstream::out );
    s2 << n;
    return s2.str();
}

std::string double2str(double n)
{
    std::ostringstream s2( std::stringstream::out );
    s2 << n;
    return s2.str();
}

bool str2int( const std::string & s , int & i )
{
    return from_string<int>(i,s,std::dec);
}

bool str2double( const std::string & s , double & i )
{
    return from_string<double>(i,s,std::dec);
}

std::vector<std::string> tokenizeLine(std::ifstream & F1)
{
    std::string sline;
    std::getline(F1,sline);
    std::string buf; 
    std::stringstream ss(sline); 
    std::vector<std::string> tokens; 
    while (ss >> buf)
	tokens.push_back(buf);
    return tokens;
}

bool fileExists( std::string f)
{
    std::ifstream inp;    
    inp.open(f.c_str(), std::ifstream::in);
    if( inp.fail() )
    {
	inp.clear(std::ios::failbit);
	inp.close();
	return false;
    }
    inp.close();
    return true;
}

std::string first2upper_str(std::string s) {
   std::string first(1, (char)std::toupper(s[0]) );
   return first + s.substr(1);
}


void stoupper(std::string& s)
{
    std::string::iterator i = s.begin();

    while (i != s.end()) {
        *i = std::toupper((unsigned char)*i);
        ++i;
    }
}

/*
std::set<Interval> get_set( std::map< std::string,Gene> & g )
{
    std::set<Interval> s;

    std::map< std::string,Gene>::iterator i = g.begin();
    while ( i != g.end() )
    {
        s.insert( i->second.gene );
        ++i;
    }
    return s;
}
*/

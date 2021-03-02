#ifndef LOG_H
#define LOG_H

#include <fstream>
#include <iostream>

class Log : public std::ofstream {
 public:
    Log() {;}

    Log( const std::string & f  )
        : std::ofstream(f.c_str(), std::ios::app) { file_name = f;}

    Log( const std::string & f, bool _open  )
        : std::ofstream(f.c_str(), std::ios::out) { file_name = f;}
    ~Log() { close(); }

    template<class T>
    Log & operator<<(const T & msg)
        {
            // output to file; mirror to stdout
            // std::ofstream & g = *this;
            std::ofstream g(file_name.c_str(), std::ios::app);
            g << msg ;
            std::cout << msg;
            return *this;
      }

    std::string file_name;
};



#endif // LOG_H

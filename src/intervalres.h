#ifndef INTERVALRES_H
#define INTERVALRES_H

#include <string>

class IntervalRes
{
public:
    IntervalRes();
    IntervalRes(std::string _desc);
    IntervalRes(std::string _id, std::string _chr, int _start, int _end);
    IntervalRes(std::string _id, std::string _desc);

    //void set_data(std::string _id, int _chr, int _start, int _pos);
    void set_data(std::string _desc);

    std::string get_id() { return id; }
    std::string get_loc();

    std::string get_interval_details() { return "Interval ID : " + id + "\nLocus : " + get_loc(); }

    std::string id;
    std::string chr;
    int start;
    int end;

    bool operator==(const IntervalRes & rhs) const {
        if( id.compare(rhs.id)==0 ) {
            return true;
        }
        else {
            return false;
        }
    }

};




#endif //

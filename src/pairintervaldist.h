#ifndef PAIRINTERVALDIST_H
#define PAIRINTERVALDIST_H

#include <string>
#include "interval.h"
#include "common.h"

static int get_chr_end_pos ( int chr ) {

    // based on
    // http://www.ornl.gov/sci/techresources/Human_Genome/posters/chromosome/faqs.shtml
    int chr_end_pos[] = {-1,
                    245203898, // chr 1
                    243315028,
                    199411731,
                    191610523,
                    180967295,  // chr 5
                    170740541,
                    158431299,
                    145908738,
                    134505819,
                    135480874, // chr 10
                    134978784,
                    133464434,
                    114151656,
                    105311216,
                    100114055, // chr 15
                    89995999,
                    81691216,
                    77753510,
                    63790860,
                    63644868, // chr 20
                    46976537,
                    49476972,
                    152634166,   // chr X
                    50961097     // chr Y
    };

    return chr_end_pos[chr];
}

class PairIntervalDist
{
public:
    PairIntervalDist() {; }


    PairIntervalDist( std::string _name1, std::string _name2, long _dist )
        : name1(_name1), name2(_name2), dist(_dist) { }

    PairIntervalDist( const Interval & reg1 , const Interval & reg2 ) {
        name1 = reg1.name;
        name2 = reg2.name;

        if(reg1.chr == reg2.chr) {
            if(reg1.bp2 < reg2.bp1 ) {
                dist = reg2.bp1 - reg1.bp2;
            }
            else if( reg2.bp2 < reg1.bp1 ) {
                dist = reg1.bp1 - reg2.bp2;
            }
        }
        else if(reg1.chr < reg2.chr) {
            dist = ( get_chr_end_pos( reg1.chr ) - reg1.bp2 ) + reg2.bp1 + MAX_PAIR_DIST;
        }
        else {
            dist = ( get_chr_end_pos( reg2.chr ) - reg2.bp2 ) + reg1.bp1 + MAX_PAIR_DIST;
        }
    }

    std::string name1;
    std::string name2;
    long dist;


    // sort-order
    bool operator<(const PairIntervalDist & rhs) const {
        if ( dist <= rhs.dist ) return true;
        else                    return false;
    }

    // overlap
    bool operator==(const PairIntervalDist & rhs) const
    {
        if ( (name1.compare(rhs.name1)==0 && name2.compare(rhs.name2)) || (name1.compare(rhs.name2) && name2.compare(rhs.name1))  ) return true;
        else return false;
    }

    // print
    friend std::ostream & operator<<( std::ostream & out , const PairIntervalDist & i )
    {
        out << i.name1 << "\t" << i.name2 << "\t" << i.dist ;
        return out;
    }

};

#endif // PAIRINTERVALDIST_H

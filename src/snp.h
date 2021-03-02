#ifndef SNP_H
#define SNP_H

#include <string>
#include <fstream>

#include "helper.h"

// =====================================================
// LD-clumping based on HapMap
//
// PL 2010-10-18
// =====================================================

class SNP
{
public:
    double p;
    std::string name;
    int bp;
    int chr;
    int left_tag;
    int right_tag;

    SNP() {; }
    SNP( double p, std::string name )
        : p(p), name(name), bp(-1), chr(-1), left_tag(-1), right_tag(-1) { }
    SNP( int chr, int bp, std::string name  )
        : p(-1), name(name), bp(bp), chr(chr) , left_tag(-1), right_tag(-1) { }
    SNP( int chr, int bp, std::string name, double p )
        : p(p), name(name), bp(bp), chr(chr) , left_tag(-1), right_tag(-1) { }
    SNP( int chr, int bp, std::string name, int left_tag, int right_tag )
        : p(-1), name(name), bp(bp), chr(chr), left_tag(left_tag), right_tag(right_tag) { }
    SNP( int chr, int bp, std::string name, double p , int left_tag, int right_tag )
        : p(p), name(name), bp(bp), chr(chr), left_tag(left_tag), right_tag(right_tag) { }

    SNP( const SNP & rhs )
    {
        p = rhs.p;
        bp = rhs.bp;
        name = rhs.name;
        chr = rhs.chr;
        left_tag = rhs.left_tag;
        right_tag = rhs.right_tag;
    }



    // sort-order
    bool operator<(const SNP & rhs) const
    {
        if ( chr < rhs.chr ) return true;
        if ( chr==rhs.chr && bp < rhs.bp ) return true;

        return false;
    }
    // overlap
    bool operator==(const SNP & rhs) const
    {
        if ( name.compare(rhs.name)==0 ) return true;

        if (chr==rhs.chr && bp==rhs.bp) return true;
        else return false;
    }

    // print
    friend std::ostream & operator<<( std::ostream & out , const SNP & i )
    {
        if ( i.name != "" ) out << i.name << ":";
        out << "chr" << i.chr << ":" << i.bp << " p=" << i.p;
        return out;
    }

    void set_tag( int left, int right ) {
        left_tag = left;
        right_tag = right;
    }


    bool no_tag() {
        if(left_tag == -1 && right_tag==-1) return true;
        else return false;
    }

    std::string get_desc() { return name + ":chr" + int2str(chr) + ":" + int2str(bp); }
    std::string get_loc() { return "chr" + int2str(chr) + ":" + int2str(bp); }

};

#endif // SNP_H

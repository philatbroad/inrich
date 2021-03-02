#ifndef INTERVAL_H
#define INTERVAL_H

#include <string>
#include <fstream>


class Interval
{
public:
    Interval() {}
    Interval( int chr , int bp1 )
            : chr(chr), bp1(bp1) , bp2(bp1), n(-1), name("n/a") { }
    Interval( int chr , int bp1 , int bp2 )
        : chr(chr), bp1(bp1) , bp2(bp2), n(-1), name("n/a") { }
    Interval( int chr , int bp1 , int bp2 , std::string name )
        : chr(chr), bp1(bp1) , bp2(bp2), n(-1), name(name) { }
    Interval( int chr , int bp1 , int bp2 , int n )
        : chr(chr), bp1(bp1) , bp2(bp2) , n(n), name("n/a") { }
    Interval( int chr , int bp1 , int bp2 , int n , std::string name )
        : chr(chr), bp1(bp1) , bp2(bp2) , n(n) , name(name) { }

    Interval( const Interval & rhs )
        {
            chr = rhs.chr;
            bp1 = rhs.bp1;
            bp2 = rhs.bp2;
            n = rhs.n;
            name = rhs.name;
        }

    std::string get_desc();

    int chr;
    int bp1;
    int bp2;
    int n;
    std::string name;

    void combine( const Interval & rhs )
        {
            // do not merge with self, if name given
            if ( name != "" && name == rhs.name ) return;

            if ( rhs.chr != chr ) return;
            if ( rhs.bp1 < bp1 ) bp1 = rhs.bp1;
            if ( rhs.bp2 > bp2 ) bp2 = rhs.bp2;

            if ( rhs.name != "" )
                       {
                           if ( name != "" ) name += ",";
                           name += rhs.name;
                       }
      }

       // sort-order
       bool operator<(const Interval & rhs) const
           {
               if ( chr < rhs.chr ) return true;
               if ( chr > rhs.chr ) return false;
               if ( bp1 < rhs.bp1 ) return true;
               if ( bp1 > rhs.bp1 ) return false;
               return bp2 < rhs.bp2;
           }

       // overlap
       bool operator==(const Interval & rhs) const
           {
               if ( chr != rhs.chr ) return false;
               return bp1 <= rhs.bp2 && bp2 >= rhs.bp1;
           }

       // print
       friend std::ostream & operator<<( std::ostream & out , const Interval & i )
           {
               if ( i.name != "" ) out << i.name << ":";
               out << "chr" << i.chr << ":" << i.bp1 << ".." << i.bp2;
               return out;
           }


};

#endif // INTERVAL_H

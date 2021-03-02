#include "lddb.h"
#include "sqlite3.h"

#include <iostream>
#include <fstream>


/*
int main( int argc, char * argv[] )
{
  std::string db_dir = argv[1];

  std::cout << "Loading LDDB ... \n";

  LDDB lddb( db_dir + "ld1.db", "/psych/genetics_data/projects/pathways/src/" );

  // commented out line that created the database
  //if( argv[0] == "create.db" ) {
  	std::cout << "Creating DB ... \n";

  	lddb.load( db_dir + "ld.table" );
  	exit(0);
  //}

  uint64_t bp1 , bp2;
  bool found;

  // true SNP = rs2220306 10 3730909 3716611 3756696

  // expect to be okay 

  found = lddb.query( 10 , 3730909 , bp1 , bp2 );  
  if ( found ) std::cout << bp1 << "\t" << bp2 << "\n";
  else std::cout << "problem\n";

  found = lddb.query( 10 , 3730909 , "rs2220306" , bp1 , bp2 );  
  if ( found ) std::cout << bp1 << "\t" << bp2 << "\n";
  else std::cout << "problem\n";


  // expect problems

  found = lddb.query( 10 , 373091 , bp1 , bp2 );  
  if ( found ) std::cout << bp1 << "\t" << bp2 << "\n";
  else std::cout << "problem\n";

  found = lddb.query( 10 , 3730909 , "rs1220306" , bp1 , bp2 );  
  if ( found ) std::cout << bp1 << "\t" << bp2 << "\n";
  else std::cout << "problem\n";

} */


bool SQL::open( const std::string & n , const std::string & scratch )
{
    rc = sqlite3_open(n.c_str(), &db);
    if ( rc )  std::cerr << db_err << "\n";

    if ( scratch != "" )
    {
        query( "PRAGMA temp_store_directory = '"
               + scratch + "';");
    }
    return rc == 0;
}

void SQL::synchronous(bool b)
{
    if ( !b )
        query( "PRAGMA synchronous=OFF;" );
    else
        query( "PRAGMA synchronous=FULL;" );
}

void SQL::close()
{
    if ( db )
    {
        sqlite3_close(db);
        db = NULL;
    }
}

bool SQL::query(const std::string & q)
{  
    char * db_err;
    rc = sqlite3_exec( db , q.c_str() , 0 , 0 , &db_err );
    if ( rc ) std::cerr << "SQLITE error: " << db_err << "\n";
    return rc == 0;
}

sqlite3_stmt * SQL::prepare(const std::string & q)
{   
    sqlite3_stmt * p;
    int rc = sqlite3_prepare_v2( db , q.c_str() , q.size() , &p , NULL );
    if ( rc ) std::cerr << "SQLITE-PREPARE ERROR: " << sqlite3_errmsg(db) << "\n";
    return rc ? NULL : p;
}

void SQL::begin()
{  
    char * db_err;
    std::string q = "BEGIN;";
    rc = sqlite3_exec( db , q.c_str() , 0 , 0 , &db_err );
    if ( rc ) std::cout << db_err << "\n";
}

void SQL::finalise(sqlite3_stmt * stmt)
{
    sqlite3_finalize(stmt);
}

bool SQL::step(sqlite3_stmt * stmt)
{

    rc = sqlite3_step( stmt );

    if ( rc == SQLITE_ERROR )
    {
        reset(stmt);
        std::cerr << "SQL error: " << sqlite3_errcode(db) << "\n";
    }

    if ( rc != SQLITE_ROW && rc != SQLITE_DONE )
        std::cerr << "STEP: " << rc << " " << sqlite3_errmsg(db) << "\n";

    return rc == SQLITE_ROW;
}

void SQL::reset( sqlite3_stmt * stmt )
{
    sqlite3_reset( stmt );
}


void SQL::bind_int( sqlite3_stmt * stmt , const std::string & index , const int value )
{
    sqlite3_bind_int( stmt ,
                      sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
                      value );
}

void SQL::bind_null( sqlite3_stmt * stmt , const std::string & index  )
{
    sqlite3_bind_null( stmt ,
                       sqlite3_bind_parameter_index( stmt , index.c_str() ) );
}

void SQL::bind_int64( sqlite3_stmt * stmt , const std::string & index , const uint64_t value )
{
    sqlite3_bind_int64( stmt ,
                        sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
                        value );
}

void SQL::bind_double( sqlite3_stmt * stmt , const std::string & index , const double value )
{
    sqlite3_bind_double( stmt ,
                         sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
                         value );
}

void SQL::bind_text( sqlite3_stmt * stmt , const std::string & index , const std::string & value )
{
    sqlite3_bind_text( stmt ,
                       sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
                       value.c_str() ,
                       value.size() ,
                       0 );
}


int SQL::get_int( sqlite3_stmt * stmt , int idx )
{
    return sqlite3_column_int( stmt , idx );
}

uint64_t SQL::get_int64( sqlite3_stmt * stmt , int idx )
{
    return sqlite3_column_int64( stmt , idx );
}

double SQL::get_double( sqlite3_stmt * stmt , int idx )
{
    return sqlite3_column_double( stmt , idx );
}

std::string SQL::get_text(  sqlite3_stmt * stmt , int idx )
{
    const unsigned char * s = sqlite3_column_text( stmt , idx );
    if ( s == NULL )
        return "";
    else
        return (const char*)s;
}

void SQL::commit()
{
    query( "COMMIT;" );
}

LDDB::LDDB () {

}

LDDB::LDDB( const std::string & f, const std::string & s )
{
    if ( ! sql.open( f,s ) ) return;

    sql.synchronous(false);

    sql.query(" CREATE TABLE IF NOT EXISTS lddat("
              "   chr       INTEGER UNSIGNED NOT NULL , "
              "   bp        INTEGER UNSIGNED NOT NULL , "
              "   varid     VARCHAR(40) NOT NULL , "
              "   bp1       INTEGER NOT NULL , "
              "   bp2       INTEGER NOT NULL , "
    	      "  PRIMARY KEY(chr,bp) ) ; " );

    sql.query( " CREATE INDEX IF NOT EXISTS i1 ON lddat(varid); " );

    q_insert = sql.prepare( "INSERT OR REPLACE INTO lddat "
                            " (varid,chr,bp,bp1,bp2) "
                            " values ( :varid, :chr, :bp, :bp1, :bp2 ) ; " );

    q_lookup = sql.prepare( "SELECT varid, bp1, bp2 FROM lddat "
                            "WHERE chr == :chr AND bp == :bp ; " );

    q_lookup2 = sql.prepare( "SELECT chr, bp FROM lddat "
                            "WHERE varid == :varid ; " );
}


void LDDB::init_lddb( const std::string & f, const std::string & s, bool newdb )
{
    if ( ! sql.open( f,s ) ) return;

    sql.synchronous(false);

    if( newdb ) {
    	sql.query(" CREATE TABLE IF NOT EXISTS lddat("
              "   chr       INTEGER UNSIGNED NOT NULL , "
              "   bp        INTEGER UNSIGNED NOT NULL , "
              "   varid     VARCHAR(40) NOT NULL , "
              "   bp1       INTEGER NOT NULL , "
              "   bp2       INTEGER NOT NULL , "
              "  PRIMARY KEY(chr,bp) ) ; " );

    	sql.query( " CREATE INDEX IF NOT EXISTS i1 ON lddat(varid); " );
    }

    q_insert = sql.prepare( "INSERT INTO lddat "
                            " (varid,chr,bp,bp1,bp2) "
                            " values ( :varid, :chr, :bp, :bp1, :bp2 ) ; " );

    q_lookup = sql.prepare( "SELECT varid, bp1, bp2 FROM lddat "
                            "WHERE chr == :chr AND bp == :bp ; " );

    q_lookup2 = sql.prepare( "SELECT chr, bp FROM lddat "
                            "WHERE varid == :varid ; " );
}


bool LDDB::query( const int chr , const uint64_t & bp , uint64_t & bp1 , uint64_t & bp2 )
{
    return query(chr,bp,"",bp1,bp2);
}

bool LDDB::query(const int chr, const uint64_t & bp, const std::string & s, uint64_t & bp1, uint64_t & bp2 )
{
    bp1 = bp2 = 0;
    sql.bind_int( q_lookup, ":chr" , chr );
    sql.bind_int64( q_lookup , ":bp" , bp );
    if ( sql.step( q_lookup ) )
    {
        if ( s != "" && s != "." )
	{
            std::string s2 = sql.get_text( q_lookup , 0 );
            if ( s2 != "" && s2 !=s )
	    {
                sql.reset( q_lookup) ;
                return false;
	    }
	}
        bp1 = sql.get_int64( q_lookup , 1 );
        bp2 = sql.get_int64( q_lookup , 2 );
    }
    else
    {
        sql.reset( q_lookup );
        return false;
    }
    sql.reset( q_lookup );
    return true;
}

bool LDDB::query( const std::string & s, int & chr, int & bp )
{
    chr = bp = 0;
    sql.bind_text( q_lookup2, ":varid", s);
    if ( sql.step( q_lookup2 ) )
    {
        chr = sql.get_int( q_lookup2 , 0 );
        bp = sql.get_int( q_lookup2 , 1 );
    }
    else
    {
        sql.reset( q_lookup2 );
        return false;
    }

    sql.reset( q_lookup2 );
    return true;
} 

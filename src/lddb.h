#ifndef __LDDB_H__
#define __LDDB_H__

#include "sqlite3.h"
#include "stdint.h"
#include <string>

class SQL {

 public:

  SQL() { db = NULL; }
  
  bool open(const std::string & n, const std::string & scratch = "" );
  void synchronous(bool);
  void close();
  bool is_open() const { return db; }
  bool query(const std::string & q);

  sqlite3_stmt * prepare(const std::string & q);
  bool step(sqlite3_stmt * stmt);
  void reset( sqlite3_stmt * stmt );
  void finalise(sqlite3_stmt * stmt);

  void begin();
  void commit();

  uint64_t last_insert_rowid()
    { return sqlite3_last_insert_rowid(db); }
  
  void bind_int( sqlite3_stmt * stmt , const std::string & index , int value );
  void bind_int64( sqlite3_stmt * stmt , const std::string & index , uint64_t value );
  void bind_double( sqlite3_stmt * stmt , const std::string & index , double value );
  void bind_text( sqlite3_stmt * stmt , const std::string & index , const std::string & value );
  void bind_null( sqlite3_stmt * stmt , const std::string & index );

  int get_int( sqlite3_stmt *, int );
  uint64_t get_int64( sqlite3_stmt *, int );
  double get_double( sqlite3_stmt *, int );
  std::string get_text( sqlite3_stmt *, int );

  static std::string header_version() 
    { return sqlite3_libversion(); }

  static std::string library_version() 
    { return SQLITE_VERSION; }
  
  sqlite3 * pointer() { return db; }

 private:  
  sqlite3 * db;
  int rc;
  char * db_err; 
};

class LDDB {
 public:
  LDDB();
  LDDB(const std::string & f, const std::string & s = "" );
  void init_lddb(const std::string & f, const std::string & s = "", bool newdb=true );

  bool query( const int , const uint64_t & bp , const std::string & s, uint64_t & bp1 , uint64_t & bp2 );
  bool query( const int , const uint64_t & bp , uint64_t & bp1 , uint64_t & bp2 );
  bool query( const std::string & s, int & chr, int & bp );

protected:
  SQL sql;
  sqlite3_stmt * q_insert;
  sqlite3_stmt * q_lookup;
  sqlite3_stmt * q_lookup2;

};



#endif

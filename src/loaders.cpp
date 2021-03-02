#include "loaders.h"
#include "lddb.h"
#include "helper.h"

#include <iostream>
#include <cstdlib>

//extern Globals globals;

bool is_valid_hum_chr(int chr) {
    if( chr<=0 || (chr>22&&chr<123) || chr>126 ) return false;
    else return true;
}

bool is_valid_chr(int chr, bool _human) {
    if(_human) return is_valid_hum_chr(chr);

    if(chr <=0 ) return false;
    else return true;
}

int get_hum_chr(std::string _chrcode) {
    int chr;
    if ( _chrcode.compare("X")==0  || _chrcode.compare("23")==0 ) chr = 123;
    else if ( _chrcode.compare("Y")==0  || _chrcode.compare("24")==0 ) chr= 124;
    else if ( _chrcode.compare("XY")==0  || _chrcode.compare("25")==0  ) chr= 125;
    else if ( _chrcode.compare("M")==0  || _chrcode.compare("26")==0  ) chr= 126;
    else chr = atoi( _chrcode.c_str() );

    return chr;
}

std::string get_hum_chrcode(int _chr) {
    if(_chr<=22) return int2str(_chr);    
    else if(_chr==123) return "X";
    else if(_chr==124) return "Y";
    else if(_chr==125) return "XY";
    else if(_chr==126) return "M";
    else return "NA";
}

int get_chr(std::string _chrcode, bool _human) {
    if(_human) return get_hum_chr(_chrcode);

    int chr;
    if ( _chrcode.compare("X")==0   ) chr = -1;
    else if ( _chrcode.compare("Y")==0    ) chr= -2;
    else if ( _chrcode == "XY"   ) chr= -3;
    else if ( _chrcode == "M"    ) chr= -4;
    else chr = atoi( _chrcode.c_str() );

    return chr;
}

int load_ranges( const std::string & rangefile, InputData & globals, std::string & err_msg  )
{
    globals.scaffold.clear();
    globals.chr_list.clear();
    globals.total_seq = 0;
    globals.use_range = true;
    globals.use_map = false;

    if ( ! fileExists( rangefile ) ) 
    {
	std::cerr << "could not find range list " << rangefile << "\n";
	exit(1);
    }
    
    std::ifstream MAP( rangefile.c_str() , std::ios::in );
    int i=0;
    int err_line=0;
    while ( !MAP.eof() )
    {	
        std::vector<std::string> lines = tokenizeLine( MAP );

	if(lines.size()==0) break;

	if(lines.size()!=3) {
		err_line++;
		continue;
	}

        int chr = get_chr(lines[0], globals.is_human);
	if ( ! ( 
	   	   ( chr >= 1 && chr <= 22 ) 
	        || ( chr >= -3 && chr <= -1 )   // not sure what encoding PH added here...
		|| ( chr >= 123 && chr <= 126 ) ) ) continue;
	
        int bp1, bp2;
        if ( ! str2int( lines[1] , bp1 ) ) continue;
        if ( ! str2int( lines[2] , bp2 ) ) continue;

        globals.scaffold[ chr ] = Interval( chr , bp1 , bp2 );
        globals.total_seq += (long unsigned int)(bp2 - bp1 + 1);
	globals.chr_list.push_back( chr );
	i++;
    }
    MAP.close();

    if(err_line>0) {
	err_msg = "  warning: 3 columns (chr bp1 bp2) required in " + rangefile + "\n" + 
	          "  " + int2str(err_line) + " lines skipped\n";
    }

    return i;
}

std::vector<Interval> load_snps( const std::string & mapfile, InputData & globals, std::string & err_msg )
{
    std::vector<Interval> snps;
    std::map<Interval, int> prev_snps;
    bool not_ordered = false;
    int x_bp = -1;
    int x_chr = 1;

    std::map<int,int> min;
    std::map<int,int> max;

    if ( ! fileExists( mapfile ) ) 
    {
	std::cerr << "could not find map file " << mapfile << "\n";
	exit(1);
    }

    try {

    std::ifstream MAP( mapfile.c_str() , std::ios::in );
    int cnt = 0;
    int err_line=0;
    while (! MAP.eof() )
    {
        std::vector<std::string> lines = tokenizeLine( MAP );
        if( lines.size()==0 ) break;
        if( lines.size()!=2 ) {
		err_line++;
		continue;;
	}

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp;
        if ( ! str2int( lines[1] , bp ) ) continue;


        // snps.push_back( Interval( chr, bp ) );
        // prev_snps.insert( Interval( chr , bp ), bp );
        prev_snps[ Interval( chr , bp ) ] = bp ;
	
	if ( ! globals.use_range )
	{
	    if ( min.find( chr ) == min.end() )
	    {
		min[chr] = bp;
		max[chr] = bp;
	    }
	    else
	    {
		if ( bp < min[chr] ) min[chr] = bp;
		if ( bp > max[chr] ) max[chr] = bp;
	    }
	}

        x_chr = chr;
        x_bp = bp;

        cnt++;

	if(cnt%100000==0) {
		std::cerr << "Reading " << cnt << " map positions\r";
		std::cout.flush();
	}
    }

    MAP.close();

    not_ordered = true;
    if(not_ordered) {
        // std::cout << "map file not ordered " << cnt << " SNPs \n";

        snps.clear();
        snps.resize(prev_snps.size());

        std::map<Interval, int>::iterator itr = prev_snps.begin();
        int i=0;
        while(itr != prev_snps.end()) {
            snps[i] = itr->first;

            itr++;
            i++;
        }
    }

    // get scaffold, if not already done
    if ( ! globals.use_range )
    {
	globals.total_seq = 0; 
	globals.scaffold.clear();
	std::map<int,int>::iterator i = min.begin();
	while ( i != min.end() )
	{
            //std::cout << "i->first " << i->first << " i->second" << i->second << " max" << max[ i->first ] << "\n";

	    globals.scaffold[ i->first ] = Interval( i->first , i->second , max[ i->first ] );
	    globals.total_seq += (long unsigned int)(max[ i->first ] - i->second + 1);
	    ++i;
	}
    }

    if(err_line>0) {
	err_msg = "  warning: 2 column (chr bp) required in " + mapfile + "\n"; 
    }

    return snps;

    }
    catch(std::bad_alloc a) {
	std::cerr << "Memory Allocation Error: reading map file\n";
	exit(1);
    }

}

std::map<std::string, std::set<Interval> > loader_targets(
        const std::string & targetfile , std::map<std::string,Gene> & genes, 
	InputData & globals, std::string & err_msg, std::string & res_msg  )
{
    
    std::map<std::string, std::set<Interval> > targets;

    if ( ! fileExists( targetfile ) ) 
    {
	std::cerr << "could not find target list " << targetfile << "\n";
	exit(1);
    }

    std::ifstream TARGET( targetfile.c_str() , std::ios::in );

    int cnt=0;
    int err_line=0;
    std::set<std::string> uniq_gene;
    std::set<std::string> uniq_gene_not_found;

    while ( ! TARGET.eof()  )
    {
	std::vector<std::string> lines = tokenizeLine( TARGET );
        if ( lines.size()==0 ) break; 
        if ( lines.size() < 2 ) { 
		err_line++;
		continue;
	}
	
	std::string genename = lines[0];
	std::string setname = lines[1];
        for (unsigned int i=2;i<lines.size(); i++) setname += " " + lines[i];

	if(lines.size()==2) {
		setname += " set_" + lines[1];
	}
	
	std::map<std::string,Gene>::const_iterator i = genes.find( genename );
	
	if ( i == genes.end() ) 
        {
            uniq_gene_not_found.insert( genename );
            continue;
        }
	
	targets[setname].insert( i->second.gene );

	uniq_gene.insert( genename );

	++cnt;
    }
    TARGET.close();

    res_msg = "  read " + int2str( uniq_gene.size()) + " unique genes/targets, in " + 
    		int2str( targets.size()) + " groups\n" + 
		"  " + int2str( cnt) + " total gene/group pairs\n";
	if( uniq_gene_not_found.size()>0 ) {	 
            res_msg = res_msg + "  " + int2str( uniq_gene_not_found.size()) +
	    " genes/targets not found in main gene-list\n";
	}
    
    // apply size filters?
    if ( globals.tsize_min || globals.tsize_max )
    {
	std::map<std::string, std::set<Interval> > targets2 = targets;
	targets.clear();
	std::map<std::string, std::set<Interval> >::iterator i = targets2.begin();
	while ( i != targets2.end() )
	{
	    bool okay = true;
	    if ( globals.tsize_min && i->second.size() < globals.tsize_min ) okay = false;
	    if ( globals.tsize_max && i->second.size() > globals.tsize_max ) okay = false;
	    if ( okay ) targets[ i->first ] = i->second;
	    ++i;
	}
	res_msg += "  after size filters, " + int2str( targets.size()) +  " groups remaining\n";
    }

    if(err_line>0) {
	err_msg = "  warning: at least 2 column (gene_id set_id) required in " + targetfile + "\n" + 
	          "  " + int2str(err_line) + " lines skipped\n";
    }

    return targets;
}



std::map<std::string,Gene> loader_genes( const std::string & genefile ,
                                         std::set<std::string> & bg, InputData & globals, 
					 std::string & err_msg )
{

    globals.largest_gene = 0;

    // CHR BP1 BP2 NAME DESC
    std::map<std::string, Gene > genes;
    if ( ! fileExists( genefile ) )
    {
        std::cerr << "could not find gene list " << genefile << "\n";
        exit(1);
    }

    std::ifstream GENES( genefile.c_str() , std::ios::in );
    int cnt=0;
    int err_line=0;
    while ( ! GENES.eof() )
    {

        std::vector<std::string> lines = tokenizeLine( GENES );

        if ( lines.size()==0 ) break; 
        if ( lines.size() < 4 ) {
		err_line++;
		continue;
	}

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp1, bp2;
        if ( ! str2int( lines[1] , bp1 ) ) continue;
        if ( ! str2int( lines[2] , bp2 ) ) continue;

        std::string name1 = lines[3];

        // are we filtering w/ a background list?
        if ( globals.background_list && bg.find( name1 ) == bg.end() )
            continue;

        // append rest as space-delim description
        std::string name2 = "";
        for (unsigned int i=4; i<lines.size(); i++) name2 += lines[i] + " ";


        bp1 -= globals.window;
        bp2 += globals.window;
        if ( bp1 < 1 ) bp1 = 1;
        if ( bp2 < bp1 ) bp2 = bp1;

        // track largest gene (+window)
        if ( bp2 - bp1 + 1 > globals.largest_gene ) globals.largest_gene = bp2 - bp1 + 1;

        genes.insert( make_pair( name1, Gene( Interval( chr , bp1 , bp2 , name1 ) , name2 ) ) );
        ++cnt;
    }
    GENES.close();

    // get 1 larger than largest gene
    ++globals.largest_gene;

    if(err_line>0) {
	err_msg = "  warning: at least 4 column (chr bp1 bp2 gene_id) required in " + genefile + "\n" + 
	          "  " + int2str(err_line) + " lines skipped\n";
    }

    return genes;
}


std::set<std::string> load_background_genes( const std::string & backgroundfile, std::string & err_msg )
{
    std::set<std::string> bg;

    if ( ! fileExists( backgroundfile ) ) 
    {
	std::cerr << "could not find background gene list " << backgroundfile << "\n";
	exit(1);
    }

    std::ifstream BG( backgroundfile.c_str() , std::ios::in );
    int err_line=0;
    while ( ! BG.eof() )
    {        
        if ( BG.eof() ) break;

        std::vector<std::string> lines = tokenizeLine( BG );
        if(lines.size()==0) break;
	/* if(lines.size()!=1) {
		err_line++; 
		continue;
	} */

        bg.insert(lines[0]);
    }
    BG.close();

    if(err_line>0) {
	err_msg = "  warning: 1 column (gene_id) required in " + backgroundfile + "\n" + 
	          "  " + int2str(err_line) + " lines skipped\n";
    }

    return bg;
}

std::map<std::string, SNP> load_ld_table( const std::string & ldfile, InputData & globals ) {

    std::map<std::string, SNP> snps;

    if ( ! fileExists( ldfile ) )  {
        std::cerr << "could not find ld table " << ldfile << "\n";
        // exit(1);
        return snps;
    }


    // SNP  CHR BP LEFT_TAG RIGHT_TAG
    std::ifstream SNPS( ldfile.c_str() , std::ios::in );
    int cnt = 0;
    while ( !SNPS.eof() )
    {
	std::vector<std::string> lines = tokenizeLine( SNPS );


	if ( lines.size() != 5 ) continue;
	
        std::string name = lines[0];

        int chr = get_chr(lines[1], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;


        int bp, left_tag, right_tag;
	if ( ! str2int( lines[2] , bp ) ) continue;	
	if ( ! str2int( lines[3] , left_tag ) ) continue;
	if ( ! str2int( lines[4] , right_tag ) ) continue;	

	
	snps.insert( make_pair( name, SNP( chr , bp, name, left_tag, right_tag ) ) );
	++cnt;
    }

    SNPS.close();
    std::cerr << "  \tread " << cnt << " snps\n";   
    return snps;    
}

std::set<Interval> load_assoc( const std::string & testfile, InputData & globals, std::string & err_msg  )
{    
    std::set<Interval> orig;
    
    // test regions
    if ( ! fileExists( testfile ) ) 
    {
	std::cerr << "could not find test interval list " << testfile << "\n";
	exit(1);
    }

    std::ifstream TEST( testfile.c_str() , std::ios::in );
    int cnt = 0;
    int err_line = 0;
    while ( !TEST.eof() )
    {
	std::vector<std::string> lines = tokenizeLine( TEST );

	if ( lines.size() == 0 ) break;
	if ( lines.size() != 3 ) { 
		err_line++;
		continue;	
	}

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;


        int bp1 , bp2;
	if ( ! str2int( lines[1] , bp1 ) ) continue;
	if ( ! str2int( lines[2] , bp2 ) ) continue;
	
	orig.insert( Interval( chr , bp1, bp2 , 0 , "int"+int2str(cnt++) ) );
	
	if ( cnt == globals.topn ) break;
    }
    TEST.close();

    if(err_line>0) {
	err_msg = "  warning: 3 columns (chr bp1 bp2) required in " + testfile + "\n" + 
	          "  " + int2str(err_line) + " lines skipped\n";
    }


    return orig;
}

std::set<SNP> load_assoc_snps( const std::string & testfile, InputData & globals , std::string & err_msg )
{
    std::set<SNP> orig;

    // data format
    // CHR BP SNP_ID P

    if ( ! fileExists( testfile ) ) 
    {
	std::cerr << "could not find test SNP list " << testfile << "\n";
	exit(1);
    }

    std::ifstream TEST( testfile.c_str() , std::ios::in );
    int err_line=0;
    int dummy_pos=1;
    while ( ! TEST.eof() )
    {
	std::vector<std::string> lines = tokenizeLine( TEST );
        if( lines.size() == 0 ) break;
	if( lines.size() != 2 && lines.size()!=1 ) { 
		err_line++;
		continue;
	}

	if( lines.size()==1 ) {
		orig.insert( SNP( -1, dummy_pos, lines[0], -1) );
		dummy_pos++;
	}
	else {
        	int chr = get_chr(lines[0], globals.is_human);
        	if( ! is_valid_chr(chr, globals.is_human) ) continue;
		int bp;
        	if ( ! str2int( lines[1] , bp ) ) continue;
        	std::string name = "NA";
		double p = -1;
		orig.insert( SNP( chr , bp, name , p) );
	}

    }
    TEST.close();

    if(err_line>0) {
	err_msg = "  warning: either (chr bp) or (snp_id) required in " + testfile + "\n" + 
	          "  " + int2str(err_line) + " lines skipped\n";
    }

    return orig;
}


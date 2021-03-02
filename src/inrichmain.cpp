#include "inrichmain.h"
#include "crandom.h"
#include "helper.h"
#include "common.h"
#include "limits.h"
#include "loaders.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

InRichMain::InRichMain()
{
    my_log = "";
    err_msg = "";
}

InRichMain::~InRichMain()
{

}


InRichMain::InRichMain(InputData & _globals) 
{
    globals = _globals;
    my_log = "";
}


int InRichMain::load_background_data(){
    background_genes = load_background_genes( globals.backgroundfile, err_msg );
    return background_genes.size();
}

int InRichMain::load_range_data() {
    return load_ranges( globals.rangefile, globals, err_msg );
}

int InRichMain::load_assoc_data() {
    orig.clear();
    orig = load_assoc( globals.testfile, globals, err_msg );
    return orig.size();
}

int InRichMain::load_assoc_snps_data() {
    orig_snp.clear();
    orig_snp = load_assoc_snps( globals.testfile, globals, err_msg );
    return orig_snp.size();
}

int InRichMain::get_closest_ref_snp_pos_2( int _snp_no, int _chr, int _pos, bool _left, int & start_pos ) {
    int x_pos = -1;
    int dist;

    // To retrieve the x_pos;
    //int start = start_pos-1 > 0 ? start_pos-1 : 0;
    if(start_pos>0 && snps[start_pos-1].chr == _chr) {
        x_pos = snps[start_pos-1].bp1;
    }
    for(int i=start_pos; i<_snp_no; i++) {

        start_pos = i;

        if(snps[i].chr < _chr ) {
            continue;
        }
        else if(snps[i].chr > _chr ) {
            break;
        }
        else {
            dist = snps[i].bp1 - _pos ;

	    // find the one right before the left tag SNP
            if(_left && dist>0) {
                start_pos--;
                return x_pos;
            }

	    // left tag SNP 
	    if(_left && dist==0) {
		return snps[i].bp1;
	    }

            if(!_left && dist>=0) {
                return snps[i].bp1;
            }

            x_pos = snps[i].bp1;
        }
    }

    // The last existing SNP precedes the extended tagging SNP region.
    if(!_left) {
        return x_pos;
    }

    return -1;
}


int InRichMain::get_closest_ref_snp_pos( int _snp_no, int _chr, int _pos, bool _left, int & start_pos ) {
    int x_pos = -1;
    int dist;

    // To retrieve the x_pos;
    //int start = start_pos-1 > 0 ? start_pos-1 : 0;
    if(start_pos>0 && snps[start_pos-1].chr == _chr) {
        x_pos = snps[start_pos-1].bp1;
    }
    for(int i=start_pos; i<_snp_no; i++) {

        start_pos = i;

        if(snps[i].chr < _chr ) {
            continue;
        }
        else if(snps[i].chr > _chr ) {
            break;
        }
        else {
            dist = snps[i].bp1 - _pos ;
            if(_left && dist>=0) {
                return snps[i].bp1;
            }
            if(!_left && dist==0) {
                return _pos;
            }
            if(!_left && dist>0) {
                start_pos--;
                return x_pos;
            }

            x_pos = snps[i].bp1;
        }
    }

    // The last existing SNP precedes the extended tagging SNP region.
    if(!_left) {
        return x_pos;
    }

    return -1;
}


int InRichMain::num_in_interval( const Interval & interval ,
                                 const std::set<Interval> & regions )
{

    Interval left( interval.chr, interval.bp1 , interval.bp1 );
    Interval right( interval.chr, interval.bp2 , interval.bp2 + globals.largest_gene );

    std::set<Interval>::iterator lwr = regions.upper_bound( left );
    std::set<Interval>::iterator upr = regions.lower_bound( right );

    int l = interval.bp1 - globals.largest_gene;
    Interval left2( interval.chr, l, l );

    std::set<Interval>::iterator lwr2 = regions.upper_bound( left2 );

    int cnt = 0 ;
    if ( lwr != lwr2 )
    {
        while ( lwr2 != lwr )
        {
            if ( lwr2->bp2 >= interval.bp1 )
            {
                ++cnt;
            }
            ++lwr2;
        }
    }

    while ( lwr != upr )
    {
        ++lwr;
        ++cnt;
    }

    return cnt;
}


std::map<Interval,std::set<Interval> >
        InRichMain::get_overlap( const std::set<Interval> & test ,
                                 const std::set<Interval> & target )
{

    std::set<Interval>::iterator t = test.begin();
    std::set<Interval>::iterator i = target.begin();
    std::set<Interval>::iterator j = target.begin();
    std::map<Interval,std::set<Interval> > overlap;

    int no=-1;
    while ( t != test.end() )
    {

        // for test region t , consider as many target overlaps as possible.
	no++;
	// std::cout << "get_overlap " << no << "\n";

        bool first_overlap = true;

        while (  ( ! ( *t < *j ) )  ||  *t == *j )
        {

            if ( *t == *j ) // overlaps?
            {
                overlap[ *t ].insert( *j );
                if ( first_overlap )
                {
                    i = j;
                    first_overlap = false;
                }

            }

            //consider next target, if one exists
            ++j;

            if ( j == target.end() ) break;
        }


        // reset to last possible overlap
        if ( ! first_overlap )
        {
            j = i;
        }
        else
            if ( j == target.end() ) break;

        // consider next test region
        ++t;
    }

    return overlap;
}



std::set<Interval> InRichMain::get_set( std::map< std::string,Gene> & g )
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

void InRichMain::cal_study_sig_hit() {
    try {

        //log( "  initiating study hit test ") ;
        org_study_hit.clear();
        p_sig_study_hit.clear();
        org_study_hit.resize( N_SIG_STUDY_HIT, 0 );
        p_sig_study_hit.resize( N_SIG_STUDY_HIT, 0 );

        org_study_gene_hit.clear();
        p_sig_study_gene_hit.clear();
        org_study_gene_hit.resize( N_SIG_STUDY_HIT, 0 );
        p_sig_study_gene_hit.resize( N_SIG_STUDY_HIT, 0 );


        std::map<int, std::set<Interval> > uniq_genes;
        for( unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
            std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
            int pway = 0;
            while ( t != targets.end() ) {
                // count the number of pathways with p-value less than a threshold
                if((double)(pcount[ pway ])/(double)(globals.nrep+1) < SIG_STUDY_HIT[ns]) {
                    org_study_hit[ns]++;                    
                    get_overlapping_genes( test, t->second, uniq_genes[ns]);
                }
                t++;
                pway++;
            }

            org_study_gene_hit[ns] = uniq_genes[ns].size();
        }

	//log( " done\n" );
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRich::cal_study_sig_hit(" + int2str(N_SIG_STUDY_HIT) + ")\n";
        exit(1);
    }
}

void InRichMain::log_main_file() {
    main_filename =  CRandom::rand_letters(filename_size) + ".main";

    std::ofstream OUT1( (main_filename).c_str() , std::ios::out );
    log_main_file(OUT1);
    OUT1.close();
}

void InRichMain::log_interval_file() {
    interval_filename =  CRandom::rand_letters(filename_size)+ ".intervals";

    std::ofstream OUT2( (interval_filename).c_str() , std::ios::out );
    log_interval_file(OUT2);
    OUT2.close();
}

void InRichMain::log_target_file(){
    target_filename =  CRandom::rand_letters(filename_size)+ ".targets";

    std::ofstream OUT3( (target_filename).c_str() , std::ios::out );
    log_target_file(OUT3);
    OUT3.close();
}
void InRichMain::log_gene_file(){
    gene_filename =  CRandom::rand_letters(filename_size)+ ".genes";

    std::ofstream OUT4( (gene_filename).c_str() , std::ios::out );
    log_gene_file(OUT4);
    OUT4.close();
}

std::string InRichMain::header() {

    
    std::ostringstream mylog;

    mylog << "\n----------------------------------------------------------------\n"
            << std::setw(6) << "T_Size"
            << std::setw(5) << ( globals.test == INTERVAL ? "Int" : "Tar" ) << "_No"
            << std::setw(15) << "Empirical_P"
            << std::setw(12) << "Corrected_P"
            << std::setw(12) << "Target\n";

    mylog << "----------------------------------------------------------------\n";

    return mylog.str();
}

void InRichMain::log_target_file(std::ofstream & OUT1){
    /*    OUT1 << std::setw(6) << "T_TARG"
            << std::setw(6) << ( globals.test == INTERVAL ? "N_INT" : "N_TARG" )
            << std::setw(12) << "P" 
            << std::setw(12) << "PCORR" 
            << std::setw(6) << "LOW_WARN" 
            << std::setw(12) << "UNSH" 
            << std::setw(4) << "" << "TARGET\n"; 
*/
    OUT1 << "_O1\t"
            << "T_TARG" << "\t"
            << ( globals.test == INTERVAL ? "N_INT" : "N_TARG" ) << "\t"
            << "P" << "\t"
            << "PCORR" << "\t"
            << "LOW_WARN" << "\t"
            << "UNSH" << "\t"
            << "TARGET" << "\t"
            << "EXP_" << ( globals.test == INTERVAL ? "N_INT" : "N_TARG" ) << "\n";

    for(unsigned int i=0; i<path_data.size(); i++) {
        OUT1 << "_O1\t" 
                << path_data.at(i).get_inrich_summary();
    }
}

bool InRichMain::load_files() {
    int file_read = 0;
    int file_tried = 0;

    if(globals.convert_mode) {
	// map file is requisite
    	file_tried++;
    	if(load_file(SNP_MAP)){
      	    file_read++;
    	}

        file_tried++;
        if(load_file(ASSOC_SNP)){
            file_read++;
        }
    }
    else {
    	if(globals.use_map) {
            file_tried++;
            if(load_file(SNP_MAP)){
                file_read++;
            }
    	}
        if(globals.background_list) {
            file_tried++;
            if(load_file(BACKGROUND)) {
                file_read++;
            }
        }
        if(globals.use_range) {
            file_tried++;
            if(load_file(RANGE)) {
                file_read++;
            }
        }

        file_tried++;
        if(load_file(GENE)){
            file_read++;
        }        
        if(!globals.exam_bp_dist) {
            file_tried++;
            if(load_file(TARGET_SET)){
                file_read++;
            }
        }
        file_tried++;
        if(load_file(ASSOC)){
            file_read++;
        }
    }

    if( file_read == file_tried ) {
        if( globals.convert_mode && file_read==2 ) return true;
        else return is_data_ok();
    }
    else {
        return false;
    }
}

void InRichMain::update_all_genes_map() {
    // ------------------------------------------------------------
    // Phil H. Lee 2010-12-15 Added
    // all_genes_set have been merged,
    // but all_genes still includes not-merged gene information
    // ------------------------------------------------------------
    std::set<Interval>::iterator itr = all_genes_set.begin();
    while(itr!=all_genes_set.end()) {
        // check whether this is a merged gene
        if(itr->name.find(",")!=std::string::npos) {
            // parse individual gene
            std::vector<std::string> gene_ids = tokenize_string(itr->name, ',');
            std::string merged_gene_symbol = "";
            std::string merged_gene_desc = "";
            for(unsigned int i=0; i<gene_ids.size(); i++) {
                // extract only the gene symbol
                std::string symbol = extract_symbol(all_genes[gene_ids[i]].desc);
                std::string desc = extract_desc(all_genes[gene_ids[i]].desc);
                if(i>0) {
                    merged_gene_symbol = merged_gene_symbol + ",";
                    merged_gene_desc = merged_gene_desc + "; ";
                }
                merged_gene_symbol = merged_gene_symbol + symbol;
                merged_gene_desc = merged_gene_desc + desc ;
            }
            merged_gene_desc = merged_gene_desc + "; ";

            all_genes.insert( make_pair( itr->name, Gene( Interval( itr->chr , itr->bp1 , itr->bp2 , itr->name ) , merged_gene_symbol + " " + merged_gene_desc ) ) );
        }

        itr++;
    }
}

std::set<Interval> InRichMain::gen_match_set(std::map< Interval , std::set<Interval> > & overlap1,
                               std::map< Interval , std::set<Interval> > & overlap2,
			       std::set<Interval> & all_target_genes,
			       std::set<Interval> & match ) {

    std::set<Interval> new_match;
    new_match.clear();

    std::set<Interval>::iterator itr = match.begin();
    int i=0;
    while ( itr != match.end() )
    {

        // Is this in the overlap list?
        std::map< Interval , std::set<Interval> >::iterator j = overlap1.find( *itr );


        // No genes overlap
        if ( j == overlap1.end() )
            new_match.insert( *itr );
        else
        {

            // for each target overlapping this gene, see if target overlaps any other gene?
            std::set<Interval> & tgt = j->second;
            std::set<Interval>::iterator s = tgt.begin();
            while ( s != tgt.end() )
            {
                std::set<Interval> & genes = overlap2[ *s ];
                if ( genes.size() == 0 )
                {
                    std::cerr << "internal logic error\n";
                    exit(1);
                }

                // okay, gene -> target -> only one gene
                if ( genes.size() == 1 )
                {
                    new_match.insert( *itr );
                }
                else
                {
                    Interval m = *genes.begin();
                    std::set<Interval>::iterator k = genes.begin();
                    while ( k != genes.end() )
                    {
                        m.combine( *k );
                        ++k;
                    }

                    new_match.insert( m );

                }
                ++s;
            }
        }
        ++itr;
        ++i;
    }

    return new_match;
}


void InRichMain::gen_test_set(std::map< Interval , std::set<Interval> > & overlap1,
                              std::map< Interval , std::set<Interval> > & overlap2) {

    std::set<Interval>::iterator itr = orig.begin();
    int i=0;
    while ( itr != orig.end() )
    {

        // Is this in the overlap list?
        std::map< Interval , std::set<Interval> >::iterator j = overlap1.find( *itr );


        // No genes overlap
        if ( j == overlap1.end() )
            test.insert( *itr );
        else
        {

            // for each target overlapping this gene, see if target overlaps any other gene?
            std::set<Interval> & tgt = j->second;
            std::set<Interval>::iterator s = tgt.begin();
            while ( s != tgt.end() )
            {
                std::set<Interval> & genes = overlap2[ *s ];
                if ( genes.size() == 0 )
                {
                    std::cerr << "internal logic error\n";
                    exit(1);
                }

                // okay, gene -> target -> only one gene
                if ( genes.size() == 1 )
                {
                    test.insert( *itr );
                }
                else
                {
                    Interval m = *genes.begin();
                    std::set<Interval>::iterator k = genes.begin();
                    while ( k != genes.end() )
                    {
                        m.combine( *k );
                        ++k;
                    }

                    test.insert( m );

                }
                ++s;
            }
        }
        ++itr;
        ++i;
    }
}


bool InRichMain::merge_data() {
    Log mylog( globals.outroot + ".out.inrich" );

    if( globals.match_genes ) {
        int drop_no = 0;

	if (globals.compact) {
	  drop_no = orig.size();
	  orig = make_compact( drop_no, orig, all_genes_set );
	  drop_no = drop_no - orig.size();
	}

        mylog << "  "
                << orig.size()
                << " test elements on gene regions\n"
                << "  "  << drop_no << " non-genic test elements dropped\n";
    }
    all_genes_set = merge( all_genes_set );
    update_all_genes_map();

    all_target_genes = get_all_targets();
    // std::set<Interval> all_target_genes = get_all_targets();


    //
    // Merge any test regions that span the same target (from *any included
    // set*), or with self, and place in test[].  This way, we are only
    // dealing with a single set of test intervals downstream (and so can
    // permute only once per replicate).
    //
    orig = merge( orig );

    std::map< Interval , std::set<Interval> > overlap2 = get_overlap( all_target_genes , orig );
    std::map< Interval , std::set<Interval> > overlap1 = get_overlap( orig , all_target_genes );

    test.clear();
    gen_test_set(overlap1, overlap2);
    test = remove_subinterval( test );
    test = merge( test );

    n_interval = test.size();
    n_set = targets.size();

    mylog << "\nAfter merging, " << all_genes_set.size() <<  " non-overlapping reference genes\n";
    mylog << "After merging, " <<  all_target_genes.size() <<  " non-overlapping genes in target sets\n";
    mylog << "After merging, " <<  n_interval <<  " non-overlapping intervals\n\n";

    return true;
}

bool InRichMain::load_file(int _mode) {

    int no=-1;
    int no1=-1;

    std::string res_msg;
    err_msg = "";

    switch(_mode) {
    case BACKGROUND:
        no = load_background_data();   // ok
	log( "  read " + int2str(no) + " background genes/targets\n");
        break;

    case RANGE:
        no = load_range_data();  // ok
	log( double2str((double)globals.total_seq/(double)1000000) +  "Mb total sequence length on " + int2str(globals.scaffold.size()) + " chromosomes\n");
	globals.precompute = false;	
        break;

    case SNP_MAP:
        snps.clear();
        snps = load_snps( globals.mapfile, globals, err_msg );   // ok
        no = snps.size();

        log( double2str(globals.total_seq/(double)1e6) + "Mb total sequence length on "
             + int2str(globals.scaffold.size()) + " chromosomes\n");
        log( "  read " + int2str(snps.size()) + " map positions\n");

        break;

    case GENE:
        all_genes.clear();   // ok
        all_genes = loader_genes( globals.genefile, background_genes, globals, err_msg );
        no = all_genes.size();
        log( "  read " + int2str(no) + " reference genes\n");
        all_genes_set = get_set( all_genes );
        break;

    case TARGET_SET:
        targets.clear();
        targets = loader_targets( globals.targetfile, all_genes, globals, err_msg, res_msg );
        no = targets.size();
	log( res_msg );
        break;

    case ASSOC:
        no = load_assoc_data();   // ok
        log( "  read " + int2str(no) + " intervals\n");

	if (globals.use_map) {
		no1 = filter_non_map_data();
		if( no1 > 0 ) {
			log( "  removed " + int2str(no1) + " intervals due to missing map information\n");
		}
	}

        break;

    case ASSOC_SNP:
        no = load_assoc_snps_data();   // ok
        log( "  read " + int2str(no) + " SNPs\n\n");

        break;

    default:
        break;
    }

    if(err_msg.size() > 1) {
	log ( err_msg );
    }

    if(no>0) {
        return true;
    }
    else {
        return false;
    }

}

int InRichMain::filter_non_map_data() {

    // identify min max index for each CHR
    std::map<int, int> min;
    std::map<int, int> max;
    int x_chr=0;

    int snp_no = snps.size(); 
    int i=0;
    for(i=0; i<snp_no; i++ ) {
        int chr = snps[i].chr;  
	if (chr==x_chr)  continue; 

	if( min.find( chr ) == min.end() ) {
		min[chr] = i;
		max[chr] = i;
		// std::cout << "min[" << chr << "] = " << i << "\n";
	}
	else {
		std::cout << "filter_non_map : ERROR\n";
	}

	if( i>0  ) {
		max[x_chr] = i-1; 		
		// std::cout << "max[" << x_chr << "] = " << (i-1) << "\n";
	}

	x_chr = chr;
    }
    max[x_chr] = i;

    // for each interval, check the existence of start bp.
    int not_found_no = 0;
    std::set<Interval> new_orig; 
    std::set<Interval>::iterator itr = orig.begin();
    while( itr != orig.end() ) {
    	if( min.find(itr->chr) == min.end()) {
	        std::string chr = get_hum_chrcode( itr->chr );
		log( "  * WARNING: no map info of chr" + chr + " for interval chr" + chr + ":" + int2str(itr->bp1) + ".." + int2str(itr->bp2) +  "\n"); 
		not_found_no++;
		itr++;
		continue;
	}
	
        bool found = bi_search( min[itr->chr], max[itr->chr], itr->bp1, itr->bp2 );
        if(!found) {
	        std::string chr = get_hum_chrcode( itr->chr );
                log( "  * WARNING: no SNP map info exists for interval chr" + chr + ":" + int2str(itr->bp1) + ".." + int2str(itr->bp2) + "\n");
		not_found_no++;
	}
	else {
		new_orig.insert(*itr);
	}

	itr++;
    }

    if(not_found_no>0) {
		orig = new_orig;
    }

    return not_found_no;
}

bool InRichMain::bi_search( int min_ix, int max_ix, int bp1, int bp2 ) {
    if(min_ix > max_ix) {
	return false;
    }

    int mid_ix = (min_ix + max_ix)/2;


    if( bp1 <= snps[mid_ix].bp1 && snps[mid_ix].bp1 <= bp2 ) {
	return true; 
    }
    else if( snps[mid_ix].bp1 > bp2 ) {
        return bi_search( min_ix, mid_ix-1, bp1, bp2 );
    }
    else {
        return bi_search( mid_ix+1, max_ix, bp1, bp2 );
    }
}


bool InRichMain::init_optional_load() {

    if(globals.exam_bp_dist) {
        log("\nInitiating positional clustering test ...  ");
        test_dist_interval = compute_dist_interval( test );
        val_dist_p_no = (int)test_dist_interval.size();
        test_dist_p = init_test_p( val_dist_p_no );

        log("DONE\n\n");
    }


    // We now have a single set of intervals, with all overlap taken care
    // of. Now permute: we need to do this twice, in order to score
    // multiple-test corrections

    //
    // If using a MAP, assign SNP counts, and check that each falls on a MAP location
    //

    if ( globals.use_map ) {
        if(!assign_snp_count( test , snps )) {
            return false;
        }
    }

    // Precompute okay positions, if using a MAP, to save time later
    // for tricky regions. This will populate the acceptable_positions()
    // value

    if ( globals.precompute )
    {
        if ( ! globals.use_map )
        {
            std::cerr << "no -m map specified, for precompute\n";
            exit(1);
        }

        if(!pre_compute( test, snps , globals.match_genes ? &all_genes_set : NULL )) {
            return false;
        }
    }

    return true;
}


void InRichMain::log(std::string _str1) {
    Log mylog(globals.outroot + ".out.inrich");
    mylog << _str1;
}

//std::map<Interval,std::vector<int> >
bool InRichMain::pre_compute( const std::set<Interval> & t ,
                              const std::vector<Interval> & pos ,
                              std::set<Interval> * genes )
{

    //std::map<Interval,std::vector<int> > m;

    // total number of SNPs in MAP
    const int s = pos.size();
    std::set<Interval>::iterator i = t.begin();
    int i_no=0;
    while ( i != t.end() )
    {

        //
        // need to match i->n SNPs
        //

        int distance = i->bp2 - i->bp1 + 1;
        double lower = distance - distance * globals.tol;
        double upper = distance + distance * globals.tol;
        if ( lower < 1 ) lower = 1;

        //
        // need to match # of genes?
        //

        int num_gene = 0;
        if ( globals.match_genes )
        {
            num_gene = num_in_interval( *i, *genes );
        }


        //
        // SNP numbers for those that match
        //

        std::vector<int> matches;

        //
        // look for all valid sets
        //

        const int s2 = s - i->n + 1;
        for (int p=0; p < s2 ; p++)
        {

            bool density_okay = ( ! globals.use_map ) || i->n == 0 ;
            bool genenum_okay = ! globals.match_genes;

            int chr=-1, bp1=-1, bp2=-1;

            // get the implied distance, moving N SNPs forwards;
            // if this spans a chromosome, we will get negative distance

            int faked;

            if(!density_okay) {
                faked = pos[ p + i->n - 1 ].bp1  - pos[p].bp1 + 1;

                if ( faked >= lower && faked <= upper )
                {
                    density_okay = true;
                    chr = pos[p].chr;
                    bp1 = pos[p].bp1;
                    bp2 = pos[ p + i->n - 1 ].bp1;
                }
            }

            //
            // Did we get a match?
            //

            if ( globals.match_genes && density_okay )
            {
                int g = num_in_interval( Interval( chr, bp1, bp2 ) , *genes );
                if ( g == num_gene ) genenum_okay = true;
            }


            //
            // Is this okay?
            //

            if ( density_okay && genenum_okay ) matches.push_back(p);


        } // next SNP position

        // If a v. large number of possible matches, then no point in
        // storing, as will be easy to find another match, and this will
        // just take up lots of memory.  Flag this by size() == 0 match,
        // i.e.  not even including self, will indicate to the match()
        // function that it doesn't have a precomputed list of positions
        // for this interval. But for now, store everything.

        globals.acceptable_positions[*i] = matches;

        //       std::cout << *i << "\t" << matches.size() << "\t"
        //  		<< (double)matches.size()/(double)s2 << "\n";

        ++i;

        //if(i_no%step==0) {
        std::cout << "Precomputing acceptable positions " << i_no+1 << "\r";
        std::cout.flush();
        //}

        ++i_no;


    } // next test interval

    log("Precomputing acceptable positions " + int2str(t.size()) + "\n");

    return true;
}


bool InRichMain::assign_snp_count(  std::set<Interval> & t ,
                                    const std::vector<Interval> & pos )
{
  
  
  const int s = pos.size();
  
  std::map<Interval,bool> seen_start;
  std::map<Interval,bool> seen_end;
  
  std::map<Interval,int> cnt;
  
  for (int p=0; p < s ; p++)
    {
      int chr = pos[p].chr;
      int bp  = pos[p].bp1;
      
      std::set<Interval>::iterator i = t.begin();
      while ( i != t.end() )
        {
	  if ( i->chr == chr )
            {
	      if ( i->bp1 == bp )
                {
		  cnt[*i]++;
		  seen_start[*i] = true;
		  // if start == end, only count 1 marker
		  if ( i->bp2 == bp ) { seen_end[*i] = true; }
                }
	      else if ( i->bp2 == bp ) { cnt[*i]++; seen_end[*i] = true; }
	      else if ( i->bp1 < bp && i->bp2 > bp ) cnt[*i]++;
            }
	  ++i;
        }
      
      if(p%10000==0) {
	std::cerr << p << " SNP counts assigned\r";
	std::cout.flush();
      }
    }
  
  log( int2str(s) +  " SNP counts assigned\n");
  
  std::set<Interval> copy_t = t;
  t.clear();
  
  std::set<Interval>::iterator i = copy_t.begin();
  while ( i != copy_t.end() )
    {
      
      if ( seen_start.find( *i ) == seen_start.end() )
        std::cerr << "did not observe EXACT start marker for interval " << *i << "\n";
      
      if ( seen_end.find( *i ) == seen_end.end() )
	std::cerr << "did not observe EXACT end marker for interval " << *i << "\n";
      
      Interval j = *i;
      
      std::map<Interval,int>::iterator ctr = cnt.find( *i );

      if( ctr != cnt.end() ) 
	{
	  j.n = ctr->second;
	}
      else 
	{
	  j.n = 0;
	  std::cerr << "no snp marker for interval " << j.name << "\n";
	}
      
      t.insert(j);
      
      ++i;
    }
  
  return true;
  
}

void InRichMain::init_data() {
    count.clear();
    actual_count.clear();
    null_tracker.clear();
    match_tracker.clear();
    pcount.clear();
    pcorr.clear();
    p_stat_sum.clear();

    org_study_hit.clear();
    p_sig_study_hit.clear();
}


void InRichMain::inrich_test() {
    my_log = "";

    init_data();

    // ------------------------------------------------------
    // 1. Load Files
    // ------------------------------------------------------
    if(!load_files()) {
        std::cout << "\n --- halting --- no input data to test : interval/genes/targets ---\n";
        return;
    }

    // ------------------------------------------------------
    // 2. Merge Data
    // ------------------------------------------------------
    if(!merge_data()) {
        log("ERROR: Data not merged");
        return;
    }


    // ------------------------------------------------------
    // 3. Optional Checking
    // ------------------------------------------------------
    if(!init_optional_load()) {
        log("ERROR: Optional Data not Loaded");
        return;
    }

    // ------------------------------------------------------
    // 4. Permutation
    // ------------------------------------------------------
    init_first_permutation();
    if(!run_first_permutation()) {
        log("ERROR: First permutation not completed");
        return;
    }


    init_second_permutation();
    if(!run_second_permutation()) {
        log("ERROR: Second permutation not completed");
        return;
    }

    // ------------------------------------------------------
    // 5. File Output
    // ------------------------------------------------------
    if(!write_output()) {
        log("ERROR: File writing not completed");
        return;
    }
}



void InRichMain::clustering_test() {
    my_log = "";

    init_data();

    // ------------------------------------------------------
    // 1. Load Files
    // ------------------------------------------------------
    if(!load_files()) {
        std::cout << "\n --- halting --- no input data to test : interval/genes ---\n";
        return;
    }


    // ------------------------------------------------------
    // 2. Drop non-genic intervals 
    // ------------------------------------------------------
    globals.match_genes = true;
    if( globals.match_genes ) {
        int drop_no = 0;

	if (globals.compact) {
	  drop_no = orig.size();
	  orig = make_compact( drop_no, orig, all_genes_set );
	  drop_no = drop_no - orig.size();
	}

	log ( "  " +  int2str(orig.size()) + " test elements on gene regions\n" + 
	      int2str(drop_no) +  " non-genic test elements dropped\n") ;
    }

    test = merge( orig );
    log( "\nAfter merging, " + int2str(test.size()) +  " non-overlapping intervals\n" ) ;
    if(test.size()<2) {
        std::cout << "\n --- halting --- no input data to test : interval ---\n";
        return;
    }
    if(test.size()>500) {
        std::cout << "\n --- halting --- too many interval data to test  ---\n";
        std::cout << "\n This test is for examining clustering of top association regions, typically selected with a p-value < 0.001 ---\n";
        return;
    }

    // ------------------------------------------------------
    // 2. Optional Checking
    // ------------------------------------------------------
    if(!init_optional_load()) {
        log("ERROR: Optional Data not Loaded");
        return;
    }

    // ------------------------------------------------------
    // 3. Permutation
    // ------------------------------------------------------
    if(!run_clustering_permutation()) {
        log("ERROR: First permutation not completed");
        return;
    }

    // ------------------------------------------------------
    // 5. File Output
    // ------------------------------------------------------
    report_clustering();
    std::cout << "\nWriting output to " + globals.outroot + ".out.inrich ...\n\n";
}

bool InRichMain::run_clustering_permutation () {

    int p=0;

    try {
        CRandom::srand( globals.seed ) ;

        for ( p = 0 ; p < globals.nrep ; p++ )
        {
            //
            // create a matched set of regions
            //
            std::set< Interval > matched = match( test ,
                                                  snps ,
                                                  match_tracker,
                                                  globals.match_genes ? &all_genes_set : NULL  );


            std::set<PairIntervalDist> random_dist_interval = compute_dist_interval( matched );
            val_dist_p_no = update_dist_p( val_dist_p_no, test_dist_p, test_dist_interval, random_dist_interval );

            std::cerr << p << " permutations         \r";
            std::cout.flush();
        }

        log( int2str(globals.nrep) + " permutations ( completed )\n");

    }
    catch(std::exception a) {
        std::cerr << "Exception: InRich::run_clustering_permutation(" << a.what() << ")\n";
        std::cerr << "rep=" << p << "\n\n";
        exit(1);
    }

    return true;
}



bool InRichMain::run_first_permutation () {

    int p=0;

    try {
        CRandom::srand( globals.seed ) ;

        std::ofstream* perms = NULL;
        if (globals.printPermutations) {
        	std::string permsFile = globals.outroot + ".out.inrich.perms";
        	perms = new std::ofstream(permsFile.c_str(), std::ios::out | std::ios::trunc);

        	*perms << "ORIG" << "\t" << 0;
        	for (std::set<Interval>::const_iterator intIt = test.begin(); intIt != test.end(); ++intIt)
        		*perms << "\t" << "chr" << get_hum_chrcode(intIt->chr) << ":" << intIt->bp1 << ".." << intIt->bp2;
        	*perms << std::endl;
        }

        for ( p = 0 ; p < globals.nrep ; p++ )
        {

            //
            // create a matched set of regions
            //
            std::cerr << p << " first-pass permutations started        \r";

            std::set< Interval > matched = match( test,
                                                  snps,
                                                  match_tracker,
                                                  globals.match_genes ? &all_genes_set : NULL  );

            if (globals.printPermutations) {
            	*perms << "PERM" << "\t" << (p+1);
            	for (std::set<Interval>::const_iterator intIt = matched.begin(); intIt != matched.end(); ++intIt)
            		*perms << "\t" << "chr" << get_hum_chrcode(intIt->chr) << ":" << intIt->bp1 << ".." << intIt->bp2;
            	if (matched.size() < test.size()) {
            		unsigned int numMissing = test.size() - matched.size();
            		for (unsigned int i = 0; i < numMissing; ++i)
            			*perms << "\t" << ".";
            	}
            	*perms << std::endl;
            }

            // =====================================================
            // PL 2010.08.30 ADDED
            // =====================================================
            if(globals.exam_bp_dist) {
                std::set<PairIntervalDist> random_dist_interval = compute_dist_interval( matched );
                val_dist_p_no = update_dist_p( val_dist_p_no, test_dist_p, test_dist_interval, random_dist_interval );
            }

            //
            // repeat test for each set
            //

            std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
            int pway = 0;

            while ( t != targets.end() )
            {
                int perm_statistic = evaluate( matched , t->second );
                if ( perm_statistic >= count[ pway ] ) ++pcount[ pway ];
                null_tracker[ pway ][ perm_statistic ]++;
                p_stat_sum[ pway ] += perm_statistic;
                ++t;
                ++pway;
            }

            std::cerr << p << " first-pass permutations         \r";
            std::cout.flush();
        }

        if (globals.printPermutations) {
        	delete perms;
        	perms = NULL;
        }

        log( int2str(globals.nrep) + " first-pass permutations ( completed )\n");

    }
    catch(std::exception a) {
        std::cerr << "Exception: InRich::run_first_permutation(" << a.what() << ")\n";
	std::cerr << "rep=" << p << "\n\n";
        exit(1);
    }

    return true;
}


void InRichMain::cal_unplaced() {
    double total_unplaced = 0;
    std::map<Interval,int>::iterator ii = match_tracker.begin();
    while ( ii != match_tracker.end() )
    {
        total_unplaced += ii->second;
        ++ii;
    }

    total_unplaced /= (double)( test.size() * globals.nrep );
    log( "Proportion unplaced/self-placed interval/permutations = " +
         double2str(total_unplaced) + "\n\n") ;
}

std::string InRichMain::get_gene_symbol( Interval & _int ) {
    std::string symbol = all_genes[_int.name].desc;
    std::string desc = "";

    if(symbol.compare("")==0) {
        std::vector<std::string> id_list = tokenize_string(_int.name, ',');
        for(unsigned int i=0; i<id_list.size(); i++) {
            std::vector<std::string> desc1 = tokenize_string(all_genes[ id_list[i] ].desc, ' ');

            if(desc1.size()>0) {
                symbol = symbol + desc1[0];
            }
            for(unsigned int j=1; j<desc1.size(); j++) {
                desc = desc + desc1[j];
            }


            if(i<id_list.size()-1) {
                symbol = symbol + ",";
                desc = desc + ",";
            }
        }

        all_genes.insert( make_pair( _int.name, Gene( Interval( _int.chr , _int.bp1 , _int.bp2 , _int.name ) , symbol + " " + desc ) ) );
    }


    return symbol + " " + desc;
}


bool InRichMain::log_target_path_data(){

    //
    // Verbose report of region stats:
    //
    if(main_data.size()>0) {
    	main_data.clear();
    }
    if(path_data.size()>0) {
    	path_data.clear();
    }

    std::map<std::string, std::set<Interval> >::iterator t = targets.begin();

    int pway = 0;
    log(header());
    bool any_results = false;

    std::ostringstream myLog;
    while ( t != targets.end() )
    {
        const std::set<Interval> & target = merge( t->second );

        // For this individual set, proportion of unplaced intervals
        //
        //
        // Get all genes in this target set that were in the test interval
        //
        std::map< Interval , std::set<Interval> > overlap = get_overlap( test, target );
        std::map< Interval , std::set<Interval> >::iterator ii = overlap.begin();

        int unplaced = 0;
        //
        // Flag if pointwise empirical p-value is at minimum bound as the
        // corrected value might be an under-estimate in this case
        //

        double pvalue = pcount[ pway ]/(double)( 1 + globals.nrep );
        bool lowest_p = pcount[ pway ] == 1;
        double meanPermStat = p_stat_sum[ pway ] / (double) globals.nrep;

        // Write to main output file
        double pvalue2 = pcorr[ pway ] / double( 1 + globals.nboot );

        while ( ii != overlap.end() )
        {
            unplaced += match_tracker[ ii->first ];

            std::set<Interval> & o1 = ii->second;
            std::set<Interval>::iterator oo = o1.begin();
            while ( oo != o1.end() )
            {
                Interval tmp = ii->first;
                Interval tmp2 = *oo;

	        std::string symbol = get_gene_symbol(tmp2);	
                main_data.push_back(AllRes(t->first,
                                           double2str(pvalue),
                                           tmp2.name,
                                           tmp2.get_desc(),
                                           symbol,
                                           tmp.get_desc()));

                ++oo;
            }
            ++ii;
        }


        const double proportion_unplaced = overlap.size() > 0 ?
                                           (double)unplaced / (double)( overlap.size() * globals.nrep ) : 0 ;


        path_data.push_back( PathwayRes ( target.size(),
                                          actual_count[ pway ],
                                          pvalue ,
                                          pvalue2 ,
                                          lowest_p ,
                                          proportion_unplaced ,
                                          t->first ,
                                          meanPermStat) );

        if ( pvalue <= globals.pthresh )
        {
            any_results = true;

            std::string mark = " ";
            if (lowest_p) {
                mark = "*";
            }
            myLog << std::setw(6) << target.size() 
                    << std::setw(8) << actual_count[ pway ]
                    << std::setw(15) << pvalue
                    << std::setw(12) << pvalue2 << mark
                    << std::setw(4) << " "
                    << t->first <<  "\n";

        }
        //
        // Output for next target group
        //


        ++t;
        ++pway;
    }

    if ( ! any_results ) {
        log( "     { -- no significant results at p < " + double2str(globals.pthresh) +  " -- }\n") ;
    }
    else {
    	log(myLog.str());
    }


    return true;
}

/*
void InRichMain::update_chr_bp_orig_snp( const std::string dbname ) {
    //init_lddb( dbname, "" );
    log( "  attached LD database " + dbname + "\n");

    std::set<SNP> new_orig_snp;
    std::set<SNP>::iterator itr = orig_snp.begin();
    bool found;
    int no=0;
    while( itr != orig_snp.end() ) {
        if(itr->chr != -1) return ;

	int chr, bp;
	found = query( itr->name, chr, bp );
	if(!found) {
            log( "  ERROR: SNP " + itr->name + " not in LD database\n") ;
            itr++;
            continue;
	}

	no++;	
	new_orig_snp.insert( SNP( chr , bp, itr->name , itr->p) );
	itr++;

        //std::cerr << itr->name << " " << chr << ":" << bp << "\n";
	if(no%10==0) {
            std::cerr << "  " << no << " SNP position retrieved\r";
            std::cout.flush();
	}
    }

    log( "  " + int2str(new_orig_snp.size()) + " SNP position retrieved\n" );
    orig_snp = new_orig_snp;

    return ;
}

bool InRichMain::snp2interval() {
    // ------------------------------------------------------
    // 1. Load Files
    // ------------------------------------------------------
    snps.clear();
    if(!load_files()) {
        log("File not loaded -- stopping");
        return false;
    }

    // ------------------------------------------------------
    // 2. Create DB
    // ------------------------------------------------------
    if(globals.create_lddb) {
        if(!make_lddb( globals.create_db_name, globals.create_db_name + ".db" )) {
            log("DB Creation not completed\n");
            return false;
        }
        globals.use_db_name = globals.create_db_name + ".db";
    }
    else {
    	init_lddb( globals.use_db_name, "", false );
    }


    log( "Converting SNPs to intervals\n");

    // ------------------------------------------------------
    // 2-2. Update chr, bp position if not available  
    // ------------------------------------------------------
    update_chr_bp_orig_snp( globals.use_db_name ); 


    // ------------------------------------------------------
    // 3. Generate and Write Intervals
    // ------------------------------------------------------
    if(globals.test == TARGET) {
    	ld_clumping( globals.use_db_name, orig_snp );
    }
    else {
    	ld_clumping_2( globals.use_db_name, orig_snp );
    }

    if(orig.size()==0) {
        log("  no interval was generated\n");
        return false;
    }
    if(!log_intervals(globals.testfile +".int")) {
        log("  interval writing not completed\n");
        return false;
    }

    log( "\n" + int2str(orig.size())
         + " intervals saved in : " + globals.testfile + ".int\n\n");

    globals.testfile = globals.testfile + ".int";
    globals.convert_mode = false;

    return true;
}

bool InRichMain::make_lddb( const std::string & dbfile, const std::string & dbname )
{
    log("Creating LD database \n");

    if( fileExists(dbname) ) {
	log("\nDB file already exits : " + dbname + "\n  To create a new DB file with the same name, first delete the current one and rerun InRich.\n" +
            "  Or, rerun InRich using the current DB file\n\n");
	return false;
    }
    init_lddb( dbname, "" );


    if ( ! sql.is_open() ) return false;
    sql.query( " DROP INDEX IF EXISTS i1; ");
    std::ifstream F( dbfile.c_str() , std::ios::in );
    sql.begin();
    int done = 0;

    while ( ! F.eof() )
    {
        std::vector<std::string> lines = tokenizeLine( F );

	if(lines.size() == 0 ) break;
	if(lines.size() < 5 ) continue;

        std::string name = lines[0];
        int chr = get_chr(lines[1], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp_, bp1_, bp2_;
        if(!str2int(lines[2], bp_)) continue;
        if(!str2int(lines[3], bp1_)) continue;
        if(!str2int(lines[4], bp2_)) continue;
        uint64_t bp = bp_;
        uint64_t bp1 = bp1_;
        uint64_t bp2 = bp2_;

        // note: no format checking -- pls add (tab-delim, 4 cols, numeric)
        // also -- assumes numeric CHR coding -- pls edit to make X,Y,XY and M
        // map onto 23,24,25,26

        //F >> name >> chr >> bp >> bp1 >> bp2;

        if(name.compare("SNP")==0) {
            continue;
        }

        if( name == "." ) name = "";
        sql.bind_text( q_insert, ":varid" , name );
        sql.bind_int( q_insert, ":chr" , chr );
        sql.bind_int64( q_insert, ":bp" , bp );
        sql.bind_int64( q_insert, ":bp1" , bp1 );
        sql.bind_int64( q_insert, ":bp2" , bp2 );
        sql.step( q_insert );
        sql.reset( q_insert );

        ++done;

	if ( done % 10000 == 0 ) {
            std::cerr << "  " << done << " SNP information loaded\r";
            std::cout.flush();
	}
    }
    F.close();

    sql.finalise( q_insert );
    sql.commit();

    std::cerr << "\n  creatting db index ...";
    sql.query( " CREATE INDEX IF NOT EXISTS i1 ON lddat(varid); " );

    log("  " + int2str(done) + " SNP information loaded\n\n");
    return true;
}

std::string InRichMain::ld_clumping_2( const std::string & lddb_name,
                                       std::set<SNP> & orig_snp ) {

    // ------------------------------------------------------------
    // Note that SNPs are ordered by their chromosomal location
    // ------------------------------------------------------------
    std::vector<Interval> interval_list;
    interval_list.clear();

    std::set<SNP>::iterator itr = orig_snp.begin();
    int x_chr=-1;
    int x_snp_bp=-1;
    int skipped=0;
    int cnt=-1;
    int done = 0;
    int snp_no = snps.size();

    int x_right_start_pos = 0;
    int x_left_start_pos = 0;
    while(itr != orig_snp.end()) {
        done++;

        SNP cur_snp = *itr;

        // --------------------------------------------------
        // Check whether the SNP information exist in ld.db
        // --------------------------------------------------
        uint64_t tag_right, tag_left;
        bool found = false;

        //std::cout << "checking " << cur_snp.chr << " " << cur_snp.bp << "\n";

        // check with the SNP position
        // found = lddb.query( cur_snp.chr, cur_snp.bp, tag_left, tag_right );
        found = query( cur_snp.chr, cur_snp.bp, tag_left, tag_right );

        //std::cout << "checking " << cur_snp.chr << " " << cur_snp.bp <<  " " << found << "\n";

        if (!found) {
            log( "  ERROR: SNP " + cur_snp.name + " at chr"  + get_chr(cur_snp.chr)  +
                 ":"  + int2str(cur_snp.bp) + " not in LD database\n") ;
            skipped++;
            //x_chr = cur_snp.chr;
            //x_snp_bp = cur_snp.bp;

            itr++;
            continue;
        }
	
        // --------------------------------------------------
        // If not the first SNP && the previous interval resides within the left tagging area of current SNP
        // then extend the previous interval
        // --------------------------------------------------
        int ext_left_pos = tag_left;
        int ext_right_pos = tag_right;
        if(ext_left_pos < cur_snp.bp) {
            ext_left_pos = get_closest_ref_snp_pos_2( snp_no, cur_snp.chr, tag_left, true, x_left_start_pos );
	}
        if(ext_right_pos > cur_snp.bp) {
            ext_right_pos = get_closest_ref_snp_pos_2( snp_no, cur_snp.chr, tag_right, false, x_right_start_pos );
        }

        if ( itr->chr == x_chr  && ext_left_pos <= x_snp_bp ) {
            Interval & prev = interval_list.back();
            // when ext_right_pos > prior interval right point
            // need to update
            if(ext_right_pos > prev.bp2) {
                prev.bp2 = ext_right_pos;
	    	if(globals.verbose) {
                    log( "  extended interval " + prev.get_desc() + "\n") ;
		}
            }
        }
        // Otherwise, register it as a interval
        else {
	    // no right or left tagging SNP exists except itself
	    if( ext_left_pos == ext_right_pos ) {
		ext_right_pos++;
	    }
            interval_list.push_back( Interval( cur_snp.chr, ext_left_pos, ext_right_pos, 0, "int"+int2str(++cnt) ) );

	    if(globals.verbose) {
            	log( "  new interval int" + int2str(cnt) + ":" + int2str( ext_left_pos) + ".." + int2str( ext_right_pos ) 
                     + " generated from SNP " + cur_snp.name + " at " + get_chr(cur_snp.chr) + ":" + int2str(cur_snp.bp) + "\n") ;
	    }
        }

        x_snp_bp = ext_right_pos;
        x_chr = cur_snp.chr;

        itr++;

	if(cnt%10==0) {
            std::cerr << "  " << cnt << " interval generated\r";
            std::cout.flush();
	}
    }


    return register_interval( interval_list );
}

std::string InRichMain::ld_clumping( const std::string & lddb_name,
                                     std::set<SNP> & orig_snp ) {

    // ------------------------------------------------------------
    // Loading Hapmap-based LD Lookup Table
    // ------------------------------------------------------------

    // ------------------------------------------------------------
    // Note that SNPs are ordered by their chromosomal location
    // ------------------------------------------------------------
    std::vector<Interval> interval_list;
    interval_list.clear();

    std::set<SNP>::iterator itr = orig_snp.begin();
    int x_chr=-1;
    int x_snp_bp=-1;
    int skipped=0;
    int cnt=-1;
    int done = 0;
    int snp_no = snps.size();

    int x_right_start_pos = 0;
    int x_left_start_pos = 0;
    while(itr != orig_snp.end()) {
        done++;

        SNP cur_snp = *itr;

        // --------------------------------------------------
        // Check whether the SNP information exist in ld.db
        // --------------------------------------------------
        uint64_t _tag_right, _tag_left;
        bool found = false;

        //std::cout << "checking " << cur_snp.chr << " " << cur_snp.bp << "\n";

        // check with the SNP position
        // found = lddb.query( cur_snp.chr, cur_snp.bp, tag_left, tag_right );
        found = query( cur_snp.chr, cur_snp.bp, _tag_left, _tag_right );

        //std::cout << "checking " << cur_snp.chr << " " << cur_snp.bp <<  " " << found << "\n";

        if (!found) {
            log( "  ERROR: SNP " + cur_snp.name + " at chr"  + get_chr(cur_snp.chr)  +
                 ":"  + int2str(cur_snp.bp) + " not in LD database\n") ;
            skipped++;
            //x_chr = cur_snp.chr;
            //x_snp_bp = cur_snp.bp;

            itr++;
            continue;
        }
	
        int tag_left = _tag_left;
        int tag_right = _tag_right;

        // --------------------------------------------------
        // If not the first SNP && the previous interval resides within the left tagging area of current SNP
        // then extend the previous interval
        // --------------------------------------------------
        int ext_left_pos;
        int ext_right_pos;
        if ( itr->chr == x_chr && tag_left <= x_snp_bp ) {
	    // for current SNP
	    if(tag_right > cur_snp.bp) {
            	ext_right_pos = get_closest_ref_snp_pos( snp_no, cur_snp.chr, tag_right, false, x_right_start_pos );
	    }
	    else {
		ext_right_pos = tag_right;
	    }

            Interval & prev = interval_list.back();
            // when ext_right_pos > prior interval right point
            // need to update
            if(ext_right_pos > prev.bp2) {
                prev.bp2 = ext_right_pos;
	    	if(globals.verbose) {
                    log( "  extended interval " + prev.get_desc() + "\n") ;
		}
            }
        }
        // Otherwise, register it as a interval
        else {
	    if(tag_left < cur_snp.bp) {
            	ext_left_pos = get_closest_ref_snp_pos( snp_no, cur_snp.chr, tag_left, true, x_left_start_pos );
	    }
	    else {
	    	ext_left_pos = tag_left;
	    }
	    if(tag_right > cur_snp.bp) {
            	ext_right_pos = get_closest_ref_snp_pos( snp_no, cur_snp.chr, tag_right, false, x_right_start_pos );
	    }
	    else {
		ext_right_pos = tag_right;
	    }

	    // no right or left tagging SNP exists except itself
	    if( ext_left_pos == ext_right_pos ) {
		ext_right_pos++;
	    }
            interval_list.push_back( Interval( cur_snp.chr, ext_left_pos, ext_right_pos, 0, "int"+int2str(++cnt) ) );

	    if(globals.verbose) {
            	log( "  new interval int" + int2str(cnt) + ":" + int2str( ext_left_pos) + ".." + int2str( ext_right_pos ) 
                     + " generated from SNP " + cur_snp.name + " at " + get_chr(cur_snp.chr) + ":" + int2str(cur_snp.bp) + "\n") ;
	    }
        }

        x_snp_bp = ext_right_pos;
        x_chr = cur_snp.chr;

        itr++;

	if(cnt%10==0) {
            std::cerr << "  " << cnt << " interval generated\r";
            std::cout.flush();
	}
    }


    return register_interval( interval_list );
}
*/


bool InRichMain::log_interval_data() {
    //
    // List all interval/all-genes pairings, post de-duping, etc
    //
    if(interval_data.size()>0) {
    	interval_data.clear();
    }

    std::map< Interval , std::set<Interval> > overlap_all = get_overlap( test , all_genes_set );
    std::set<Interval>::iterator k = test.begin();
    while ( k != test.end() )
    {


        std::map< Interval , std::vector<int> >::iterator ip = globals.acceptable_positions.find( *k );
        int npos = -1;
        if ( ip != globals.acceptable_positions.end() )
            npos = ip->second.size();

        std::set<Interval> & o2 = overlap_all[ *k ];
        Interval interval = *k;
        if ( o2.size() == 0 )
        {
            interval_data.push_back( Interval2Gene( interval.get_desc() , k->n, npos, (int)o2.size() ));
        }
        else
        {
            std::set<Interval>::iterator jj = o2.begin();
            jj = o2.begin();
            while ( jj != o2.end() )
            {
                interval_data.push_back( Interval2Gene( interval.get_desc() , k->n, npos, (int)o2.size(), jj->name , all_genes[ jj->name ].desc ));
                ++jj;
            }
        }
        ++k;
    }


    return true;
}

bool InRichMain::write_output() {

    log("\nTotal of " + int2str(n_interval) + " intervals tested for " + int2str(n_set) + " target sets\n\n");

    cal_unplaced();

    log_target_path_data();
    log_interval_data();

    report_sig_study_hit( );
    log("\n\n\n" + sig_hit + "\n\n");

    report_clustering();

    log( "Writing output to " + globals.outroot + ".out.inrich ...\n\n");
    log_summary_file(globals.outroot + ".out.inrich", true);

    return true;
}


bool InRichMain::report_clustering() {
    if(globals.exam_bp_dist) {
        if( globals.exam_bp_dist_top > 0 && globals.exam_bp_dist_top < val_dist_p_no ) {
            report_dist_p( globals.exam_bp_dist_top, test_dist_p, test_dist_interval, test, globals.outroot );
        }
        else {
            report_dist_p( val_dist_p_no, test_dist_p, test_dist_interval, test, globals.outroot );
        }
        return true;
    }
    else {
        return false;
    }
}

bool InRichMain::run_second_permutation() {
    for ( int p = 0 ; p < globals.nboot ; p++ )
    {
        //
        // Create a matched set of regions
        //

        // ignore these -- already scored above
        std::map<Interval, int> dummy_match_tracker;

        std::set< Interval > matched = match( test ,
                                              snps ,
                                              dummy_match_tracker,
                                              globals.match_genes ? &all_genes_set : NULL  );

        // Goal is to obtain the minimum p-value over all sets for this
        // permutation (i.e. assuming this was the actual dataset). We can use
        // null_*_tracker[] to facilitate this.

        int min_test = globals.nrep + 1;


        //
        // Repeat test for each set
        //

        std::map<std::string, std::set<Interval> > ::iterator t = targets.begin();
        int pway = 0;



        // =====================================================
        // PL 2010.09.07 UPDATED
        // =====================================================
        std::vector<int> cur_study_hit(N_SIG_STUDY_HIT, 1);
        std::map< int, std::set<Interval> > cur_study_gene_hit;
        while ( t != targets.end() )
        {
            const int count = evaluate( matched , t->second );

            // How many times a larger or equal value observed? --> emp. p-val
            if ( null_tracker[ pway ][ count ] < min_test )
                min_test = null_tracker[ pway ][ count ];


            // =====================================================
            // PL 2010.09.07 UPDATED
            // =====================================================
            //double rep_emp_p = (double)(null_tracker[ pway ][ count ] +1) / (double)(globals.nrep+1);
            double rep_emp_p = (double)(null_tracker[ pway ][ count ]) / (double)(globals.nrep+1);

            for(unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
                if( rep_emp_p < SIG_STUDY_HIT[ ns ] ) {                    
                    get_overlapping_genes( matched, t->second, cur_study_gene_hit[ns]);
                    cur_study_hit[ ns ]++;
                }
            }

            ++t;
            ++pway;
        }

        // =====================================================
        // PL 2010.09.07 UPDATED
        // =====================================================
        for(unsigned  int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
            if( cur_study_hit[ ns ] >= (int)org_study_hit[ ns ]) {
                p_sig_study_hit[ ns ]++;
            }

            if( cur_study_gene_hit[ns].size() >= org_study_gene_hit[ ns ]) {
                p_sig_study_gene_hit[ ns ]++;
            }
        }

        //
        // We now have, for this replicate, what the min p-value would have
        // been; use this to score pcorr
        //

        for (int pway=0 ; pway < n_set; pway++) {
            if ( min_test <= pcount[pway] ) pcorr[pway]++;
        }

        std::cerr << p << " second-pass permutations         \r";
        std::cout.flush();
    }


    log(int2str(globals.nboot) + " second-pass permutations ( completed )\n");


    return true;
}


void InRichMain::log_interval_file(std::ofstream & OUT2) {
    OUT2 << "_O3\t" 
            << "INTERVAL" << "\t"
            << "N_SNP" << "\t"
            << "N_ALTLOC" << "\t"
            << "N_GENE" << "\t"
            << "GENE_ID" << "\t"
            << "GENE_DESC" << "\n";

    for(unsigned int i=0; i<interval_data.size(); i++) {
        OUT2 << "_O3\t" << interval_data.at(i).get_inrich_summary();
    }
}

void InRichMain::log_main_file(std::ofstream & OUT3){
    OUT3 << "_O2\t" 
            << "INTERVAL" << "\t"
            // << "N_UNSH" << "\t"
            << "GENE_LOC" << "\t"
            << "GENE_ID" << "\t"
            << "GENE_DESC" << "\t"
            << "TARGET" << "\t"
            << "P\n";

    for(unsigned int i=0; i<main_data.size(); i++) {
        OUT3 << "_O2\t" 
                << main_data.at(i).get_inrich_summary();
    }
}

void InRichMain::log_gene_file(std::ofstream & OUT4){
    OUT4 << "GENE_ID" << "\t"
            << "INTERVAL" << "\t"
            << "GENE_DESC" << "\t"
            << "P" << "\t"
            << "TARGET" << "\n";

    std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
    int   pway = 0;

    while ( t != targets.end() )
    {
        const std::set<Interval> & target = merge( t->second );

        //
        // Get all genes in this target set that were in the test interval
        //

        std::map< Interval , std::set<Interval> > overlap = get_overlap( test , target );

        std::map< Interval , std::set<Interval> >::iterator ii = overlap.begin();
        double pvalue = pcount[ pway ]/(double)( 1 + globals.nrep ); // pathway p-value

        while ( ii != overlap.end() ) // intervals overlapping with the target pathway
        {
            std::set<Interval> & o1 = ii->second;  // genes overlapping with each interval
            std::set<Interval>::iterator oo = o1.begin();
            while ( oo != o1.end() ) // for each gene
            {
                if ( pvalue < globals.pthresh )  // only target p-value is less than a threshold
                {
                    of_interest[ oo->name ] = ii->first ;   // for gene -> link an interval
                    std::map< std::string, double >::iterator k = of_interest_bestp.find( oo->name );
                    if ( k == of_interest_bestp.end() || pvalue < k->second )
                    {
                        of_interest_bestp[ oo->name ] = pvalue;
                        of_interest_bestset[ oo->name ] = t->first;  // register the target set with the best p-value
                    }
                }
                ++oo;
            }
            ++ii;
        }
        ++t;
        ++pway;

    }

    std::map<std::string,Interval>::iterator j = of_interest.begin();
    while ( j != of_interest.end() )
    {
        OUT4 << j->first << "\t"
                << j->second << "\t"
                << all_genes[ j->first ].desc << "\t"
                << of_interest_bestp[ j->first ] << "\t"
                << of_interest_bestset[ j->first ] << "\n";
        ++j;
    }
}


void InRichMain::log_summary_file(bool _append){
    //log("Writing the following output files:\n");

    if(globals.aligator) {
        //summary_filename = QDir::tempPath().toStdString() + "/" + globals.outroot + ".tmp.aligator";
        summary_filename =  "./" + globals.outroot + ".tmp.aligator";
    }
    else {
        //summary_filename = QDir::tempPath().toStdString() + "/" +  globals.outroot + ".tmp.inrich";
        summary_filename = "./" +  globals.outroot + ".tmp.inrich";
    }


    log_summary_file(summary_filename, _append);
}


void InRichMain::log_summary_file(std::string _filename, bool _append) {
    //log("Writing the following output files:\n");

    //summary_filename = CRandom::rand_letters(filename_size) + ".out.inrich";

    if(!_append) {
        remove(_filename.c_str());
    }

    std::ofstream OUT( (_filename).c_str() , std::ios::app );

    OUT << "----------------------------------------------------------------\n";
    OUT << "I. Main Analysis Results\n";
    OUT << "----------------------------------------------------------------\n";
    log_target_file(OUT);

    OUT << "\n\n\n";
    OUT << "----------------------------------------------------------------\n";
    OUT << "II. Interval-Gene-Target Summary\n";
    OUT << "----------------------------------------------------------------\n";
    log_main_file(OUT);

    /*
    OUT << "\n\n\n";
    OUT << "----------------------------------------------------------------\n";
    OUT << "III. Associated Gene Summary\n";
    OUT << "----------------------------------------------------------------\n";
    log_gene_file(OUT);
    */

    OUT << "\n\n\n";
    OUT << "----------------------------------------------------------------\n";
    OUT << "III. Interval Summary\n";
    OUT << "----------------------------------------------------------------\n";
    log_interval_file(OUT);
}


// =======================================================
// PL 2010.10.22 UPDATE
//
// To speed up, this check up is done later
// =======================================================

std::set<Interval>
        InRichMain::match( const std::set<Interval> & t ,
                           const std::vector<Interval> & pos ,
                           std::map<Interval,int> & match_tracker ,
                           std::set<Interval> * genes
                           )
{

    std::set<Interval> m;
    const int s = pos.size();


    std::set<Interval>::iterator i = t.begin();
    while ( i != t.end() )
    {
        // std::cout << "\t" << *i << "\n";
        match_interval(m, i, pos, s, match_tracker, genes);
        // std::cout << "\tmatch size=" << m.size() << "\n";
        ++i;
    }

    int no=0;
    const int mx = 10000;
    while ( 1 ) {
        std::set<Interval>::iterator itr;

        if( overlap_interval( m, itr )) {
            match_interval(m, itr, pos, s, match_tracker, genes);

            //std::cout << "overlapping match [" << no << "] deleted " << *itr << "\n";
            no++;

            m.erase(itr);
            //std::cout << "random match size=" << m.size() << "\n";
        }
        else {
            break;
        }

        if(no==mx) break;
    }

    // =====================================================================================================
    // PL 2011.10.02 UPDATE
    // to restrict the number of maximum intervals matching to each target gene to one as the original model
    // =====================================================================================================
    m = merge( m );

    std::map< Interval , std::set<Interval> > overlap2 = get_overlap( all_target_genes , m );
    std::map< Interval , std::set<Interval> > overlap1 = get_overlap( m , all_target_genes );

    m = gen_match_set(overlap1, overlap2, all_target_genes, m);
    m = remove_subinterval( m );
    m = merge( m );

    // =====================================================================================================
    // PL 2011.10.09 Compensate when the number of random intervals is smaller than the original 
    // =====================================================================================================
    bool added = false;
    int max_additional_no = 1;
    int added_no = 0;
    i = t.begin();
    int org_size = t.size();
    while ( m.size() < org_size && i != t.end() && added_no < max_additional_no ) {
        //std::cout << "\ttest org_size=" << org_size << " != match size=" << m.size() << "\n";
        match_interval(m, i, pos, s, match_tracker, genes);
	added = true;
        ++i;
	added_no++;
    }

    if(added) {
    	m = merge( m );
    }

    return m;
}


void
        InRichMain::match_interval( std::set<Interval> & m,
                                    std::set<Interval>::iterator & i,
                                    const std::vector<Interval> & pos ,
                                    const int s,
                                    std::map<Interval,int> & match_tracker ,
                                    std::set<Interval> * genes
                                    )
{
    // for each interval
    // need to match i->n SNPs
    //

    int distance = i->bp2 - i->bp1 + 1;
    double lower = distance - distance * globals.tol;
    double upper = distance + distance * globals.tol;
    if ( lower < 1 ) lower = 1;

    //
    // need to match # of genes?
    //

    int num_gene = 0;
    if ( globals.match_genes )
    {
        num_gene = num_in_interval( *i, *genes );
    }


    //
    // max # of trys to find match
    //

    const int mx = 1000000;
    int cnt = 0;

    bool density_okay = ( ! globals.use_map ) || i->n == 0 ;
    bool genenum_okay = ! globals.match_genes;

    while ( 1 )
    {

        int chr=-1, bp1=-1, bp2=-1;

        // get the implied distance, moving N SNPs forwards;
        // if this spans a chromosome, we will get negative distance

        // Use map positions to match on SNP density?

        if ( globals.use_map )
        {

            // select only from blessed, pre-computed positions?
            int r = 0;

            if ( globals.precompute
                 && globals.acceptable_positions.find( *i ) != globals.acceptable_positions.end() )
            {
                const std::vector<int> & ppos = globals.acceptable_positions[ *i ];
                if( ppos.size()>0 ) {
                    r = ppos[ CRandom::rand( (int)ppos.size() ) ] ;
                }
                else {
                    // select a SNP at random
                    r = CRandom::rand( s - i->n );
                }
            }
            else
            {
                // select a SNP at random
                r = CRandom::rand( s - i->n );
            }



            if ( globals.precompute )
            {
                density_okay = genenum_okay = true;

                if(i->n > 0) {
                    chr = pos[r].chr;
                    bp1 = pos[r].bp1;
                    bp2 = pos[ r + i->n - 1 ].bp1;
                }
                else {
                    // so that this region won't have an overlapping SNP
                    chr = pos[r].chr;
                    bp1 = pos[r].bp1 + 1;
                    bp2 = pos[ r + 1 ].bp1 - 1;
                }
            }
            else
            {
                int faked;

                if(i->n==0) {
                    faked = ( pos[ r + 1 ].bp1 - 1) - ( pos[ r ].bp1 + 1 ) + 1;
                }
                else {
                    faked = pos[ r + i->n - 1 ].bp1  - pos[r].bp1 + 1;
                }
                // the lenght of the randomly selected interval
                if ( faked >= lower && faked <= upper )
                {
                    density_okay = true;
                    if(i->n == 0) {
                        chr = pos[r].chr;
                        bp1 = pos[r].bp1 + 1;
                        bp2 = pos[ r + 1 ].bp1 - 1;
                    }
                    else {
                        chr = pos[r].chr;
                        bp1 = pos[r].bp1;
                        bp2 = pos[ r + i->n - 1 ].bp1;
                    }
                }
                else
                    density_okay = false;
            }
        }
        else
        {
            // Select interval from anywhere in the legal range
            // of fixed length in bp
            Interval i2 = globals.random();
            chr = i2.chr;
            bp1 = i2.bp1;
            bp2 = bp1 + (distance-1);


	    //std::cerr << "Here comes ! chr=" << chr << " bp1=" << bp1 << " bp2=" << bp2 << "\n";
        }


        //
        // Did we get a match?
        //

        if ( (!globals.precompute) && globals.match_genes && density_okay )
        {
            int g = num_in_interval( Interval( chr, bp1, bp2 ) , *genes );
            if ( g == num_gene ) genenum_okay = true;
        }


        //
        // Is this okay?
        //
        if ( density_okay && genenum_okay )
        {
            Interval i1( chr, bp1 , bp2, i->n );
            m.insert( i1 );

            // check not same by chance? okay to leave if is, but track this (i.e. indicative of small space)
            if ( i1 == *i && ! globals.exam_bp_dist )
            {
                match_tracker[*i]++;
            }

            break;
        }

        // If can't find anything, use original (i.e. will be conservative)
        if ( ++cnt > mx )
        {
            m.insert( *i );
            if(! globals.exam_bp_dist ) {
                match_tracker[*i]++;
            }
            break;
        }
    }
}


bool InRichMain::overlap_interval( const std::set<Interval> & s,
                                   std::set<Interval>::iterator & jtr )
{
    std::set<Interval>::iterator i = s.begin();
    std::set<Interval>::iterator last = i;
    ++i;
    while ( i != s.end() )
    {
        if ( *last == *i )
        {
            jtr = i;
            return true;
        }
        last = i;
        ++i;
    }
    return false;
}



bool InRichMain::overlap( const std::set<Interval> & s )
{
    std::set<Interval>::iterator i = s.begin();
    std::set<Interval>::iterator last = i;
    ++i;
    while ( i != s.end() )
    {
        if ( *last == *i )
        {
            return true;
        }
        last = i;
        ++i;
    }
    return false;
}



std::set<Interval> InRichMain::merge( const std::set<Interval> & intervals )
{
    std::set<Interval> r;
    std::set<Interval>::iterator i = intervals.begin();

    Interval current = *intervals.begin();

    int no=0;
    while ( i != intervals.end() )
    {
        //std::cout << no << " " << *i << "\n";

        // if overlap then merge
        if ( *i == current )
        {
            current.combine( *i );
        }
        else
        {
            r.insert( current );
            current = *i;
        }
        ++i;
        no++;
    }
    //std::cout << no << " intervals\n";

    // set, so will not dupe last element
    r.insert(current);
    return r;
}

std::set<Interval>  InRichMain::get_all_targets() {
    std::set<Interval> all_targets;
    std::map< std::string, std::set<Interval> >::iterator t = targets.begin();
    while ( t != targets.end() )
    {
        // =====================================================
        // PL 2010.08.31 ADDED --- Why commented?
        // =====================================================
        t->second = merge( t->second );

        std::set<Interval>::iterator i = t->second.begin();
        while ( i != t->second.end() )
        {
            all_targets.insert( *i );
            ++i;
        }
        ++t;
    }

    return all_targets;
}


// ,====================================================
// PL 2010.07.18 ADDED
// ====================================================
std::set<Interval> InRichMain::remove_subinterval( const std::set<Interval> & _org_test) {

    std::set<Interval> new_test;
    std::set<Interval>::iterator i = _org_test.begin();

    while ( i != _org_test.end() ) {
        bool exclude=false;

        std::set<Interval>::iterator j = _org_test.begin();

        while ( j != _org_test.end() && j->chr <= i->chr ) {

            // not self but *i is included in *j
            if( j->chr == i->chr && j->name.compare(i->name)!=0 && j->bp1<=i->bp1 && j->bp2>=i->bp2) {
                exclude=true;
                break;
            }

            ++j;
        }

        if(!exclude) {
            new_test.insert(*i);
        }

        ++i;
    }

    return new_test;
}

std::vector< std::vector<int> > InRichMain::init_study_hit( int nrep )
{
    std::vector< std::vector<int> > sig_study_hit;

    try{
        // ---------------------------------------------------
        // init the vector for significant study hits
        // ---------------------------------------------------
        sig_study_hit.resize(nrep);
        for(int i=0; i<nrep; i++) {
            for(unsigned int j=0; j<N_SIG_STUDY_HIT; j++) {
                sig_study_hit[i][j] = 0;
            }
        }
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRich::init_study_hit\n";
        exit(1);
    }

    return sig_study_hit;
}

void InRichMain::count_sig_study_hit_all ( int rep,
                                           int nrep,
                                           std::vector< std::vector<int> > & sig_study_hit,
                                           const std::vector<int> & pcount )
{
    nrep++;

    for(unsigned int i=0; i<pcount.size(); i++) {
        double pc = (double)pcount[i]/(double)nrep;

        for(unsigned int j=0; j<N_SIG_STUDY_HIT; j++) {
            if( pc <= SIG_STUDY_HIT[j] ) {
                sig_study_hit[rep][j]++;
            }
        }
    }
}


void InRichMain::count_sig_study_hit ( int rep,
                                       int nrep,
                                       std::vector< std::vector<int> > & sig_study_hit,
                                       const int pcount )
{
    nrep++;

    double pc = (double)pcount/(double)nrep;
    for(unsigned int j=0; j<N_SIG_STUDY_HIT; j++) {
        if( pc <= SIG_STUDY_HIT[j] ) {
            sig_study_hit[rep][j]++;
        }
    }
}

void InRichMain::report_sig_study_hit ()
{
    std::ostringstream myLog;

    /*
    myLog << "----------------------------------------------------------------\n"
            << "Target_P_Threshold " << "\t"
            << "Number_of_Targets" << "\t"
            << "Significance" << "\n"
            << "----------------------------------------------------------------\n";

    for(int i=0; i<N_SIG_STUDY_HIT ; i++) {
        myLog   << SIG_STUDY_HIT[ i ] << "\t"
                << org_study_hit[ i ] << "\t"
                << (double)(p_sig_study_hit[ i ]+1)/(double)(globals.nboot+1) << "\n" ;

    }

    myLog << "\n\n"; */

    myLog <<  "----------------------------------------------------------------\n"
            << std::setw( 18 ) <<  "Target_P_Threshold" 
            << std::setw( 28) <<  "Uniq_Gene_No_in_Targets" 
            << std::setw( 18) <<  "Significance" << "\n"
            << "----------------------------------------------------------------\n";

    for(unsigned int i=0; i<N_SIG_STUDY_HIT ; i++) {
        myLog   << std::setw( 18) <<  SIG_STUDY_HIT[ i ] 
                << std::setw( 28) <<  org_study_gene_hit[ i ] 
                << std::setw( 18) <<  (double)(p_sig_study_gene_hit[ i ]+1)/(double)(globals.nboot+1) << "\n" ;

    }

    myLog << "\n\n";

    sig_hit = myLog.str();
}

// ====================================================
// PL 2010.08.30 ADDED
// ====================================================
std::set<PairIntervalDist> InRichMain::compute_dist_interval( std::set<Interval> & regions ) {
    if( CAL_DIST_IF_ADJACENT ) {
        return compute_dist_interval_adjacent( regions );
    }
    else {
        return compute_dist_interval_all( regions );
    }
}

std::set<PairIntervalDist> InRichMain::compute_dist_interval_all ( std::set<Interval> & regions ) {
    std::set<PairIntervalDist> dist_int;

    std::set<Interval>::iterator i = regions.begin();

    while ( i != regions.end() ) {
        //debugging
        // std::cout << "interval i=" << *i << "\n";

        std::set<Interval>::iterator j = regions.begin();
        while ( j != regions.end()) {
            //debugging
            // std::cout << "interval j=" << *j << "\n";

            if( j->name.compare(i->name) == 0 ) {
                ++j;
                break;
            }
            ++j;
        }

        // examine the pair distances
        while ( j != regions.end() ) {
            if ( CAL_DIST_IF_SAME_CHR ) {
                if(j->chr == i->chr) {
                    dist_int.insert( PairIntervalDist( *i, *j ) );
                    /* if( j->name.compare(i->name) == 0 ) {
                                        std::cout << "New " << " PairIntervalDist( " << *i << " and " << *j << "\n";
                                        std::cout << j->name << "," << i->name << "\n";
                                } */
                }
            }
            else {
                dist_int.insert( PairIntervalDist( *i, *j ) );
                /* if( j->name.compare(i->name) == 0 ) {
                                        std::cout << "New Here " << " PairIntervalDist( " << *i << " and " << *j << "\n";
                                        std::cout << j->name << "," << i->name << "\n";
                        } */
            }
            ++j;
        }

        ++i;
    }

    return dist_int;
}

std::set<PairIntervalDist> InRichMain::compute_dist_interval_adjacent ( std::set<Interval> & regions ) {
    std::set<PairIntervalDist> dist_int;

    std::set<Interval>::iterator i = regions.begin();

    while ( i != regions.end() ) {
        Interval cur_i = *i;

        ++i;
        if( i != regions.end() ) {
            if ( CAL_DIST_IF_SAME_CHR ) {
                if(cur_i.chr == i->chr) {
                    dist_int.insert( PairIntervalDist( cur_i, *i ) );
                }
            }
            else {
                dist_int.insert( PairIntervalDist( cur_i, *i ) );
            }
        }
    }

    return dist_int;
}

std::vector<int> InRichMain::init_test_p ( int nsize ) {
    std::vector<int> dist_p;

    try {
        dist_p.resize(nsize);
        for(int i=0; i<nsize; i++) {
            dist_p[i] = 0;
        }
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRich::init_test_p(" + int2str(nsize) + ")\n";
        exit(1);
    }


    return dist_p;
}

void InRichMain::report_dist_p( int val_p_no,
                                std::vector<int> & dist_p,
                                std::set<PairIntervalDist> & test_dist,
                                std::set<Interval> & test,
                                std::string filename ) {


    log( "\n----------------------------------------------------------------\n" );
    log( " Interval Clustering Analysis Results\n" );
    log( "----------------------------------------------------------------\n" );
    log( "RANK\tINT1\tINT2\tDIST\tP\tINT1_POS\tINT2_POS\n" );

    std::set<PairIntervalDist>::iterator itr = test_dist.begin();
    for(int i=0; i<val_p_no; i++) {
        if( itr == test_dist.end() ) {
            // this cannot happen!
            std::cerr << "ERROR: \n" ;
            break;
        }

        double p = (double)(dist_p[i]+1)/(globals.nrep+1);

        char info1[50] = "N/A";
        char info2[50] = "N/A";
        std::set<Interval>::iterator jtr = test.begin();
        int all_find = 0;
        while( jtr != test.end() ) {
	    std::string chr_code = get_chr(jtr->chr);
            if(jtr->name.compare(itr->name1)==0) {
                sprintf(info1, "chr%s:%d..%d", chr_code.c_str(), jtr->bp1 , jtr->bp2);
                all_find++;
            }
            if(jtr->name.compare(itr->name2)==0) {
                sprintf(info2, "chr%s:%d..%d", chr_code.c_str(), jtr->bp1 , jtr->bp2);
                all_find++;
            }
            if(all_find==2) {
                break;
            }
            ++jtr;
        }

        if(p < 0.05) {
            log( int2str(i+1) + "\t" + itr->name1 + "\t" + itr->name2 + "\t" + int2str(itr->dist) + "*\t" + double2str(p) + "\t" + info1 + "\t" + info2 + "\n");
        }
        else {
            log( int2str(i+1) + "\t" + itr->name1 + "\t" + itr->name2 + "\t" + int2str(itr->dist) + "\t" + double2str(p) + "\t" + info1 + "\t" + info2 + "\n");
        }

        ++itr;
    }     

}

int InRichMain::update_dist_p( int val_p_no,
                               std::vector<int> & dist_p,
                               std::set<PairIntervalDist> & org_dist,
                               std::set<PairIntervalDist> & rand_dist ) {

    int np = dist_p.size();

    std::set<PairIntervalDist>::iterator org = org_dist.begin();
    std::set<PairIntervalDist>::iterator ran = rand_dist.begin();
    int i=0;
    while(i<np) {
        if( org != org_dist.end() && ran != rand_dist.end() ) {
            if( org->dist >= ran->dist ) {
                dist_p[i]++;
            }

            // std::cout << "i=" << i << " org->dist=" << org->dist << " ran->dist=" << ran->dist << "\n";
        }
        else {
            break;
        }

        ++org;
        ++ran;
        ++i;
    }

    if( val_p_no < i ) {
        return val_p_no;
    }
    else {
        return i;
    }
}




int InRichMain::evaluate( const std::set<Interval> & test ,
                          const std::set<Interval> & target ) {


    std::map< Interval , std::set<Interval> > overlap = get_overlap( test , target );
    //std::cerr << "test.size=" << test.size() << " overlap=" << overlap.size() << " target.size=" << target.size() <<  "\n";

    // number of test regions with 1+ target
    if ( globals.test == INTERVAL )
        return overlap.size() >= globals.min_cnt ? overlap.size() : 0 ;


    // Alternatively, get number of targets with 1+ test

    std::set<Interval> otarget;
    std::map< Interval , std::set<Interval> >::iterator i = overlap.begin();
    while ( i != overlap.end() )
    {
        std::set<Interval>::iterator j = i->second.begin();
        while ( j != i->second.end() )
        {
            otarget.insert( *j );
            ++j;
        }
        ++i;
    }
    return otarget.size() >= globals.min_cnt ? otarget.size() : 0 ;
}


void InRichMain::get_overlapping_genes( const std::set<Interval> & test ,
                                        const std::set<Interval> & target ,
                                        std::set<Interval> & otarget ) {

    std::map< Interval , std::set<Interval> > overlap = get_overlap( test , target );
    std::map< Interval , std::set<Interval> >::iterator i = overlap.begin();
    while ( i != overlap.end() )
    {
        std::set<Interval>::iterator j = i->second.begin();
        while ( j != i->second.end() )
        {
            otarget.insert( *j );
            ++j;
        }
        ++i;
    }
    return;
}




void InRichMain::init_second_permutation() {
    cal_study_sig_hit();

    //
    // We now have pcount[pathway] scored with 1..(R+1) hits
    // and we've also constructed the null distribution of scores
    // for each pathway (null_tracker[]).  Now we can loop through
    // a second time, repeating the exact same series of permutations,
    // and score the multiple test correction
    //

    //
    // Make cumulative distributions of null_*_tracker
    //

    for (int p=0;p<n_set;p++)
    {
        int t = 0;
        std::map<int,int>::reverse_iterator i = null_tracker[p].rbegin();
        while ( i != null_tracker[p].rend() )
        {
            t += i->second;
            i->second = t;
            ++i;
        }
    }

    try {
        pcorr.resize( n_set, 1 );
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRich::init_second_permutation(" + int2str(n_set) + ")\n";
        exit(1);
    }
}

void InRichMain::init_first_permutation() {
    // We also need to make a table for quick subsequent calculation of the
    // "nested permutations", to implement the multiple test correction For
    // each pathway, the distribution of scores (for t(interval) only, for
    // now)

    int pway = 0;

    try {
        null_tracker.resize( n_set );
        count.resize( n_set, 0 );
        actual_count.resize( n_set, 0 );
        pcount.resize(n_set, 1);
        p_stat_sum.resize(n_set, 0);

        //
        // Score originals, and save, for all sets
        //


        std::map<std::string, std::set<Interval> >::iterator t = targets.begin();

        while ( t != targets.end() )
        {
            count[ pway ] = evaluate( test , t->second );
            null_tracker[ pway ][ count[ pway] ]++;

            ++t;
            ++pway;
        }

        //
        // Because we set the count to 0 if below min_cnt, we want to get the true
        // count, if different from above, for output purposes
        //

        if ( globals.min_cnt == 1 ) actual_count = count;
        else
        {
            int tmp = globals.min_cnt;
            globals.min_cnt = 1;

            t = targets.begin();
            int pway = 0;
            while ( t != targets.end() )
            {
                actual_count[ pway ] = evaluate( test , t->second );
                ++t;
                ++pway;
            }

            globals.min_cnt = tmp;
        }
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRich::init_first_permutation(" + int2str(n_set) + ")\n";
        exit(1);
    }
    catch(std::exception a) {
        std::cerr << "Exception: InRich::init_first_permutation(" << a.what() << ")\n";
	std::cerr << "pway=" << pway << "\n\n";
        exit(1);
    }

}

std::set<Interval> InRichMain::make_compact( int cnt, std::set<Interval> & orig, std::set<Interval> & genes) {

    std::set<Interval>::iterator itr = orig.begin();

    std::set<Interval> compact_orig;

    int drop_no = 0;
    while( itr!= orig.end()) {

        std::set<Interval>::iterator jtr = genes.begin();
        bool on_gene=false;

        while (jtr != genes.end()) {
            // if input interval overlaps with a gene
            if(*itr == *jtr) {
                compact_orig.insert( *itr );
                on_gene = true;
                break;
            }
            jtr++;

        }

        if(! on_gene) {
            drop_no++;
        }

        itr++;
    }

    return compact_orig;
}

bool InRichMain::log_intervals(std::string _file) {
    std::ofstream OUT ( _file.c_str() , std::ios::out );
    std::set<Interval>::iterator itr = orig.begin();
    while( itr!= orig.end()) {
        OUT << itr->chr << " " << itr->bp1 << " " << itr->bp2 << "\n";
        //OUT << *itr << "\n";
        itr++;
    }
    OUT.close();
    return true;
}

std::string InRichMain::register_interval ( std::vector<Interval> & int_list ) {
    orig.clear();

    double avg_size = 0;
    double max_size = 0;
    double min_size = 10000000;
    int interval_no = int_list.size();

    std::string warning_msg = ""; //Your SNP MAP file doesn't match with your assoc SNP file!\n\n";
    for(int i=0; i< interval_no; i++) {
	if( int_list[i].bp2<=0 || int_list[i].bp1 <= 0 ) {
	    continue;
	}

        int size = int_list[i].bp2 - int_list[i].bp1;

        if( size<=0 ) {
            std::cerr << "\nERROR: Null Interval Detected " << int_list[i] ;
            warning_msg += "\nERROR: Null Interval Detected " + int_list[i].get_desc();
	    continue;
        }

        avg_size += size;

        if(max_size < size) {
            max_size = size;
        }
        if(min_size > size) {
            min_size = size;
        }

        orig.insert ( int_list[i] );
    }

    avg_size /= interval_no*1000;
    max_size /= 1000;
    min_size /= 1000;

    log( "  " + int2str(orig.size()) + " interval generated\n\n" +
         "  average interval size : " + double2str(avg_size) + " kb\n" +
         "  maximum interval size : " + double2str(max_size) + " kb\n"
         "  minimum interval size : " + double2str(min_size) + " kb\n" +
         warning_msg);

    std::string dummy = "";
    return dummy;
}



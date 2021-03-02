#include "inrichmultimain.h"
#include "loaders.h"
#include "crandom.h"

#include <stdlib.h>

InRichMultiMain::InRichMultiMain()
{
}

InRichMultiMain::~InRichMultiMain()
{
}



InRichMultiMain::InRichMultiMain(InputData & _globals) : InRichMain(_globals)
{

}


void InRichMultiMain::inrich_multi_test() {

    init_data();

    // ------------------------------------------------------
    // 1. Load Files : Use the original version
    // ------------------------------------------------------
    if(!load_files()) {
        std::cout << "\n --- halting --- no input data to test : interval/genes/targets ---\n";
        return;
    }

    // ------------------------------------------------------
    // 2. Merge Data
    // ------------------------------------------------------
    if(!merge_data()) {
        log("ERROR: Data not merged\n");
        return;
    }


    // ------------------------------------------------------
    // 3. Optional Checking
    // ------------------------------------------------------
    if(!init_optional_load()) {
        log("ERROR: Optional Data not Loaded\n");
        return;
    }

    // ------------------------------------------------------
    // 4. Permutation
    // ------------------------------------------------------
    log("\n");
    init_first_permutation();
    if(!run_first_permutation()) {
        log("ERROR: First permutation not completed\n");
        return;
    }


    init_second_permutation();
    if(!run_second_permutation()) {
        log("ERROR: Second permutation not completed\n");
        return;
    }

    // ------------------------------------------------------
    // 5. File Output
    // ------------------------------------------------------
    if(!write_output()) {
        log("ERROR: File writing not completed\n");
        return;
    }
}


bool InRichMultiMain::load_files() {
    int file_read = 0;
    int file_tried = 0;

    if(globals.convert_mode) {
        file_tried++;
        if(load_file(ASSOC_SNP)){
            file_read++;
            log("\n");
        }
    }
    else {
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
    }

    if(globals.use_map) {
    	file_tried++;
    	if(load_file(SNP_MAP)){
        	file_read++;
    	}
    }

    if(!globals.convert_mode) {
        file_tried++;
        if(load_file(GENE)){
            file_read++;
        }
        file_tried++;
        if(load_file(TARGET_SET)){
            file_read++;
        }
        file_tried++;
        if(load_file(ASSOC)){
            file_read++;
        }
    }

    if( file_read == file_tried ) {
        return is_data_ok();
    }
    else {
        return false;
    }
}

bool InRichMultiMain::is_data_ok() {

    std::map<std::string, std::set<Interval> >::iterator itr = multi_orig.begin();
    while(itr!=multi_orig.end()) {
        if(itr->second.size()==0) {
            return false;
        }
        itr++;
    }

    if ( globals.use_map ) {
    	std::map<std::string, std::vector<Interval> >::iterator jtr = multi_snps.begin();
    	while(jtr!=multi_snps.end()) {
        	if(jtr->second.size()==0) {
            		return false;
        	}
        	jtr++;
    	}
    }

    return InRichMain::is_data_ok();
}

bool InRichMultiMain::load_file(int _mode) {
    my_log = "";

    int no = -1;
    std::string res_msg;

    switch(_mode) {

    case BACKGROUND:
        no = load_background_data();   // ok
        break;

    case RANGE:
        no = load_range_data();  // ok
	log( double2str((double)globals.total_seq/(double)1000000) +  "Mb total sequence length on " + int2str(globals.scaffold.size()) + " chromosomes\n");
        break;

    case SNP_MAP:
        if(globals.multi_mapfile) {
            //log( "  reading multi-map file positions\n");
            multi_snps.clear();
            no = load_multi_snps( globals.mapfile );
            //log( "  read " + int2str(no) + " map positions in total\n");
        }
        else {            
            snps.clear();
            snps = load_snps( globals.mapfile, globals, err_msg );   // ok
            no = snps.size();

            log( double2str(globals.total_seq/(double)1e6) + "Mb total sequence length on "
                 + int2str(globals.scaffold.size()) + " chromosomes\n");
            log( "  read " + int2str(no) + " map positions\n");
        }
        break;

    case GENE:
        all_genes.clear();   // ok
        all_genes = loader_genes( globals.genefile, background_genes, globals , err_msg );
        no = all_genes.size();
        all_genes_set = get_set( all_genes );
        log( "  read " + int2str(no) + " reference genes\n");
        break;

    case TARGET_SET: // ok
        targets.clear();
        targets = loader_targets( globals.targetfile, all_genes, globals, err_msg, res_msg );
        no = targets.size();
	log (res_msg);
        break;

    case ASSOC:
        no = load_multi_assoc();   // need to be defined
        break;

    case ASSOC_SNP:
        no = load_multi_assoc_snps();   // need to be defined
        break;


    default:
        break;
    }

    if(err_msg.size()>0) {
	log(err_msg);
    }

    if(no>0) {
        return true;
    }
    else {
        return false;
    }
}


int InRichMultiMain::load_multi_snps( std::string _list_file ) {

    if ( ! fileExists( _list_file ) )
    {
        std::cerr << "could not find map list file " << _list_file << "\n";
        exit(1);
    }

    std::ifstream TEST( _list_file.c_str() , std::ios::in );
    int cnt = 0;
    int err_line=0;
    while ( !TEST.eof() ) {
    	std::vector<std::string> lines = tokenizeLine( TEST );

	if(lines.size() == 0) break;
		
	if(lines.size() != 2) {
		err_line++;
	 	continue;
	}

        std::string set = lines[0];
        stoupper(set);

	std::string mapfile = lines[1];

        if ( ! fileExists( mapfile ) )
        {
            std::cerr << "could not find snp map file " <<  mapfile << "\n";
            exit(1);
        }

        int no = load_snp_map( set, mapfile );
        cnt += no;

        log( "  read " + int2str(no) + " map positions (" + set + ")\n");
    }
    TEST.close();


    return cnt;
}

int InRichMultiMain::load_snp_map( std::string _set, std::string _mapfile )
{
    err_msg = "";
    multi_snps[_set] = load_snps(_mapfile, globals, err_msg);
    multi_scaffold[_set] = globals.scaffold;

    if(err_msg.size() > 0) {
    	log( err_msg );
    }
    log( double2str(globals.total_seq/(double)1e6) + "Mb total sequence length on "
         + int2str(globals.scaffold.size()) + " chromosomes\n");

    globals.scaffold.clear();

    return multi_snps[_set].size();
}

int InRichMultiMain::load_multi_assoc( )
{
    multi_orig.clear();

    // test regions

    if ( ! fileExists( globals.testfile ) )
    {
        std::cerr << "could not find test interval list " << globals.testfile << "\n";
        exit(1);
    }

    std::ifstream TEST( globals.testfile.c_str() , std::ios::in );
    int cnt = 0;
    int err_line = 0;
    while ( ! TEST.eof() )
    {
        std::vector<std::string> lines = tokenizeLine( TEST );

        if ( lines.size() == 0 ) break;
        if ( lines.size() != 4 ) { 
		err_line++;
		continue;
	}

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp1 , bp2;
        if ( ! str2int( lines[1] , bp1 ) ) continue;
        if ( ! str2int( lines[2] , bp2 ) ) continue;

        stoupper(lines[3]);
        multi_orig[lines[3]].insert( Interval( chr , bp1, bp2 , 0 , "int"+int2str(cnt++) ) );

        if ( cnt == globals.topn ) break;
    }

    TEST.close();

    if( err_line > 0 ) {
	log ( "  warning: 4 column (chr bp1 bp2 set_name) required in " + globals.testfile + "\n");
    }

    log( "\n  read " + int2str(cnt) + " intervals in total\n");
    n_multi_set_no = multi_orig.size();
    std::map<std::string, std::set<Interval> >::iterator itr = multi_orig.begin();
    while(itr != multi_orig.end()) {
        log ( "  read " + int2str(itr->second.size()) + " intervals (" + itr->first  + ")\n" );
        itr++;
    }
    log( "\n");

    return cnt;
}

int InRichMultiMain::load_multi_assoc_snps( )
{
    multi_orig_snps.clear();

    // data format
    // CHR BP SNP_ID P

    if ( ! fileExists( globals.testfile ) )
    {
        std::cerr << "could not find test SNP list " << globals.testfile << "\n";
        exit(1);
    }

    std::ifstream TEST( globals.testfile.c_str() , std::ios::in );
    int cnt = 0;
    int err_line = 0;
    while ( ! TEST.eof() )
    {

        std::vector<std::string> lines = tokenizeLine( TEST );

        int size = lines.size();
        if ( size == 0 ) break;
        if ( size != 5  ) { 
		err_line++;
		continue;
	}

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp;
        if ( ! str2int( lines[1] , bp ) ) continue;

        std::string name = lines[2];

        double p;
        if ( ! str2double( lines[3] , p ) ) continue;

        stoupper(lines[4]);
        multi_orig_snps[lines[4]].insert( SNP( chr , bp, name , p) );

        cnt++;
    }

    TEST.close();

    if( err_line > 0 ) {
	log ( "  warning: 5 column (chr bp snp_id p set_name) required in " + globals.testfile + "\n");
    }
    log( "\n  read " + int2str(cnt) + " SNPs in total\n");
    n_multi_set_no = multi_orig_snps.size();
    std::map<std::string, std::set<SNP> >::iterator itr = multi_orig_snps.begin();
    while(itr != multi_orig_snps.end()) {
        log ( "  set " + itr->first + " : " + int2str(itr->second.size()) + " SNPs\n" );
        itr++;
    }

    return cnt;
}


bool InRichMultiMain::merge_data() {
    Log mylog( globals.outroot + ".out.inrich" );

    std::map<std::string, std::set<Interval> >::iterator itr = multi_orig.begin();

    if( globals.match_genes ) {
        while(itr != multi_orig.end()) {
            int drop_no = itr->second.size();
            itr->second = make_compact( drop_no, itr->second, all_genes_set );
            drop_no = drop_no - itr->second.size();
            mylog << "  " << itr->second.size()
                    << " test elements on gene regions (" + itr->first + ")\n"
                    << "  "  << drop_no << " non-genic test elements dropped\n";
            itr++;
        }
    }
    all_genes_set = merge( all_genes_set );

    update_all_genes_map();

    std::set<Interval> all_target_genes = get_all_targets();



    //
    // Merge any test regions that span the same target (from *any included
    // set*), or with self, and place in test[].  This way, we are only
    // dealing with a single set of test intervals downstream (and so can
    // permute only once per replicate).
    //

    n_interval = 0;
    itr = multi_orig.begin();
    while(itr != multi_orig.end()) {
        std::map< Interval , std::set<Interval> > overlap2 = get_overlap( all_target_genes , itr->second );

        std::map< Interval , std::set<Interval> > overlap1 = get_overlap( itr->second , all_target_genes );

        orig = itr->second;
        test.clear();
        gen_test_set(overlap1, overlap2);

        test = remove_subinterval( test );

        multi_test[itr->first].clear();
        multi_test[itr->first] = merge( test );

        n_interval += multi_test[itr->first].size();

        itr++;
    }

    n_set = targets.size();

    mylog << "\nAfter merging, " << all_genes_set.size() <<  " non-overlapping reference genes\n";
    mylog << "After merging, " <<  all_target_genes.size() <<  " non-overlapping genes in target sets\n";

    itr = multi_test.begin();
    while(itr != multi_test.end()) {
        mylog << "After merging, " + int2str(itr->second.size()) + " non-overlapping intervals (" + itr->first + ")\n";

        itr++;
    }
    mylog << "\n";

    return true;
}


bool InRichMultiMain::init_optional_load() {
    std::map<std::string, std::set<Interval> >::iterator itr = multi_test.begin();

    //
    // If using a MAP, assign SNP counts, and check that each falls on a MAP location
    //

    if ( globals.use_map ) {
        log("Assigning snp counts ");
        while(itr != multi_test.end()) {
            if( globals.multi_mapfile ) {
                if(!assign_snp_count( itr->second, multi_snps[itr->first] )) {
                    return false;
                }
            }
            else {
                if(!assign_snp_count( itr->second, snps )) {
                    return false;
                }
            }            
            itr++;
        }
        log(" : DONE\n\n");
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

        itr = multi_test.begin();
        while(itr != multi_test.end()) {
            //log("  set " + itr->first + "\n");
            if(globals.multi_mapfile) {
                if(!pre_compute( itr->second, multi_snps[itr->first] , globals.match_genes ? &all_genes_set : NULL )) {
                    return false;
                }
                //log("  (" + itr->first + "\n");
            }
            else {
                if(!pre_compute( itr->second, snps , globals.match_genes ? &all_genes_set : NULL )) {
                    return false;
                }
            }            
            itr++;
        }
    }

    return true;
}


void InRichMultiMain::init_first_permutation() {
    // We also need to make a table for quick subsequent calculation of the
    // "nested permutations", to implement the multiple test correction For
    // each pathway, the distribution of scores (for t(interval) only, for
    // now)

    try {
        std::map<std::string, std::set<Interval> >::iterator ttr = multi_test.begin();

        while(ttr != multi_test.end()) {
            multi_null_tracker[ttr->first].clear();
            multi_count[ttr->first].clear();
            multi_count[ttr->first].resize( n_set, 0 );

            ttr++;            
        }

        // actual_count.resize( n_set, 0 );
        pcount.resize(n_set, 1);

        //
        // Score originals, and save, for all sets
        //
        ttr = multi_test.begin();
        while(ttr != multi_test.end()) {
            std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
            int pway = 0;

            while ( t != targets.end() )
            {
                multi_count[ttr->first][ pway ] = evaluate( ttr->second , t->second );
                multi_null_tracker[ttr->first][ pway ].push_back( multi_count[ttr->first][ pway ] );

                ++t;
                ++pway;
            }

            ttr++;
        }


    }
    catch(std::bad_alloc a) {
        exit(1);
    }
}



int InRichMultiMain::evaluate( const std::set<Interval> & test ,
                               const std::set<Interval> & target ) {

    std::map< Interval , std::set<Interval> > overlap = get_overlap( test , target );


    // number of test regions with 1+ target
    if ( globals.test == INTERVAL ) {
        return overlap.size();
        //return overlap.size() >= globals.min_cnt ? overlap.size() : 0 ;
    }

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

    return otarget.size();
    //return otarget.size() >= globals.min_cnt ? otarget.size() : 0 ;
}


bool InRichMultiMain::run_first_permutation () {

    CRandom::srand( globals.seed ) ;

    for ( int p = 0 ; p < globals.nrep ; p++ )
    {
        std::map< std::string, std::set< Interval > > multi_matched;
        std::map<std::string, std::set<Interval> >::iterator ttr = multi_test.begin();

        // -------------------------------------------------------
        // select matched random set for each original assoc list
        // -------------------------------------------------------
        while(ttr != multi_test.end()) {
	    if(globals.use_map) {
            	globals.scaffold = multi_scaffold[ttr->first];
	    }
	    
            if(globals.multi_mapfile) {
                multi_matched[ttr->first] = match( ttr->second,
                                                   multi_snps[ttr->first],
                                                   match_tracker,
                                                   globals.match_genes ? &all_genes_set : NULL  );
            }
            else {
                multi_matched[ttr->first] = match( ttr->second,
                                                   snps,
                                                   match_tracker,
                                                   globals.match_genes ? &all_genes_set : NULL  );
            }

            ttr++;
        }

        // -------------------------------------------------------
        // repeat test for each target
        // -------------------------------------------------------
        std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
        int pway = 0;
        while ( t != targets.end() )
        {
            int perm_statistic;
            bool increase = false;

            //std::cout << t->first << " " << t->second.size() << "\n";
            std::map< std::string, std::set< Interval > >::iterator itr = multi_matched.begin();
            while(itr != multi_matched.end()) {
                perm_statistic = evaluate( itr->second, t->second );
                multi_null_tracker[itr->first][ pway ].push_back( perm_statistic );

                if( multi_count[itr->first][ pway ] < (int)globals.min_cnt || perm_statistic >= multi_count[itr->first][ pway ] ) {
                    increase = true;
                }
                itr++;

            }
            if ( increase ) ++pcount[ pway ];

            ++t;
            ++pway;
        }
        std::cout << p << " first-pass permutations         \r";
        std::cout.flush();
    }

    log( int2str(globals.nrep) + " first-pass permutations ( completed )\n");

    return true;
}


void InRichMultiMain::cal_study_sig_hit() {
    try {
        org_study_hit.clear();
        p_sig_study_hit.clear();
        org_study_hit.resize( N_SIG_STUDY_HIT, 0 );
        p_sig_study_hit.resize( N_SIG_STUDY_HIT, 0 );

        org_study_gene_hit.clear();
        p_sig_study_gene_hit.clear();
        org_study_gene_hit.resize( N_SIG_STUDY_HIT, 0 );
        p_sig_study_gene_hit.resize( N_SIG_STUDY_HIT, 0 );


        std::map<int, std::set<Interval> > uniq_genes;
        for(unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
            std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
            int pway = 0;
            while ( t != targets.end() ) {
                // count the number of pathways with p-value less than a threshold
                if((double)(pcount[ pway ])/(double)(globals.nrep+1) < SIG_STUDY_HIT[ns]) {
                    org_study_hit[ns]++;

                    // retrieve the overlapping genes with all the sets
                    std::map< std::string, std::set< Interval > >::iterator itr = multi_test.begin();
                    while(itr != multi_test.end()) {
                        get_overlapping_genes( itr->second, t->second, uniq_genes[ns]);
                        itr++;
                    }
                }
                t++;
                pway++;
            }

            org_study_gene_hit[ns] = uniq_genes[ns].size();
        }
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRich::cal_study_sig_hit(" + int2str(N_SIG_STUDY_HIT) + ")\n";
        exit(1);
    }
}

void InRichMultiMain::init_second_permutation() {
    cal_study_sig_hit();
    //
    // We now have pcount[pathway] scored with 1..(R+1) hits
    // and we've also constructed the null distribution of scores
    // for each pathway (null_tracker[]).  Now we can loop through
    // a second time, repeating the exact same series of permutations,
    // and score the multiple test correction
    //

    try {
        pcorr.resize( n_set, 1 );
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: InRichMulti::init_second_permutation(" + int2str(n_set) + ")\n";
        exit(1);
    }
}

bool InRichMultiMain::run_second_permutation() {

    // Here, we have nrep Bootstrap samples
    for ( int p = 0 ; p < globals.nboot ; p++ )
    {
        std::map<Interval, int> dummy_match_tracker;

        std::map< std::string, std::set< Interval > > multi_matched;
        std::map<std::string, std::set<Interval> >::iterator ttr = multi_test.begin();
        while(ttr != multi_test.end()) {            
            // This is the assumed real sample
            if(globals.multi_mapfile) {
	    	if(globals.use_map) {
            		globals.scaffold = multi_scaffold[ttr->first];
	    	}

                multi_matched[ttr->first] = match( ttr->second,
                                                   multi_snps[ttr->first],
                                                   dummy_match_tracker,
                                                   globals.match_genes ? &all_genes_set : NULL  );
            }
            else {
                multi_matched[ttr->first] = match( ttr->second,
                                                   snps,
                                                   dummy_match_tracker,
                                                   globals.match_genes ? &all_genes_set : NULL  );
            }
            ttr++;
        }

        // Goal is to obtain the minimum p-value over all sets for this
        // permutation (i.e. assuming this was the actual dataset). We can use
        // null_*_tracker[] to facilitate this.
        std::map<std::string, std::map<int,int> > boot_org_count;
        std::vector<unsigned int> boot_study_hit(N_SIG_STUDY_HIT, 0);
        std::map<std::string, std::set<Interval> > ::iterator t = targets.begin();
        int pway = 0;

	// evalute the number of assumed real sets
        while ( t != targets.end() )
        {
            // calculate the overlapping number between the assumed bootstrap sample and pathway
            ttr = multi_test.begin();
            while(ttr != multi_test.end()) {
                boot_org_count[ttr->first][pway] = evaluate( multi_matched[ttr->first], t->second );

                ttr++;
            }
	    t++;
	    pway++;
	}

	// calculate the empirical p-values for all pathways for this replicate
        std::vector<int> boot_pcount;
	boot_pcount.resize( n_set, 0 );
        for ( int p2 = 0 ; p2 < globals.nrep ; p2++ )
        {
                // assume that I get one randomly matched sample each step
		int r = CRandom::rand(globals.nrep);

        	t = targets.begin();
		pway = 0;
        	while ( t != targets.end() )
        	{
			// for each pathway
                	ttr = multi_test.begin();
                	bool increase = false;
                	while(ttr != multi_test.end()) {
                                if( (int)multi_null_tracker[ttr->first][pway][r] >= boot_org_count[ttr->first][pway] ) {
                        		increase = true;
                        		break;
                    		}
                    		ttr++;
			}

                	if(increase) {
                    		boot_pcount[pway]++;
                	}

			t++;
			pway++;
                }
	}

	// find the minimum p-value among all pathway-level empirical p-values 	
        std::map<int, std::set<Interval> > cur_study_gene_hit;

        double min_test = globals.nrep;
       	t = targets.begin();
	pway = 0;
       	while ( t != targets.end() ) {
		if( min_test > boot_pcount[pway] ) {
			min_test = boot_pcount[pway];
		}

                // empirical p-value for this pathway
                double boot_emp_p = (double)(boot_pcount[pway]+1)/(globals.nrep+1);

                // counting the number of pathways with emp_p < 0.05, 0.01, 0.001
                for( unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
                	if( boot_emp_p < SIG_STUDY_HIT[ ns ] ) {
                    		boot_study_hit[ ns ]++;

                                // retrieve the overlapping genes with all the sets
                                std::map< std::string, std::set< Interval > >::iterator itr = multi_matched.begin();
                                while(itr != multi_matched.end()) {
                                    get_overlapping_genes( itr->second, t->second, cur_study_gene_hit[ns]);
                                    itr++;
                                }

                	}
            	}
		t++;
		pway++;
	}



        //
        // We now have, for this replicate, what the min p-value would have
        // been; use this to score pcorr
        //
       	t = targets.begin();
	pway = 0;
        min_test++;  // pcount[pway] starts from 1 to nboot+1
       	while ( t != targets.end() )
       	{
            if ( min_test <= pcount[pway] ) {
                pcorr[pway]++;
            }
            ++t;
            ++pway;
        }


        for( unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
            if( boot_study_hit[ ns ] >= org_study_hit[ ns ]) {
                p_sig_study_hit[ ns ]++;
            }

            if( cur_study_gene_hit[ns].size() >= org_study_gene_hit[ ns ]) {
                p_sig_study_gene_hit[ ns ]++;
            }
        }

        std::cout << p << " second-pass permutations         \r";
        std::cout.flush();
    }

    log(int2str(globals.nboot) + " second-pass permutations ( completed )\n");

    return true;
}



int InRichMultiMain::get_random_match_no( std::string _set_name, int _pway ) {
    // sample the overlapping number of a randomly matched set
    // based on the distribution in multi_null_tracker

    int r = CRandom::rand(globals.nrep);
    return multi_null_tracker[_set_name][_pway][r];
}

std::map<std::string, int> InRichMultiMain::get_random_match_no_list( int _pway ) {
    // sample the overlapping number of a randomly matched set
    // based on the distribution in multi_null_tracker

    int r = CRandom::rand(globals.nrep);
    std::map<std::string, int> no_list;
    std::map<std::string, std::set<Interval> >::iterator ttr = multi_test.begin();
    while(ttr != multi_test.end()) {
        no_list[ttr->first] = multi_null_tracker[ttr->first][_pway][r];
	ttr++;
    }

    return no_list;
}


bool InRichMultiMain::write_output() {

    log("\nTotal of " + int2str(n_interval) + " intervals tested for " + int2str(n_set) + " target sets\n\n");

    log_path_data();

    report_sig_study_hit( );
    log("\n\n\n" + sig_hit + "\n\n");

    std::cout << "Writing output to " + globals.outroot + ".out.inrich ... ";
    log_summary_file(globals.outroot + ".out.inrich");
    std::cout << " DONE\n\n";

    return true;
}

std::string InRichMultiMain::header(bool _line) {

    std::ostringstream mylog;
    std::string mode = globals.test == INTERVAL ? "Interval" : "Target" ;

    if(_line) {
        mylog << "\n----------------------------------------------------------------\n";
    }
    mylog << "N_Target\t";

    std::map<std::string, std::vector<int> >::iterator itr = multi_count.begin();
    while( itr != multi_count.end() ) {
        mylog << itr->first << "_" << mode << "_No\t";
        itr++;
    }

    mylog << "P\t"
            << "PCorr\t"
            << "Target\n";

    if(_line) {
        mylog << "\n----------------------------------------------------------------\n";
    }

    return mylog.str();
}


bool InRichMultiMain::log_path_data() {

    //
    // Verbose report of region stats:
    //
    std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
    int pway = 0;
    log(header());
    bool any_results = false;
    path_summary.clear();
    while ( t != targets.end() )
    {
        const std::set<Interval> & target = merge( t->second );

        //
        // Flag if pointwise empirical p-value is at minimum bound as the
        // corrected value might be an under-estimate in this case
        //

        double pvalue = pcount[ pway ]/(double)( 1 + globals.nrep );
        bool lowest_p = pcount[ pway ] == 1;
        double pvalue2 = pcorr[ pway ] / double( 1 + globals.nboot );

	if(pvalue2 < pvalue) {
		pvalue2 = pvalue;
	}

        std::string str = "";
        std::string mark = " ";
        if (lowest_p) {
            mark = "*";
        }

        str = mark + int2str(target.size()) +  "\t" ;

        std::map<std::string, std::vector<int> >::iterator itr = multi_count.begin();
        while( itr != multi_count.end() ) {
            str = str + int2str(itr->second[ pway ]) +  "\t";
            itr++;
        }
        str = str + double2str(pvalue) + "\t" +
              double2str(pvalue2) + mark + "\t" +
              t->first +  "\n";

        if ( pvalue <= globals.pthresh )
        {
            any_results = true;
            log(str);
        }
        path_summary.push_back(str);

        ++t;
        ++pway;
    }

    if ( ! any_results ) {
        log( "     { -- no significant results at p < " + double2str(globals.pthresh) +  " -- }\n") ;
    }


    return true;
}


void InRichMultiMain::log_summary_file(std::string _filename, bool _append) {
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


void InRichMultiMain::log_target_file(std::ofstream & OUT1){

    OUT1 << "T_TARG" << "\t";


    std::string mode = globals.test == INTERVAL ? "N_INT" : "N_TARG" ;
    std::map<std::string, std::vector<int> >::iterator itr = multi_count.begin();
    while( itr != multi_count.end() ) {
        OUT1 << mode << "_" << itr->first << "\t";
        itr++;
    }

    OUT1   << "P" << "\t"
            << "PCORR" << "\t"
            << "TARGET" << "\n";

    for(unsigned int i=0; i<path_summary.size(); i++) {
        OUT1 << path_summary[i];
    }
}


void InRichMultiMain::log_main_file(std::ofstream & OUT ) {
       OUT  << "SET" << "\t"
            << "INTERVAL" << "\t"
            << "GENE_LOC" << "\t"
            << "GENE_ID" << "\t"
            << "GENE_DESC" << "\t"
            << "TARGET" << "\t"
            << "P\n";

       std::map<std::string, std::vector<int> >::iterator itr = multi_count.begin();
       while( itr != multi_count.end() ) {
           log_main_data (OUT, itr->first);
           itr++;
       }
}


void InRichMultiMain::log_main_data(std::ofstream & OUT , std::string _set){

    //
    // Verbose report of region stats:
    //

    std::map<std::string, std::set<Interval> >::iterator t = targets.begin();
    int pway = 0;
    while ( t != targets.end() )
    {
        const std::set<Interval> & target = merge( t->second );

        // For this individual set, proportion of unplaced intervals
        //
        //
        // Get all genes in this target set that were in the test interval
        //
        std::map< Interval , std::set<Interval> > overlap = get_overlap( multi_test[_set], target );
        std::map< Interval , std::set<Interval> >::iterator ii = overlap.begin();

        //
        // Flag if pointwise empirical p-value is at minimum bound as the
        // corrected value might be an under-estimate in this case
        //

        double pvalue = pcount[ pway ]/(double)( 1 + globals.nrep );
        //bool lowest_p = pcount[ pway ] == 1;

        while ( ii != overlap.end() )
        {
            std::set<Interval> & o1 = ii->second;
            std::set<Interval>::iterator oo = o1.begin();
            while ( oo != o1.end() )
            {
                Interval tmp = ii->first;
                Interval tmp2 = *oo;

                OUT  << _set << "\t"
                     << tmp.get_desc() << "\t"
                     << tmp2.get_desc() << "\t"
                     << tmp2.name << "\t"
                     << all_genes[ tmp2.name ].desc << "\t"
                     << t->first << "\t"
                     << pvalue << "\n";

                ++oo;
            }
            ++ii;
        }

        ++t;
        ++pway;

    }
}

void InRichMultiMain::log_interval_file(std::ofstream & OUT2) {
    OUT2 << "SET" << "\t"
            << "INTERVAL" << "\t"
            << "N_SNP" << "\t"
            << "N_GENE" << "\t"
            << "GENE_ID" << "\t"
            << "GENE_DESC" << "\n";

    std::map<std::string, std::vector<int> >::iterator itr = multi_count.begin();
    while( itr != multi_count.end() ) {
        log_interval_data (OUT2, itr->first);
        itr++;
    }
}

void InRichMultiMain::log_interval_data(std::ofstream & OUT2, std::string _set) {
    //
    // List all interval/all-genes pairings, post de-duping, etc
    //
    interval_data.clear();

    std::map< Interval , std::set<Interval> > overlap_all = get_overlap( test , all_genes_set );
    std::set<Interval>::iterator k = test.begin();
    while ( k != test.end() )
    {
        /*
        std::map< Interval , std::vector<int> >::iterator ip = globals.acceptable_positions.find( *k );
        int npos = -1;
        if ( ip != globals.acceptable_positions.end() )
            npos = ip->second.size(); */

        std::set<Interval> & o2 = overlap_all[ *k ];
        Interval interval = *k;
        if ( o2.size() == 0 )
        {
            OUT2 << _set << "\t"
                    << interval.get_desc() << "\t"
                    << k->n << "\t"
                    << "0\t\t\n";

            //interval_data.push_back( Interval2Gene( interval.get_desc() , k->n, npos, (int)o2.size() ));
        }
        else
        {
            std::set<Interval>::iterator jj = o2.begin();
            jj = o2.begin();
            while ( jj != o2.end() )
            {
                OUT2 << _set << "\t"
                        << interval.get_desc() << "\t"
                        << k->n << "\t"
                        << (int)o2.size() << "\t"
                        << jj->name << "\t"
                        << all_genes[ jj->name ].desc << "\n";

                //interval_data.push_back( Interval2Gene( interval.get_desc() , k->n, npos, (int)o2.size(), jj->name , all_genes[ jj->name ].desc ));
                ++jj;
            }
        }
        ++k;
    }
}

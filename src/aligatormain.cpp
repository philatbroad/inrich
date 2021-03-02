#include "aligatormain.h"
#include "log.h"

#include <iostream>
#include <iomanip>

AligatorMain::AligatorMain()
{
    fake_name = false;
}

void AligatorMain::set_globals(InputData & _globals)
{
    globals = _globals;
}


AligatorMain::~AligatorMain()
{
}

void AligatorMain::make_cumul() {
        for(unsigned int p=0; p<n_targets; p++) {
		int t=0;
		std::map<int,int>::reverse_iterator i = rep_dist[p].rbegin();
		while ( i != rep_dist[p].rend() ) {
			t += i->second;
			i->second = t;
			++i;
		}

	}
}

void AligatorMain::aligator_test() {
    init_data();

    // ------------------------------------------------------
    // 1. Load Files
    // ------------------------------------------------------
    if(!load_files()) {
        log("ERROR: File not loaded\n");
        return;
    }

    // ------------------------------------------------------
    // 2. Pre-processing
    // ------------------------------------------------------
    if(!gen_snp2gene()) {
        log("ERROR: snp2gene not completed\n");
        return;
    }

    if(fake_name) {
        //log("  checking SNP Names ... \n");
        if(!update_fake_name()) {
            log("Checking names not completed");
            return;
        }

        //log("  generating sig gene set ...\n");
        if(!gen_sig_genes()) {
            log("Counting genes not completed\n");
            return;
        }
    }
    else {
        //log("  generating sig gene set ...\n");
        if(!gen_sig_genes( )) {
            log("ERROR: generating sig gene set not completed\n");
            return;
        }
    }

    if(n_sig_genes==0) {
        log("\n\nERROR: No significant genes were mapped to test!\n\n");
        return;
    }

    // ------------------------------------------------------
    // 3. Count the number of genes for each pathway
    // ------------------------------------------------------
    //log("\n\nTesting enrichment\n  counting number of genes for target sets ...\n");
    if(!cal_org_count( )) {
        log("Counting genes not completed");
        return ;
    }

    CRandom::srand( globals.seed );

    //log("  calculating empirical P-values ...\n");
    if(!first_permutation( )) {
        log("\nFirst-pass permutation not completed");
        return ;
    }

    make_cumul();


    //log("  calculating crrected P-values ... \n");
    if(!second_permutation( )) {
        log("Second-pass permutation not completed");
        return;
    }


    if(!write_output( )) {
        log("File writing not completed");
        return;
    }

}


bool AligatorMain::load_files() {
    int fileRead=0;
    int file_tried=0;

    file_tried++;
    if(load_file(file_tried, ALI_SNP_MAP )) {
        fileRead++;
    }

    file_tried++;
    if(load_file(file_tried, ALI_ASSOC )){
        fileRead++;
    }

    file_tried++;
    if(load_file(file_tried, ALI_GENE )){
        fileRead++;
    }

    file_tried++;
    if(load_file(file_tried, ALI_TARGET_SET )){
        fileRead++;
    }

    if( fileRead == file_tried ) {
        return true;
    }
    else {
        return false;
    }
}


bool AligatorMain::load_file(int _modeNum, int _mode) {

    int no=0;

    switch(_mode) {
    case ALI_ASSOC:
        no = load_snps(  "significant", globals.testfile, sig_snps );
        break;

    case ALI_SNP_MAP:
        no = load_snps(  "reference", globals.mapfile, all_snps );
        break;

    case ALI_GENE:
        no = load_genes_with_chr_genes( );
        break;

    case ALI_TARGET_SET:
        no = load_targets( );
        break;

    default:
        break;
    }

    if(no>=0) {
        return true;
    }
    else {
        return false;
    }
}



int AligatorMain::load_targets(  )
{

    std::ifstream TARGET( globals.targetfile.c_str() , std::ios::in );

    int cnt = 0;
    std::set<std::string> uniq_gene_not_found;

    while ( ! TARGET.eof() )
    {

        std::vector<std::string> lines = tokenizeLine( TARGET );
        if( lines.size()==0 ) break;
        if ( lines.size() < 2 ) continue;

        std::string genename = lines[0];
        std::string setname = lines[1];
        for (unsigned int i=2;i<lines.size(); i++) setname += " " + lines[i];

        std::map<std::string,Gene>::const_iterator i = all_genes.find( genename );

        if ( i == all_genes.end() )
        {
            uniq_gene_not_found.insert( genename );
            continue;
        }

        targets[setname].insert( genename );
        target_genes.insert( genename );

        ++cnt;
    }

    TARGET.close();

    log( "  read " + int2str( target_genes.size()) + " unique genes in "
         + int2str(targets.size()) + " targets\n"
         + "  " + int2str( cnt ) + " total gene/target pairs\n"
         + "  " + int2str( uniq_gene_not_found.size()) + " genes not found in reference gene-list\n" );

    // apply size filters?
    if ( globals.tsize_min || globals.tsize_max )
    {
        std::map<std::string, std::set<std::string> > targets2 = targets;
        targets.clear();
        std::map<std::string, std::set<std::string> >::iterator i = targets2.begin();
        while ( i != targets2.end() )
        {
            bool okay = true;
            if ( globals.tsize_min && i->second.size() < globals.tsize_min ) okay = false;
            if ( globals.tsize_max && i->second.size() > globals.tsize_max ) okay = false;
            if ( okay ) targets[ i->first ] = i->second;
            ++i;
        }

        log ( "  after size filters, "  + int2str(targets.size()) + " targets remaining\n") ;
    }

    n_targets = targets.size();

    return targets.size();
}


int AligatorMain::load_genes_with_chr_genes(   )
{

    globals.largest_gene = 0;

    // CHR BP1 BP2 NAME DESC

    //std::map<std::string, Gene > genes;
    if ( ! fileExists( globals.genefile ) )
    {
        std::cerr << "could not find gene list " << globals.genefile << "\n";
        exit(1);
    }

    std::ifstream GENES( globals.genefile.c_str() , std::ios::in );
    int cnt = 0;
    while (! GENES.eof() )
    {

        std::vector<std::string> lines = tokenizeLine( GENES );

        if( lines.size()==0 ) break;
        if ( lines.size() < 5 ) continue;

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp1, bp2;
        if ( ! str2int( lines[1] , bp1 ) ) continue;
        if ( ! str2int( lines[2] , bp2 ) ) continue;
        if( bp1 < 1 || bp2 < 1 ) continue;

        std::string name1 = lines[3];

        // are we filtering w/ a background list?
        if ( globals.background_list && background_genes.find( name1 ) == background_genes.end() )
            continue;

        // append rest as space-delim description
        std::string name2 = lines[4];
        for (unsigned int i=5; i<lines.size(); i++) name2 += " " + lines[i];


        bp1 -= globals.window;
        bp2 += globals.window;
        if ( bp1 < 1 ) bp1 = 1;
        if ( bp2 < bp1 ) bp2 = bp1;


        all_genes.insert( make_pair( name1, Gene( Interval( chr , bp1 , bp2 , name1 ) , name2 ) ) );

        chr_genes[chr].insert(Interval(chr, bp1, bp2, name1));

        ++cnt;
    }

    // get 1 larger than largest gene
    ++globals.largest_gene;

    GENES.close();


    log( "  read " + int2str(cnt) + " reference genes\n");

    return all_genes.size();
}

int AligatorMain::load_snps( std::string _snp_type, std::string _file, std::set<SNP> & snps )
{
    std::map<int,int> min;
    std::map<int,int> max;


    std::ifstream MAP( _file.c_str() , std::ios::in );
    int i = 0;
    int snp_no=0;
    while ( ! MAP.eof() )
    {
        i++;

	if(globals.topn>0 && globals.topn<i) break;

        std::vector<std::string> lines = tokenizeLine( MAP );
        if(lines.size()==0) break;
        if(lines.size()<2) {
            continue;
        }

        int chr = get_chr(lines[0], globals.is_human);
        if( ! is_valid_chr(chr, globals.is_human) ) continue;

        int bp;
        if( ! str2int(lines[1], bp) ) continue;
        if ( bp < 1 ) continue;

        std::string name;
        double p;

        if(lines.size()==2) {
            snp_no++;
            name = "snp" + int2str(snp_no);
            snps.insert( SNP (  chr, bp, name ) );

            if( i==1 ) {
                fake_name = !fake_name;
            }
        }
        else if(lines.size()==3) {
            snps.insert( SNP (  chr, bp, lines[2] ) );
        }
        else if(lines.size()==4) {
            if(str2double(lines[3], p)) {
                snps.insert( SNP (  chr, bp, lines[2], p ) );
            }
            else {
                snps.insert( SNP (  chr, bp, lines[2] ) );
            }
        }

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

	if(i%50000==0) {
		std::cerr << "Reading " << i << " SNP positions\r";
		std::cout.flush();
	}
    }


    MAP.close();

    if(_snp_type.compare("reference")==0) {
        globals.total_seq = 0;
        globals.scaffold.clear();
        std::map<int,int>::iterator i = min.begin();
        while ( i != min.end() )
        {
            globals.scaffold[ i->first ] = Interval( i->first , i->second , max[ i->first ] );
            globals.total_seq += (long unsigned int)(max[ i->first ] - i->second + 1);
            ++i;
        }

        log( double2str(globals.total_seq/(double)1e6) + "Mb total sequence length on "
             + int2str(globals.scaffold.size()) + " chromosomes\n");
    }

    log( "  read " + int2str(snps.size()) + " " + _snp_type + " SNPs\n");


    return snps.size();
}


bool AligatorMain::update_fake_name( ) {
    int cnt = 0;
    std::set<SNP>::iterator itr = sig_snps.begin();
    while ( itr != sig_snps.end() )
    {
        std::set<SNP>::iterator jtr = all_snps.find(*itr);
        if(jtr != all_snps.end()) {
            snp2genes[itr->name] = snp2genes[jtr->name];
        }
        ++cnt;
        ++itr;
    }
    return true;
}


bool AligatorMain::gen_snp2gene(  )
{
    //_msg = "temporary set to null";

    int cnt = 0;
    int uniq_gene_no = 0;
    int gene_no_in_targets = 0;

    // This is done for all SNPs. Note that SNP is ordered by their chromosomal location.
    std::set<SNP>::iterator itr = all_snps.begin();
    int x_chr = -1;
    std::set<Interval> chr_gene_list;
    std::set<Interval>::iterator jtr;
    bool start = false;
    std::string x_snp_id;
    int chr, bp;
    int snps2genes_no = 0;
    while ( itr != all_snps.end() )
    {


        // chromosomal location
        chr = itr->chr;
        bp = itr->bp;

        // retrieve the list of chromosome genes if this is the first SNP on the chromosome
        if( x_chr != chr ) {
            chr_gene_list = chr_genes[chr];
            jtr = chr_gene_list.begin();
            start = true;
        }
        else {
            start = false;
        }

        // if there are proceding SNPs, examine this SNP also matches to the previously examined genes
        if(!start) {
            if(snp2genes.end() != snp2genes.find(x_snp_id)) {
                std::vector<Gene> prev_genes = snp2genes[x_snp_id];
                for( unsigned int p=0; p< prev_genes.size(); p++ ) {
                    Gene ptr = prev_genes.at(p);
                    if(ptr.gene.bp1<=bp && ptr.gene.bp2>=bp) {
                        //snp2genes.insert(make_pair(itr->name, ));
                        snp2genes[itr->name].push_back(ptr);
                        snps2genes_no++;
                    }
                }
            }

        }

        while ( jtr != chr_gene_list.end() )
        {
            // if a SNP matches to the gene
            if(jtr->bp1 > bp) {
                break;
            }
            if(jtr->bp1 <= bp && jtr->bp2 >= bp) {
                snp2genes[itr->name].push_back( all_genes[jtr->name] );
                snps_in_genes.push_back(itr->name);
                uniq_gene_no++;
                snps2genes_no++;

                if(target_genes.find(jtr->name) != target_genes.end()) {
                    gene_no_in_targets++;
                }
            }

            ++jtr;
        }

        x_snp_id = itr->name;
        x_chr = chr;

        ++cnt;
        ++itr;
    }

    log( "\n  " + int2str( uniq_gene_no ) + " genes mapped from "
         + int2str( snps2genes_no ) + " reference SNPs\n  "
         + int2str( gene_no_in_targets ) + " genes found in targets\n" );


    return true;
}


bool AligatorMain::gen_sig_genes( )
{

    // This is done for all significant SNPs
    std::set<SNP>::iterator itr = sig_snps.begin();
    while ( itr != sig_snps.end() )
    {

        // identify mapping genes
        // std::map<std::string, std::vector<Gene> >::const_iterator jtr = snp2genes.find(itr->name);
        std::map<std::string, std::vector<Gene> >::const_iterator jtr = snp2genes.find(itr->name);
        if ( jtr != snp2genes.end() )
        {
            std::vector<Gene> genes = jtr->second;
            for(unsigned int i=0; i<genes.size(); i++) {
                sig_genes[ genes[i].gene.name ].insert( *itr );
            }
        }
        else {
            // This is problematic!!!
        }

        ++itr;
    }

    n_sig_genes = sig_genes.size();

    log( "  " + int2str(sig_snps.size()) + " significant SNPs mapped to "
         + int2str( n_sig_genes ) + " genes\n" );


    return true;
}

void AligatorMain::init_data() {
    org_count.clear();
    rep_dist.clear();
    emp_p.clear();
    corr_p.clear();
    org_study_hit.clear();
    p_sig_study_hit.clear();
    exp_hit_count.clear();
}


bool  AligatorMain::cal_org_count( ) {


    std::map<std::string, std::set<std::string> >::iterator itr = targets.begin();
    int pway = 0;

    try {
        org_count.resize( n_targets, 0);
	rep_dist.resize( n_targets );


        // for each pathway
        while (itr != targets.end()) {
            org_count[ pway ] = count_overlap( sig_genes, itr->second );
	    rep_dist[ pway ] [ org_count[ pway] ]++;

            ++itr;
            ++pway;
        }


        log( "\nTotal of " + int2str(sig_genes.size()) + " genes tested for " + int2str(n_targets) + " target sets\n\n");

    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: AligatorMain::cal_org_count(" + int2str(n_targets) + ")\n";
        exit(1);
    }

    return true;
}

int AligatorMain::count_overlap(  std::map<std::string, std::set<SNP> > & genes,
                                  std::set<std::string> & pathway_genes )
{
    std::set<std::string>::iterator itr = pathway_genes.begin();

    int cnt = 0;
    while ( itr != pathway_genes.end() ) {

        /*
        std::map<std::string, std::set<SNP> >::iterator jtr = genes.begin();
        while( jtr !=  genes.end() ) {
            if( (*itr).compare( jtr->first )==0 ) {
                //std::string symbol = find_name( all_genes, j->gene);
                cnt++;
                break;
            }
            ++jtr;
        } */
        std::map<std::string, std::set<SNP> >::iterator jtr = genes.find(*itr);
        if( jtr != genes.end() ) {
            cnt++;
        }

        ++itr;
    }

    return cnt;
}

std::map<std::string, std::set<SNP> > AligatorMain::gen_random_set( ) {

    std::map<std::string, std::set<SNP> > rand_list;

    //std::map<std::string, std::vector<Gene> > snp2genes;
    //std::vector<std::string> snps_in_genes;


    int snp_no = snps_in_genes.size();
    std::set<SNP> dummy_snps;
    while( rand_list.size() < n_sig_genes ) {
        int r = CRandom::rand( snp_no );

        std::vector<Gene> genes = snp2genes[ snps_in_genes[r] ] ;
        for (unsigned int i=0; i<genes.size(); i++) {
            rand_list[ genes[i].gene.name ] = dummy_snps;

            if(rand_list.size() == n_sig_genes) {
                break;
            }
        }
    }

    return rand_list;
}



bool AligatorMain::first_permutation( ) {

    try {
        emp_p.resize( n_targets );

        unsigned int pway;


        for( int i=0; i < globals.nrep; i++) {
            std::map<std::string, std::set<SNP> > rand_set = gen_random_set();

            std::map<std::string, std::set<std::string> >::iterator itr = targets.begin();
            pway = 0;
            while (itr != targets.end()) {
                if( i==0 ) {
                    emp_p[ pway ] = 1;
                }

                unsigned int rand_count = count_overlap( rand_set, itr->second );

                if ( org_count[ pway ] < globals.min_cnt) {
                    if( i==0 ) {
                        emp_p[ pway ] = globals.nrep + 1;
                    }
                }
                else if ( rand_count >= org_count[ pway ] ) {
                    ++emp_p[ pway ];
                }

                rep_dist[ pway ][ rand_count ]++;

                ++itr;
                ++pway;
            }


                std::cerr << i << " first-pass permutations \r" ;
                std::cout.flush();

        }

        double d = globals.nrep + 1;
        org_study_hit.resize( N_SIG_STUDY_HIT, 0 );
        for(unsigned int pway=0; pway<n_targets; pway++) {
            emp_p[ pway ] /=  d;

            for(unsigned  int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
                if( emp_p[ pway ] < SIG_STUDY_HIT[ ns ] ) {
                    org_study_hit[ ns ]++;
                }
            }
        }
        log(  int2str(globals.nrep) + " first-pass permutations (completed) \n" );

    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: AligatorMain::first_permutation(" + int2str(globals.nrep) + ")\n";
        exit(1);
    }
    return true;
}





void AligatorMain::log(std::string _msg) {
    Log mylog(globals.outroot + ".out.aligator");
    mylog << _msg;
}

std::vector<int>  AligatorMain::gen_random_sets( int _observed ) {
    std::vector<int> rand_sets;

    while( (int)rand_sets.size() < globals.nrep ) {
        int r = CRandom::rand( globals.nrep );
        if( r != _observed ) {
            rand_sets.push_back(r);
        }
    }

    return rand_sets;
}


bool AligatorMain::second_permutation( ) {

    try {
        exp_hit_count.resize( n_targets );
        corr_p.resize( n_targets );
        p_sig_study_hit.resize( N_SIG_STUDY_HIT );

        int pway;


        for( int i=0; i < globals.nboot; i++) {

            // This is the replicate random gene list selected to be the ``observed'' data.
            //int fake_observed = CRandom::rand( globals.nrep );
            //std::vector<int> rand_sets = gen_random_sets( fake_observed );
	    
            std::map<std::string, std::set<SNP> > rand_set = gen_random_set();

            std::vector<int> cur_study_hit(N_SIG_STUDY_HIT, 0);

            std::map<std::string, std::set<std::string> >::iterator itr = targets.begin();



            // Goal is to obtain the minimum p-value over all sets for this
            // permutation (i.e. assuming this was the actual dataset).

            double min_p = 1;
            pway = 0;

            // ammong all pathway p-values, select the minimum
            while (itr != targets.end()) {
		 const int count = count_overlap( rand_set, itr->second );


                // calculate the empirical p-value for this pathway
                double cur_emp_p = 1;

                // if number of sig genes is less than globals.min_cnt, automatically set p-value=1
                // this doesn't change min_test
                // so, we skip
                if( count < (int)globals.min_cnt ) {
                    cur_emp_p = 1;
                }
                else {
               	    cur_emp_p = (double)(rep_dist[ pway ][ count ]) / (double)(globals.nrep+1);
                }

                // update the min count
                if( cur_emp_p < min_p ) {
                    min_p = cur_emp_p;
                }

                // update the expected number of categories  with at least this p-value
                for(unsigned  int p=0; p< n_targets; p++) {
                    if( cur_emp_p <= emp_p[p]  ) {
                        exp_hit_count[ p ]++;
                    }
                }

                // update the number of categories
                for(unsigned  int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
                    if( cur_emp_p < SIG_STUDY_HIT[ ns ] ) {
                        cur_study_hit[ ns ]++;
                    }
                }

                ++itr;
                ++pway;
            }


            // for all original targets
            // update a corrected p-value
            for( unsigned int pway=0; pway < n_targets; pway++) {
                if( min_p <= emp_p[ pway ] ) {
                    corr_p[ pway ]++;
                }
            }

            for( unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++) {
                if( cur_study_hit[ ns ] >= (int)org_study_hit[ ns ]) {
                    p_sig_study_hit[ ns ]++;
                }
            }


            std::cerr << i << " second-pass permutations \r" ;
            std::cout.flush();

        }

        // calculate the expected number of hits for each category
        // as the average number of categories per bootstrap replicate with p values less than or euqal to the uncorrected p-value from the original data
        for( unsigned int pway=0; pway < n_targets; pway++) {
            exp_hit_count[ pway ] = ( exp_hit_count[ pway ]+1 ) / ( globals.nboot + 1);
            corr_p[ pway ] = ( corr_p[ pway ]+1 ) / ( globals.nboot + 1);
        }
        for( unsigned int ns=0; ns < N_SIG_STUDY_HIT; ns++ ) {
            p_sig_study_hit[ ns ] = ( p_sig_study_hit[ ns ]+1 ) / ( globals.nboot + 1);
        }

        log( int2str(globals.nboot) + " second-pass permutations (completed) \n" );
    }
    catch(std::bad_alloc a) {
        std::cerr << "Memory Allocation Error: AligatorMain::second_permutation(" + int2str(globals.nboot) + ")\n";
        exit(1);
    }

    return true;
}


bool AligatorMain::write_output( ) {
    summary_filename = globals.outroot + ".out.aligator";

    log("\n\nWriting output to : " + summary_filename + "\n\n");

    if(!log_target_path_data( )) {
        return false;
    }

    if(!log_snp_data( )) {
        return false;
    }

    report_sig_study_hit();

    log_summary_file();


    return true;
}


void AligatorMain::log_summary_file(){

    std::ofstream OUT( (summary_filename).c_str() , std::ios::app );

    OUT << "----------------------------------------------------------------\n";
    OUT << "I. Main Analysis Results\n";
    OUT << "----------------------------------------------------------------\n";
    log_target_file(OUT);

    OUT << "\n\n\n";
    OUT << "----------------------------------------------------------------\n";
    OUT << "II. Target-Gene-SNP Summary\n";
    OUT << "----------------------------------------------------------------\n";
    log_main_file(OUT);


    OUT << "\n\n\n";
    OUT << "----------------------------------------------------------------\n";
    OUT << "III. Significant Gene Summary\n";
    OUT << "----------------------------------------------------------------\n";
    log_gene_file(OUT);
}


void AligatorMain::log_target_file(std::ofstream & OUT1){
    OUT1 << "O1_" << "\t"
            << "TAR_SIZE" << "\t"
            << "SIG_GENE" << "\t"
            << "EXP_GENE" << "\t"
            << "P" << "\t"
            << "PCORR" << "\t"
            << "EXP_HIT" << "\t"
            << "TARGET" << "\n";

    for(unsigned int i=0; i<path_data.size(); i++) {
        OUT1 << "O1_" << "\t" << path_data.at(i).get_aligator_summary();
    }
}


void AligatorMain::log_gene_file(std::ofstream & OUT2) {
    OUT2 << "O3_" << "\t" << "GENE_ID" << "\t"
            << "GENE_LOC" << "\t"
            << "GENE_DESC" << "\t"
            << "SNP" << "\t"
            << "SNP_LOC" << "\t"
            << "SNP_P" << "\n";

    for(unsigned int i=0; i<snp_data.size(); i++) {
        OUT2 << "O3_" << "\t" << snp_data.at(i).get_summary();
    }
}

void AligatorMain::log_main_file(std::ofstream & OUT3){
    OUT3 << "O2_" << "\t" << "TARGET" << "\t"
            << "P" << "\t"
            << "GENE_ID" << "\t"
            << "GENE_LOC" << "\t"
            << "GENE_DESC" << "\t"
            << "SNP" << "\n";

    for(unsigned int i=0; i<main_data.size(); i++) {
        OUT3 << "O2_" << "\t" << main_data.at(i).get_aligator_summary();
    }
}



bool AligatorMain::log_target_path_data( ) {
    //
    // Verbose report of region stats:
    //
    if( main_data.size() > 0 ) {
    	main_data.clear();
    }
    if( path_data.size() > 0 ) {
    	path_data.clear();
    }

    std::map<std::string, std::set<std::string> >::iterator target = targets.begin();
    int pway = 0;
    double min_p = 1/(globals.nrep+1);
    bool any_results = false;
    log_header();

    std::ostringstream myLog;
    double sig_ratio = (double)sig_genes.size() / (double)all_genes.size();

 try {
    while ( target != targets.end() )
    {

        bool lowest_p = emp_p[ pway ] == min_p;

        // Write to main output file
        std::set<std::string>  genes = target->second;
        std::set<std::string>::iterator gene = genes.begin();
        std::string p = double2str(emp_p[ pway ]);
        std::string target_name = target->first;

        while ( gene != genes.end() )
        {
	    std::map<std::string, Gene>::iterator gtr = all_genes.find(*gene);
	    std::map<std::string, std::set<SNP> >::iterator str = sig_genes.find(*gene);
	    if(gtr != all_genes.end() && str!=sig_genes.end()) {
            Gene cur_gene = all_genes[*gene];
            std::string gene_id = *gene;
            std::string gene_loc = cur_gene.get_loc();
            std::string gene_desc = cur_gene.get_desc();

            std::set<SNP> snps = sig_genes[*gene];
            std::set<SNP>::iterator snp = snps.begin();
            while( snp != snps.end() ) {
                SNP cur_snp = *snp;
                std::string snp_desc = cur_snp.get_desc();

                main_data.push_back( AllRes(target_name,
                                               p,
                                               gene_id,
                                               gene_loc,
                                               gene_desc,
                                               snp_desc) );
                ++snp;
            }
	    }

            ++gene;
        }

        double exp_gene_no = ((double)target->second.size()) * sig_ratio;
	std::string dummy = "";
        path_data.push_back( PathwayRes ( int2str(target->second.size()),
                                             int2str( org_count[ pway ]),
                                             emp_p[ pway ],
                                             corr_p[ pway ],
                                             int2str(lowest_p),
					     dummy,
                                             target_name,
                                             exp_hit_count[ pway ],
                                             exp_gene_no) );

        if ( emp_p[ pway ] <= globals.pthresh )
        {
            any_results = true;

            std::string mark = " ";
            if (lowest_p) {
                mark = "*";
            }

            myLog << std::setw(6) << target->second.size()
                  << std::setw(12) << org_count[ pway ]
                  << std::setw(15) << emp_p[ pway ]
                  << std::setw(12) << corr_p[ pway ] << mark
                  << std::setw(4) << " " << target_name <<  "\n" ;
        }
        //
        // Output for next target group
        //

        ++target;
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
catch(std::exception & e) {
	std::cerr << "Exception log_target_path_data " << e.what() << "\n";
	std::cerr << " pway=" << pway << "\n\n";
	return false;
}

}



bool AligatorMain::log_snp_data( ) {
    snp_data.clear();

    std::map<std::string, std::set<SNP> >::iterator itr = sig_genes.begin();
    while ( itr != sig_genes.end() )
    {   std::set<SNP> snps = itr->second;
        std::set<SNP>::iterator jtr = snps.begin();
        while( jtr != snps.end() ) {
            Gene gene = all_genes[itr->first];
            SNP snp = *jtr;
            snp_data.push_back( SigGene(itr->first, gene.get_chr_loc(), gene.get_desc(), snp.name, snp.get_loc(), snp.p ));
            ++jtr;
        }

        ++itr;
    }

    return true;
}

void AligatorMain::log_header() {

    std::ostringstream mylog;

    mylog << "\n----------------------------------------------------------------\n"
            << std::setw(6) << "T_Size"
            << std::setw(12) << "Sig_Gene_No"
            << std::setw(15) << "Empirical_P"
            << std::setw(12) << "Corrected_P"
            << std::setw(12) << "Target\n";

    mylog << "----------------------------------------------------------------\n";

    log(mylog.str());
}


void AligatorMain::report_sig_study_hit()
{
    std::ostringstream myLog;

    myLog << "\n----------------------------------------------------------------\n"
            <<  std::setw(18) << "Target_P_Threshold " 
            <<  std::setw(28) << "Number_of_Targets"  
            <<  std::setw(18) << "Significance" << "\n"
            << "----------------------------------------------------------------\n";

    for(unsigned int i=0; i<N_SIG_STUDY_HIT; i++) {
        myLog   <<  std::setw(18) << SIG_STUDY_HIT[ i ] 
                <<  std::setw(28) << org_study_hit[ i ]
                <<  std::setw(18) << p_sig_study_hit[ i ] << "\n";
    }

    myLog << "\n\n";

    sig_hit = myLog.str();

    log("\n\n\n" + sig_hit + "\n\n");
}

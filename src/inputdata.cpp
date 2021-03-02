#include "inputdata.h"
#include "crandom.h"
#include "helper.h"
#include "getopt.h"
#include "log.h"

#include <ostream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <time.h>

InputData::InputData()
{
    verbose = false;
    precompute = true;
    is_snps = false;
    tsize_min = 2;
    topn = 0;
    tsize_max = 0;
    test = INTERVAL;
    largest_gene = 5000000;
    match_genes = true;
    window = 0;
    nrep = -1;
    nboot = -1;

    printPermutations = false;

    srand( (unsigned)time ( NULL ) );
    seed = rand();

    pthresh = 0.05;
    tol = 0.1;
    use_map = false;   // =======================>
    use_range = false;
    background_list = false;
    min_cnt = 2;
    multi_set = false;
    multi_mapfile = false;
    //multi_testmode = 0;

    is_human = true;

    exam_bp_dist = false;
    exam_bp_dist_top = 0;

    mapfile = "";
    testfile = "";
    genefile = "";
    targetfile = "";
    rangefile = "";
    backgroundfile = "";
    create_db_name = "";
    use_db_name = "";
    outroot = "";

    compact = true;
    aligator = false;
    //create_lddb = false;
    //use_lddb = false;
    convert_mode = false;
}



void InputData::set_data( int argc, char* argv[] ) {
    parse_arg( argc, argv );
    check_arg();
    set_arg();

    if(aligator) {
        Log mylog( outroot + ".out.aligator", true );
	mylog << get_aligator_cmd();
    }
    else {
    	log_arg();
    }
}

int InputData::parse_arg( int argc, char* argv[] ) {

    int c = -1;
    opterr = 0;

    while ((c = getopt (argc, argv, "2a:b:cd:ef:g:h:i:j:klm:n:o:p:q:r:t:uvw:yx:z:s")) != -1) {
        switch (c)
        {
        case '2':
            test = TARGET;
            break;

        case 'a':
            testfile = optarg;
            break;

        case 'b':
            backgroundfile = optarg ;
            background_list = true;
            break;

        case  'c':
            precompute = false;
            break;

        case 'd':
            tol = atof( optarg );
            break;

        case 'e':
            match_genes = false;
            break;

        case 'f':
            seed = atoi(optarg);
            break;

        case 'g':
            genefile = optarg;
            break;

        case 'h':
            exam_bp_dist = true;
            exam_bp_dist_top = atoi( optarg );
            break;

        case 'i':
            tsize_min = atoi( optarg );
            break;

        case 'j':
            tsize_max = atoi( optarg );
            break;

        case 'k':
            compact = false;
            break;

        case 'l':        
            is_snps = true;
            aligator = true;
            //nboot = atoi( optarg );
            break;

        case 'm':
            use_map = true;
            mapfile = optarg;
            multi_mapfile = is_multi_mapfile();
            break;

        case 'n':
            topn = atoi( optarg);
            break;

        case 'o':
            outroot = optarg;
            break;

        case 'p':
            pthresh = atof( optarg );
            break;

        case 'q':
            nboot = atoi(optarg);
            break;

        case 'r':
            nrep = atoi(optarg);
            break;

        case 's':
        	printPermutations = true;
            break;

        /*
        case 's':
            is_snps = true;
            convert_mode = true;
            use_lddb = true;
            use_db_name = optarg;

            if(use_db_name.find(".db")==std::string::npos) {
                create_lddb = true;
                use_lddb = false;
                create_db_name = use_db_name;
                use_db_name = create_db_name + ".db";
            }

            break;
         */


        case 't':
            targetfile = optarg;
            break;

        case 'u':
            is_human = false;
            break;

        case 'v':
            verbose = true;
            break;

        case 'w':
            window = atoi( optarg );
            break;

        case 'x':
            rangefile = optarg;
            break;

        case 'y':
            multi_set = true;
            //multi_testmode = atoi( optarg );
            break;

        case 'z':
            min_cnt = atoi(optarg);
            break;

        case '?':
            if ( optopt == 'a' ||
                 optopt == 'b' ||
                 optopt == 'f' ||
                 optopt == 'g' ||
                 optopt == 't' ||
                 optopt == 'h' ||
                 optopt == 'i' ||
                 optopt == 'j' ||
                 optopt == 'm' ||
                 optopt == 'x' ||
                 optopt == 'o' ||
                 optopt == 'q' ||
                 optopt == 'u' ||
                 optopt == 'w' ||
                 optopt == 's' ||
                 optopt == 'r' ||
                 optopt == 'p' ||
                 optopt == 'd' ||
                 optopt == 'z' ||
                 optopt == 'y' ||
                 optopt == 'n'
                 )
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr,
                         "Unknown option character `\\x%x'.\n",
                         optopt);
            return 1;

        default:
            abort ();
        }
    }


    if(outroot.compare("")==0) {
        outroot = "test";
    }

    return 0;
}

int InputData::check_arg() {

    // --------------------------------
    // Checking up essential parameters
    // --------------------------------
    if ( testfile == "" ) { std::cerr << "no associated region file specified, -a\n"; exit(1); }

    if( convert_mode ) {
    	if ( mapfile == "" ) { std::cerr << "no map file specified, -m\n"; exit(1); }
    }
    else {
    	if ( targetfile == "" ) { std::cerr << "no target file specified, -t\n"; exit(1); }
    	if ( genefile == "" ) { std::cerr << "no gene file specified, -g\n"; exit(1); }

    	use_map = mapfile != "";
    	use_range = rangefile != "";
    	if ( ! ( use_range || use_map ) )
    	{
            std::cout << " --- halting --- must specify either -x {ranges} or -m {snp-map} ---\n";
            exit(1);
    	}
        if ( use_range && use_map )
        {
            std::cout << " --- halting --- must specify only one of -x {ranges} or -m {snp-map} ---\n";
            exit(1);
        }
    }

    if(nrep==-1) {
        nrep=5000;
    }
    if(nboot==-1) {
	nboot = nrep;
    }

    return 0;
}

void InputData::set_arg() {
    if ( seed == 0 ) {
        seed = time(0);
    }

    CRandom::srand( seed );

    if(multi_set) {
	if(nboot==0) {
            nboot = nrep;
	}
    }
}

void InputData::log_arg() {
    Log mylog( outroot + ".out.inrich", true );

    mylog << "----------------------------------------------------------------\n"
            << inrich_ver  << " : " + cur_date()
            << "http://atgu.mgh.harvard.edu/inrich\n"
            << "----------------------------------------------------------------\n";

    if(convert_mode) {
        mylog   << "project-title  (-o)  :  " << outroot << "\n"
                << " test-regions  (-a)  :  " << testfile << "\n"
                << "     map-file  (-m)  :  " << ( use_map ? mapfile : "--no map--" ) << "\n"
                << "    ld-db-file (-s)  :  " << use_db_name << "\n"
          	<< "----------------------------------------------------------------\n\n";
    }
    else {
        mylog  << "       project-title  (-o)  :  " << outroot << "\n"
                << "        test-regions  (-a)  :  " << testfile << "\n"
                << "           gene-list  (-g)  :  " << genefile << "\n"
                << "     background-genes (-b)  :  " << ( background_list ? backgroundfile : "--no background set--" ) << "\n"
                << "          range-file  (-x)  :  " << ( !use_range ? "--no ranges--" : rangefile ) << "\n"
                << "         target-file  (-t)  :  " << targetfile << "\n"
	        << "         compact      (-k)  :  " << ( compact ? "YES" : "NO" ) << "\n";
        mylog << " target size-filter  (-i,j) :  ";
        if ( tsize_min ) mylog << tsize_min; else mylog << 1;
        mylog << "..";
        if ( tsize_max ) mylog << tsize_max; else mylog << "all";;
        mylog << "\n";
        mylog << "    min-obs threshold (-z)  :  " << min_cnt << "\n"
                << "           test-type  "
                << ( test == TARGET ? "(-2)  :  TARGETS\n" : "(-1)  :  INTERVALS\n");
        if ( topn == 0 ) mylog << "       top-N-regions  (-n)  :   --all--\n";
        else mylog << "       top-N-regions  (-n)  :  " << topn << "\n";
        mylog << "           kb-window  (-w)  :  " << window << "\n";

        if ( multi_mapfile ) {
            mylog  << "      multi-map-file  (-m)  :  " << mapfile << "\n";
        }
        else {
            mylog  << "            map-file  (-m)  :  " << ( use_map ? mapfile : "--no map--" ) << "\n";
        }

        if ( multi_set ) {
            mylog << "      multi-set-mode  (-y)  :  yes\n";
        }

        if ( use_map ) mylog << "       match-density  (-d)  :  " << tol << "\n";
        else mylog << "       match-density  (-d)  :  " << "--no map--\n";

        if(	exam_bp_dist ) mylog << "     pairwise-distance (-h) :  " << exam_bp_dist_top << "\n";

        mylog   << "         pre-compute  (-c)  :  " << ( precompute ? "YES" : "NO" ) << "\n"
                << "         match-genes  (-e)  :  " << ( match_genes ? "YES" : "NO" ) << "\n"
                << "      num-replicates  (-r)  :  " << nrep << "\n"
                << "      num-bootstraps  (-q)  :  " << nboot << "\n"
                << "         random-seed  (-f)  :  " << seed << "\n"
                << "           display-p  (-p)  :  " << pthresh << "\n"
                << "   printPermutations  (-s)  :  " << ( printPermutations ? "YES" : "NO" ) << "\n"
                << "----------------------------------------------------------------\n\n";

    }
}

bool InputData::is_multi_mapfile( ) {
    if ( ! fileExists( mapfile ) )
    {
        std::cerr << "could not find map list file " << mapfile << "\n";
        exit(1);
    }

    std::ifstream MAP( mapfile.c_str() , std::ios::in );
    while ( !MAP.eof() ) {
        std::vector<std::string> lines = tokenizeLine( MAP );
	if( lines.size()<2 ) break;

        if( is_all_numeric(lines[1]) ) {
            // if the second column is base bp, this is not a multi map file
            return false;
        }
    }
    MAP.close();

    return true;
}

std::string InputData::get_cmd() {
    if(aligator) {
        return get_aligator_cmd();
    }
    else {
        return get_inrich_cmd();
    }
}


std::string InputData::get_qt_cmd() {
    if(aligator) {
        return get_aligator_qt_cmd();
    }
    else {
        return get_inrich_qt_cmd();
    }
}

std::string InputData::get_aligator_cmd() {

    std::ostringstream ss;

    ss      << "----------------------------------------------------------------\n"
            << "ALIGATOR v1.0 : " + cur_date()
            << "Reference : Peter Holmans et al. 2010. Am J Hum Genet  \n"
            << "            http://x004.psycm.uwcm.ac.uk/~peter/\n"
            << "Program   : http://atgu.mgh.harvard.edu/inrich\n"
            << "----------------------------------------------------------------\n"
            << "        Project Title  (-o) :  " << outroot << "\n"
            << "             SNP File  (-a) :  " << testfile << "\n"
            << "             Map File  (-m) :  " << ( use_map ? mapfile : "--no map--" ) << "\n"
            << "            Gene File  (-g) :  " << genefile << "\n"
            << "          Target File  (-t) :  " << targetfile << "\n"
            << "      Background File  (-b) :  " << ( background_list ? backgroundfile : "--no range--" ) << "\n";

    ss      << "Target Size Filter  (-i,-j) :  ";
    if ( tsize_min ) ss << tsize_min; else ss << 1;
    ss      << "..";
    if ( tsize_max ) ss << tsize_max; else ss << "all";
    ss << "\n";
    ss      << "     Min-obs Threshold (-z) :  " << min_cnt << "\n";

    ss      << "            Kb-Window  (-w) :  " << window << "\n"
            << "            Display-P  (-p) :  " << pthresh << "\n";


    ss      << "       Num Replicates  (-r) :  " << nrep << "\n"
            << "       Num Bootstraps  (-q) :  " << nboot << "\n"
            << "          Random Seed  (-f) :  " << seed << "\n";

    ss <<  "----------------------------------------------------------------\n\n";


    return ss.str();
}

std::string InputData::get_inrich_cmd() {

    std::ostringstream ss;



    ss      << "----------------------------------------------------------------\n"
            << inrich_ver << " : " + cur_date()
            << "http://atgu.mgh.harvard.edu/inrich\n"
            << "----------------------------------------------------------------\n";

    if(convert_mode) {
        ss      << "Project Title  (-o) :  " << outroot << "\n"
                << " Test Regions  (-a) :  " << testfile << "\n"
                << "     Map File  (-m) :  " << ( use_map ? mapfile : "--no map--" ) << "\n"
                << "    LD DB File (-s) :  " << use_db_name << "\n";
    }
    else {
        ss      << "        Project Title  (-o) :  " << outroot << "\n"
                << "         Test Regions  (-a) :  " << testfile << "\n"
                << "             Map File  (-m) :  " << ( use_map ? mapfile : "--no map--" ) << "\n"
                << "            Gene List  (-g) :  " << genefile << "\n"
                << "          Target File  (-t) :  " << targetfile << "\n"
                << "      Background Genes (-b) :  " << ( background_list ? backgroundfile : "--no background set--" ) << "\n"
                << "           Range File  (-x) :  " << ( !use_range ? "--no ranges--" : rangefile ) << "\n";


        ss      << "            Test Type  "
                << ( test == TARGET ? "(-2) :  TARGETS\n" : "(-1) :  INTERVALS\n");
        ss      << "Target Size Filter  (-i,-j) :  ";
        if ( tsize_min ) ss << tsize_min; else ss << 1;
        ss      << "..";
        if ( tsize_max ) ss << tsize_max; else ss << "all";
        ss << "\n";
        ss      << "     Min-obs Threshold (-z) :  " << min_cnt << "\n";

        ss      << "            Kb-Window  (-w) :  " << window << "\n"
                << "            Display-P  (-p) :  " << pthresh << "\n";

        if ( topn == 0 ) ss << "        Top-N Regions  (-n) :   --all--\n";
        else ss << "        Top-N Regions  (-n) :  " << topn << "\n";

        ss      << "       Num Replicates  (-r) :  " << nrep << "\n"
                << "       Num Bootstraps  (-q) :  " << nboot << "\n"
                                << "          Random Seed  (-f) :  " << seed << "\n";

        if ( use_map ) ss << "        Match Density  (-d) :  " << tol << "\n";
        else           ss << "        Match Density  (-d) :  " << "--no map--\n";

        ss      << "           Precompute  (-c) :  " << ( precompute ? "YES" : "NO" ) << "\n"
                << "          Match Genes  (-e) :  " << ( match_genes ? "YES" : "NO" ) << "\n";

        if(	exam_bp_dist ) ss << "Pairwise Distance (-h) :  " << exam_bp_dist_top << "\n";

    }

    ss <<  "----------------------------------------------------------------\n\n";


    return ss.str();
}


std::string InputData::get_aligator_qt_cmd() {

    std::ostringstream ss;

    ss      << "----------------------------------------------------------------\n"
            << "ALIGATOR v1.0 : " + cur_date()
            << "Reference : Peter Holmans et al. 2010. Am J Hum Genet  \n"
            << "http://x004.psycm.uwcm.ac.uk/~peter/\n"
            << "Program  : http://atgu.mgh.harvard.edu/inrich\n"
            << "----------------------------------------------------------------\n"
            << "Project Title  (-o) \t:  " << outroot << "\n"
            << "SNP File  (-a) \t:  " << testfile << "\n"
            << "Map File  (-m) \t:  " << ( use_map ? mapfile : "--no map--" ) << "\n"
            << "Gene File  (-g) \t:  " << genefile << "\n"
            << "Target File  (-t) \t:  " << targetfile << "\n"
            << "Background File  (-b) \t:  " << ( background_list ? backgroundfile : "--no range--" ) << "\n";

    ss      << "Target Size Filter  (-i,-j) \t:  ";
    if ( tsize_min ) ss << tsize_min; else ss << 1;
    ss      << "..";
    if ( tsize_max ) ss << tsize_max; else ss << "all";
    ss << "\n";
    ss      << "Min-obs Threshold (-z) \t:  " << min_cnt << "\n";

    ss      << "Kb-Window  (-w) \t:  " << window << "\n"
            << "Display-P  (-p) \t:  " << pthresh << "\n";


    ss      << "Num Replicates  (-r) \t:  " << nrep << "\n"
            << "Num Bootstraps  (-q) \t:  " << nboot << "\n"
            << "Random Seed  (-f) \t:  " << seed << "\n";

    ss <<  "----------------------------------------------------------------\n\n";


    return ss.str();
}

std::string InputData::get_inrich_qt_cmd() {

    std::ostringstream ss;



    ss      << "----------------------------------------------------------------\n"
            << inrich_ver << " : " + cur_date()
            << "http://atgu.mgh.harvard.edu/inrich\n"
            << "----------------------------------------------------------------\n";

    if(convert_mode) {
        ss      << "Project Title  (-o) \t:  " << outroot << "\n"
                << "Test Regions  (-a) \t:  " << testfile << "\n"
                << "Map File  (-m) \t:  " << ( use_map ? mapfile : "--no map--" ) << "\n"
                << "LD DB File (-s) \t:  " << use_db_name << "\n";
    }
    else {
        ss      << "Project Title  (-o) \t:  " << outroot << "\n"
                << "Test Regions  (-a) \t:  " << testfile << "\n"
                << "Map File  (-m) \t:  " << ( use_map ? mapfile : "--no map--" ) << "\n"
                << "Gene List  (-g) \t:  " << genefile << "\n"
                << "Target File  (-t) \t:  " << targetfile << "\n"
                << "Background Genes (-b) \t:  " << ( background_list ? backgroundfile : "--no background set--" ) << "\n"
                << "Range File  (-x) \t:  " << ( !use_range ? "--no ranges--" : rangefile ) << "\n";


        ss      << "Test Type  "
                << ( test == TARGET ? "(-2) \t:  TARGETS\n" : "(-1) \t:  INTERVALS\n");
        ss      << "Target Size Filter  (-i,-j) \t:  ";
        if ( tsize_min ) ss << tsize_min; else ss << 1;
        ss      << "..";
        if ( tsize_max ) ss << tsize_max; else ss << "all";
        ss << "\n";
        ss      << "Min-obs Threshold (-z) \t:  " << min_cnt << "\n";

        ss      << "Kb-Window  (-w) \t:  " << window << "\n"
                << "Display-P  (-p) \t:  " << pthresh << "\n";

        if ( topn == 0 ) ss << "Top-N Regions  (-n) \t:   --all--\n";
        else ss << "Top-N Regions  (-n) \t:  " << topn << "\n";

        ss      << "Num Replicates  (-r) \t:  " << nrep << "\n"
                << "Num Bootstraps  (-q) \t:  " << nboot << "\n"
                << "Random Seed  (-f) \t:  " << seed << "\n";

        if ( use_map ) ss << "Match Density  (-d) \t:  " << tol << "\n";
        else           ss << "Match Density  (-d) \t:  " << "--no map--\n";

        ss      << "Precompute  (-c) \t:  " << ( precompute ? "YES" : "NO" ) << "\n"
                << "Match Genes  (-e) \t:  " << ( match_genes ? "YES" : "NO" ) << "\n";

        if(	exam_bp_dist ) ss << "Pairwise Distance (-h) \t:  " << exam_bp_dist_top << "\n";

    }

    ss <<  "----------------------------------------------------------------\n\n";


    return ss.str();
}

std::string  InputData::check_file(std::string _file, std::string _filetype, bool _should ) {
    if(_should) {
        if(_file.compare("")==0) {
            return "* Input Required: " + _filetype + "\n\n";
        }
        else if ( ! fileExists( _file ) )
        {
            return "* No " + _filetype + "\n\n" + _file + "\n\n";
        }
    }
    else if(_file.compare("")!=0 && ! fileExists( _file ) ) {
        return "* No " + _filetype + "\n\n" + _file + "\n\n";
    }
    return "";
}

std::string  InputData::check_input_data() {
    std::string msg = "";

    if(outroot.compare("")==0) {
        msg = msg + "* Input Required: Project Title\n\n";
    }    


    if(aligator) {
        msg = msg + check_file(testfile, "Associated Interval File");
        msg = msg + check_file(mapfile, "Reference SNP File");
        msg = msg + check_file(genefile, "Reference Gene File");
        msg = msg + check_file(targetfile, "Target Set File");
    }
    /*
    else if(convert_mode) {
        msg = msg + check_file(testfile, "Associated SNP File");
        msg = msg + check_file(mapfile, "Reference SNP Map File");

        if(create_lddb) {            
            msg = msg + check_file(create_db_name, "LD Information File");
        }
        else if(use_lddb) {
            msg = msg + check_file(use_db_name, "LD DB File");
        }
    }
    */
    else {
        msg = msg + check_file(testfile, "Associated Interval File");
        msg = msg + check_file(genefile, "Reference Gene Map File");
        msg = msg + check_file(targetfile, "Target Set File");

	msg = msg + check_file(mapfile, "Reference SNP Map File", false);
        msg = msg + check_file(rangefile, "Range File", false);
        msg = msg + check_file(backgroundfile, "Background Gene File", false);

        if(use_map && use_range) {
            msg = msg + "ERROR: only one of map file or range file can be used";
        }
        if(!use_map && !use_range) {
            msg = msg + "ERROR: either one of map file or range file should be used";
        }

        if(tol<=0 || tol>=1.0 ) {
            msg = msg + "Note: 0.0 < Random Interval Length Mapping Ratio  < 1.0\n";
        }

        if(tsize_min && tsize_max && tsize_min > tsize_max) {
            msg = msg + "Note: Minimum Target Set Size <= Maximum Target Set Size\n";
        }

        if(exam_bp_dist && (exam_bp_dist_top<0 || exam_bp_dist_top>100)) {
            msg = msg + "Note: 1 <= Interval Clustering Diplay Top N <= 100";
        }
    }

    return msg;
}

void InputData::set_convert_param(
        std::string _outroot,
        std::string _testfile,
        std::string _mapfile,
        bool _create_db,
        std::string _create_dbName,
        std::string _useDBName ) {

    convert_mode = true;

    outroot = _outroot;
    testfile = _testfile;
    use_map = true;
    mapfile = _mapfile;
    //create_lddb = _create_db;
    //use_lddb = !create_lddb;
    //create_db_name = _create_dbName;
    //use_db_name = _useDBName;
}

void InputData::set_default_param(std::string _project_title,
                                  bool _aligator,
                                  std::string _testfile,
                                  std::string _genefile,
                                  std::string _mapfile,
                                  std::string _targetfile ) {

    outroot = _project_title;
    aligator = _aligator;

    testfile =  _testfile;
    genefile = _genefile;
    targetfile = _targetfile;


    mapfile = _mapfile;
    if(mapfile.compare("")!=0) {
        use_map = true;
    }
    else {
        use_map = false;
    }

    nrep = 5000;
    nboot = 1000;
    tsize_min = 2;
}

void InputData::set_advanced_param(std::string _project_title,
                                   bool _aligator,
                                   std::string _testfile,
                                   std::string _genefile,
                                   std::string _mapfile,
                                   std::string _targetfile,
                                   std::string _rangefile,
                                   std::string _backfile,
                                   int _min_pathway_size,
                                   int _max_pathway_size,
                                   int _min_cnt,
                                   int _window,
                                   double _p,
                                   int _nrep,
                                   int _nboot,
                                   int _seed,
                                   bool _intervalMode,
                                   bool _precompute,
                                   bool _match_gene,
                                   double _tol,
                                   int _topn,
                                   bool _exam_bp_dist,
                                   int _exam_bp_dist_topn) {

    set_default_param(_project_title,
                      _aligator,
                      _testfile,
                      _genefile,
                      _mapfile,
                      _targetfile );

    // ---------------------------------
    // Inputfile
    // ---------------------------------
    rangefile = _rangefile;
    if(  rangefile.compare("")!=0 ) {
        use_range = true;
    }    
    else {
        use_range = false;
    }

    backgroundfile = _backfile;
    if(  backgroundfile.compare("")!=0 ) {
        background_list = true;
    }

    // ---------------------------------
    // Parameters
    // ---------------------------------
    tsize_min = _min_pathway_size;
    tsize_max = _max_pathway_size;
    min_cnt = _min_cnt;
    window = _window;
    pthresh = _p;
    nrep = _nrep;
    nboot = _nboot;
    seed = _seed;
    precompute = _precompute;
    match_genes = _match_gene;
    tol = _tol;

    topn = _topn;
    // largest_gene = 5000000;

    if(_intervalMode) {
        test = INTERVAL;
    }
    else {
        test = TARGET;
    }


    exam_bp_dist = _exam_bp_dist;
    exam_bp_dist_top = _exam_bp_dist_topn;

}

void InputData::set_aligator_param(std::string _project_title,
                                   std::string _testfile,
                                   std::string _mapfile,
                                   std::string _genefile,
                                   std::string _targetfile,
                                   int _min_pathway_size,
                                   int _max_pathway_size,
                                   int _min_cnt,
                                   int _window,
                                   int _nrep,
                                   int _nboot,
                                   int _seed) {

    set_default_param(_project_title,
                      true,
                      _testfile,
                      _genefile,
                      _mapfile,
                      _targetfile );


    // ---------------------------------
    // Parameters
    // ---------------------------------
    tsize_min = _min_pathway_size;
    tsize_max = _max_pathway_size;
    min_cnt = _min_cnt;
    window = _window;
    nrep = _nrep;
    nboot = _nboot;
    seed = _seed;
}

/*
bool InputData::included( const Interval & i ) const
{
    std::map<int,Interval>::const_iterator k = scaffold.find( i.chr );
    if ( k == scaffold.end() ) return false;
    return i.bp1 >= k->second.bp1 && i.bp2 <= k->second.bp2;
}
*/

Interval InputData::random() const
{
    long int i = CRandom::rand( (int)chr_list.size()  );
    int chr = chr_list[i];

    //std::cerr << "i=" << i << " chr=" << chr << " scaffold.size=" << scaffold.size() << "\n";

    std::map<int,Interval>::const_iterator j = scaffold.find(chr);
    if ( j != scaffold.end() ) {
	int bp1 = j->second.bp1 + CRandom::rand( j->second.bp2 - j->second.bp1 + 1 ); 

	// only the starting position matters
	return Interval( chr, bp1, bp1 + 1 );
    }

    // this cannot happen
    return Interval(1,0);
}


/* 
// original version: 2011-Feb-28 replaced
Interval InputData::random() const
{
    long int r = CRandom::rand( total_seq );
    std::map<int,Interval>::const_iterator j = scaffold.begin();
    long int c = 0;
    while ( j != scaffold.end() )
    {
        c += j->second.bp2 - j->second.bp1 + 1 ;

	// only the starting position matters
        if ( r < c )
            return Interval( j->second.chr , j->second.bp1 + ( c - r ) );

        ++j;
    }
    return Interval(1,0);
}
*/


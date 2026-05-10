#include <vector>
#include <iostream>
#include <getopt.h>

void printUsage(){
    std::cerr<<"Usage :\n";
    std::cerr<<"     lfr_per_barcode_53 [options] < long_fragment.details.txt\n\n";
    std::cerr<<"Options :\n"
             <<"     -h/--help                      print this usage.\n"
             <<"     -l/--valid_lfr_readpair        min valid lfr read pair threshold.\n"
             <<"                                        default 3.\n"
             <<"     -b/--valid_barcode_readpair    min valid lfr read pair threshold.\n"
             <<"                                        default 5.\n";
}

int main(int argc , char ** argv)
{
    static struct option long_options[] = {
        {"valid_lfr_readpair",  required_argument, NULL, 'l'},
        {"valid_barcode_readpair",  required_argument, NULL, 'b'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "l:b:h";
    int valid_lfr_readpair = 3;
    int valid_barcode_readpair = 5 ;
    while(1){
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'l':
                valid_lfr_readpair = atoi(optarg);
                break;
            case 'b':
                valid_barcode_readpair = atoi(optarg);
                break;
            case 'h':
                printUsage();
                return 0;
            default:
                printUsage();
                return 1 ;
        }
    }
    if( valid_barcode_readpair < 0 || valid_lfr_readpair < 0){
        std::cerr<<"error : valid_barcode_readpair < 0 || valid_lfr_readpair < 0 . exit ..."<<std::endl;
        return 1 ;
    }
    if( valid_barcode_readpair < valid_lfr_readpair ) {
        std::cerr<<"error : valid_barcode_readpair < valid_lfr_readpair . exit ..."<<std::endl;
        return 1 ;
    }
    int lfr_freq[20];
    for( int i = 0 ; i < 20 ; i++ ) { lfr_freq[i] = 0 ;}
    int prev_barcode = -1 ;
    std::string line ;
    std::getline(std::cin,line); // remove header
    std::vector< std::vector<int> > cache ;
    while( ! std::cin.eof() ) {
        int tmp , barcode , rp ;
        std::cin>>barcode>>tmp>>rp>>tmp>>tmp;
        if( barcode != prev_barcode ) {
            prev_barcode = barcode ;
            cache.push_back(std::vector<int>());
        }
        cache[cache.size()-1].push_back(rp);
    }
    for ( int i = 1 ; i<= valid_barcode_readpair ; i++ ){
        std::cerr<<"##############\t################"<<std::endl;
        std::cerr<<"LFR-filter\tread-pair<"<<i<<std::endl;
        // total count
        int b = 0 ;int bm = 0; int l = 0; int ul = 0 ;
        // filter  bacode-read-pair<valid_barcode_read_pair
        int b5 = 0 ;int bm5 = 0; int l5 = 0; int ul5 = 0;
        // filter  bacode-read-pair<valid_lfr_pair
        int bi = 0 ;int bmi = 0; int li = 0; int uli = 0;
        // count   bacode-read-pair=1
        int bib = 0 ; int biv=0 ; int bil=0;
        for( const auto & d : cache ) {
            int br = 0 ; // read-pair in barcode
            int v  = 0 ; // valid barcode or not
            int lc = 0 ; // valid lfr count
            int uc = 0 ; // uncertain lfr count
            if ( d.size() == 1 ) {
                lc=1 ;
                v =1 ;
                br=d[0];
            } else {
                for( int r : d ) {
                    br += r ;
                    if (r >=i) { lc ++; v=1; } else { uc ++ ; }
                }
            }
            l += lc;
            b += v ;
            ul += uc ;
            if ( br >=valid_barcode_readpair ) {
                b5 ++ ;
                bm5 +=v ;
                l5 +=lc;
                ul5+=uc ;
            }
            if ( br >=i ) {
                bi ++ ;
                bmi +=v ;
                li +=lc;
                uli += uc ;
            }
            if( br > 1 )
                bm +=v ;
            if( i == 1 ){
                if( br > 1 &&br < valid_lfr_readpair ){
                    biv += v ;
                    bib += 1;
                    bil +=lc;
                }
            }
            if( i == valid_lfr_readpair && br >=valid_barcode_readpair ) {
                if( lc > 0 && lc < 20 ) lfr_freq[lc-1]++;
                if( lc > 20 ) lfr_freq[19]++;
            }
        }
        std::cerr<<"Totoal-data\t---------"<<std::endl;
        std::cerr<<"total_barcode\t"<<cache.size()<<std::endl;
        std::cerr<<"uncertain_barcode\t"<<cache.size()-b<<std::endl;
        std::cerr<<"valid_barcode\t"<<b<<std::endl;
        std::cerr<<"valid_lfr\t"<<l<<std::endl;
        std::cerr<<"uncertain_lfr\t"<<ul<<std::endl;
        std::cerr<<"valid_lfr/valid_barcode\t"<<float(l)/float(b)<<std::endl;

        std::cerr<<"filter-barcode-read-pair<2\t---------"<<std::endl;
        std::cerr<<"total_barcode\t"<<cache.size()-(b-bm)<<std::endl;
        std::cerr<<"uncertain_barcode\t"<<cache.size()-b<<std::endl;
        std::cerr<<"valid_barcode\t"<<bm<<std::endl;
        std::cerr<<"valid_lfr\t"<<l-(b-bm)<<std::endl;
        std::cerr<<"(valid_lfr/valid_barcode)\t"
                                     <<float(l-(b-bm))/float(bm)<<std::endl;
        if( i >2 && i < valid_barcode_readpair ) {
            std::cerr<<"filter-barcode-read-pair-filete<"<<i<<"\t---------"<<std::endl;
            std::cerr<<"total_barcode\t"<<bi<<std::endl;
            std::cerr<<"uncertain_barcode\t"<<bi-bmi<<std::endl;
            std::cerr<<"valid_barcode\t"<<bmi<<std::endl;
            std::cerr<<"valid_lfr\t"<<li<<std::endl;
            std::cerr<<"uncertain_lfr\t"<<uli<<std::endl;
            std::cerr<<"(valid_lfr/valid_barcode)\t"
                <<float(li)/float(bmi)<<std::endl;
        }
        if ( valid_barcode_readpair > 2 ) {
            std::cerr<<"filter-barcode-read-pair<"<<valid_barcode_readpair<<"\t---------"<<std::endl;
            std::cerr<<"total_barcode\t"<<b5<<std::endl;
            std::cerr<<"uncertain_barcode\t"<<b5-bm5<<std::endl;
            std::cerr<<"valid_barcode\t"<<bm5<<std::endl;
            std::cerr<<"valid_lfr\t"<<l5<<std::endl;
            std::cerr<<"uncertain_lfr\t"<<ul5<<std::endl;
            std::cerr<<"(valid_lfr/valid_barcode)\t"
                <<float(l5)/float(bm5)<<std::endl;
        }
        if( i == 1 ){
            std::cerr<<"total_barcode(if-1<barcode-read-pair<)"<<valid_barcode_readpair<<"\t"<<bib<<std::endl;
            //std::cerr<<"uncertain_barcode(if-1<barcode-read-pair<=5))\t"<<bib-biv<<std::endl;
            std::cerr<<"valid_barcode(if-1<barcode-read-pair<)"<<valid_barcode_readpair<<"\t"<<biv<<std::endl;
            std::cerr<<"valid_lfr(if-1<barcode-read-pair<)"<<valid_barcode_readpair<<"\t"<<bil<<std::endl;
            std::cerr<<"valid_lfr/valid_barcode(if-1<barcode-read-pair<)"<<valid_barcode_readpair<<"\t"<<float(bil)/float(biv)<<std::endl;
        }
    }
    std::cerr<<"\nvalid-lfr-per-barcode-count:\n";
    for( int i = 0 ; i < 20 ; i++ ){
        std::cerr<<i+1<<'\t'<<lfr_freq[i]<<'\n';
    }
    std::cerr<<"notice : >20 count in 20!\n";
}

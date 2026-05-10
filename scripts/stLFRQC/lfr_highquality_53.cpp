#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <getopt.h>

void printUsage(){
    std::cerr<<"Usage :\n";
    std::cerr<<"     lfr_length_53 [options] < long_fragment.details.txt\n\n";
    std::cerr<<"Options :\n"
             <<"     -h/--help                      print this usage.\n"
             <<"     -l/--valid_lfr_readpair        min valid lfr read pair threshold.\n"
             <<"                                        default 3.\n"
             <<"     -b/--valid_barcode_readpair    min valid lfr read pair threshold.\n"
             <<"                                        default 5.\n"
             <<"     -g/--good_lfr_readpair         default 5\n"
             <<"     -m/--good_lfr_length           default 10,000\n";
}

int main(int argc , char ** argv)
{
    static struct option long_options[] = {
        {"valid_lfr_readpair",      required_argument, NULL, 'l'},
        {"valid_barcode_readpair",  required_argument, NULL, 'b'},
        {"good_lfr_readpair"    ,   required_argument, NULL, 'g'},
        {"good_lfr_length"    ,     required_argument, NULL, 'm'},
        {"help",  no_argument,       NULL, 'h'},
        {0, 0, 0, 0}
    };
    static char optstring[] = "l:b:g:m:h";
    int valid_lfr_readpair = 3;
    int valid_barcode_readpair = 5 ;
    int good_lfr_readpair = 5 ;
    int good_lfr_length = 10000;
    while(1) {
        int c = getopt_long(argc, argv, optstring, long_options, NULL);
        if (c<0) break;
        switch (c){
            case 'l':
                valid_lfr_readpair = atoi(optarg);
                break;
            case 'b':
                valid_barcode_readpair = atoi(optarg);
                break;
            case 'g':
                good_lfr_readpair = atoi(optarg);
                break;
            case 'm':
                good_lfr_length = atoi(optarg);
                break;
            case 'h':
                printUsage();
                return 0;
            default:
                printUsage();
                return 1 ;
        }
    }
    if( valid_barcode_readpair < 0 || valid_lfr_readpair < 0) {
        std::cerr<<"error : valid_barcode_readpair < 0 || valid_lfr_readpair < 0 . exit ..."<<std::endl;
        return 1 ;
    }
    if( valid_barcode_readpair < valid_lfr_readpair ) {
        std::cerr<<"error : valid_barcode_readpair < valid_lfr_readpair . exit ..."<<std::endl;
        return 1 ;
    }
    if( good_lfr_readpair < valid_lfr_readpair ) {
        std::cerr<<"error : good_lfr_readpair < valid_lfr_readpair . exit ..."<<std::endl;
        return 1 ;
    }
    // cache all lfr info
    std::vector<std::vector<std::pair<int,int>> > cache ;

    // 1. load long_fragment_detail.txt
    int prev_barcode = -1 ;
    std::string line ;
    std::getline(std::cin,line); // remove header
    while( ! std::cin.eof() ) {
        int tmp , barcode , rp , length;
        std::cin>>barcode>>tmp>>rp>>tmp>>length;
        if(barcode != prev_barcode ) {
            prev_barcode = barcode ;
            cache.push_back(std::vector<std::pair<int,int>>());
        }
        assert(length>0);
        cache[cache.size()-1].push_back(std::make_pair(rp,length));
    }

    // 2. filter valid lfr information
    int min_brp = valid_barcode_readpair ; // min read-pair-number of valid barcode
    int min_lrp = valid_lfr_readpair; // min read-pair-number of valid lfr
    std::cerr<<"##############\t################"<<std::endl;
    std::cerr<<"valid-LFR-filter\tread-pair<"<<min_lrp<<std::endl;
    std::cerr<<"valid-barcode-filter\tread-pair<"<<min_brp<<std::endl;
    std::cerr<<"good-lfr-filter\tread-pair<"<<good_lfr_readpair<<std::endl;
    std::cerr<<"good-lfr-filter\tlength<"<<good_lfr_length<<std::endl;
    int total_barcode = 0 ;
    int total_lfr = 0 ;
    int valid_barcode = 0 ;
    int valid_lfr = 0 ;
    int good_lfr = 0 ;
    for(const auto & d : cache ) {
        total_barcode ++ ;
        int br = 0 ; // read-pair in barcode
        int v = 0  ; // valid barcode or not
        int lc = 0 ;
        int ulc = 0 ;
        if( d.size() == 1 ) {
            v =1;
            br=d.at(0).first;
            lc = 1 ;
        } else {
            for( const std::pair<int,int> & p : d ) {
                br += p.first ;
                if ( p.first >=min_lrp )
                {    v=1; lc ++ ; } else { ulc ++ ; }
            }
        }
        total_lfr += ( lc + ulc );
        valid_lfr += lc ;
        valid_barcode += v ;
        if( v==1 && br >=min_brp ) {
            if( d.size() == 1 && good_lfr_readpair == 1 && d.at(0).second >= good_lfr_length )
                good_lfr++;
            else {
                for( const std::pair<int,int> & p : d )
                    if ( p.first >=good_lfr_readpair && p.second >=good_lfr_length )
                        good_lfr ++ ;
            }
        }
    }
    std::cerr<<"--------------\t----------------"<<std::endl;
    std::cerr<<"total-barcode\t"<<total_barcode<<'\n';
    std::cerr<<"valid-barcode\t"<<valid_barcode<<'\n';
    std::cerr<<"total-lfr\t"<<total_lfr<<'\n';
    std::cerr<<"valid-lfr\t"<<valid_lfr<<'\n';
    std::cerr<<"good-lfr\t"<<good_lfr<<'\n';
    std::cerr<<"good-lfr/valid-lfr\t"<<float(good_lfr)/float(valid_lfr)<<'\n';
    std::cerr<<"##############\t################"<<std::endl;
}

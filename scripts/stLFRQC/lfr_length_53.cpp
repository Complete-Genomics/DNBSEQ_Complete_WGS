#include <vector>
#include <iostream>
#include <map>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <getopt.h>

// lfr-length-distrubution , bined by 1000 bp
std::map<int , int > lfr_length_freq_1000;
// lfr-length-cache : for min ,max and median
std::vector<long long> length_cache;

void incr(int length) {
    assert(length > 0 );
    int b_index=length/1000;
    if( lfr_length_freq_1000.find(b_index) == lfr_length_freq_1000.end() )
        lfr_length_freq_1000[b_index]=1;
    else
        lfr_length_freq_1000[b_index]++;
    length_cache.push_back(length);
}

void printUsage(){
    std::cerr<<"Usage :\n";
    std::cerr<<"     lfr_length_53 [options] < long_fragment.details.txt\n\n";
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
    std::cerr<<"Barcode-filter\tread-pair<"<<min_brp<<std::endl;
    std::cerr<<"LFR-filter\tread-pair<"<<min_lrp<<std::endl;
    for(const auto & d : cache ) {
        int br = 0 ; // read-pair in barcode
        int v = 0  ; // valid barcode or not
        if( d.size() == 1 ) {
            v =1;
            br=d.at(0).first;
        }else {
            for( const std::pair<int,int> & p : d ) {
                br += p.first ;
                if ( p.first >=min_lrp )
                    v=1; 
            }
        }
        if( v==1 && br >=min_brp ) {
            if( d.size() == 1 ) 
                incr(d.at(0).second);
            else {
                for( const std::pair<int,int> & p : d )
                    if ( p.first >=min_lrp ) 
                        incr(p.second);
            }
        }
    }
    std::sort(length_cache.begin(),length_cache.end());
    long long sum = 0 ;
    for( long x : length_cache) sum +=x ;
    std::cerr<<"Average-lfr-length\t"<<sum/length_cache.size()<<std::endl;
    std::cerr<<"Min-lfr-length\t"<<length_cache.at(0)<<std::endl;
    std::cerr<<"Max-lfr-length\t"<<length_cache.at(length_cache.size()-1)<<std::endl;
    std::cerr<<"Median-lfr-length\t"<<length_cache.at(length_cache.size()/2)<<std::endl;
    std::cerr<<"LFR-length/1000\tcount\n";
    for(const std::pair<long,long> & p : lfr_length_freq_1000 ){
        std::cerr<<p.first<<"\t"<<p.second<<"\n";
    }
}

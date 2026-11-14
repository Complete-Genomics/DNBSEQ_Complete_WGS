BEGIN{
    BS["0_0_0"]=0;
    bindex=1;
    RB[""]=0;
    rindex=1;
    RF[""]=0;
    findex=1;
}
{
    if( $1 in RN ){
        ri = RN[$1]
    }else {
        ri = rindex ;
        RN[$1] = rindex;
        rindex = rindex + 1;
    }
    if(length(RN) > 300000 ) {
        delete RN;
        RN[""]=0;
    }
    if( $2 in BS ){
        bs = BS[$2];
    }else {
        bs = bindex ;
        BS[$2] = bindex;
        bindex = bindex + 1;
    }
    if( $4 in RF ){
        rf = RF[$4];
    }else {
        rf = findex ;
        RF[$4] = findex;
        findex = findex + 1;
    }
    if( $2 != "0_0_0" && ((0+$3) < 1024) && $4 != "*" )
        printf("%s\t%s\t%s\t%s\n",ri,bs,rf,$5);
}

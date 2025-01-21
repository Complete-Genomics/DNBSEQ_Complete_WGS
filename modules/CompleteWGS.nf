nextflow.enable.dsl=2

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.gatk_version != "v3" && params.gatk_version != "v4") { exit 1, 'wrong gatk version!'}
if (!params.use_megabolt && params.dv_version == "v0.7" ) { exit 1, 'dv v0.7 can only be used in megabolt!'}

include {toCsv                                     } from "${params.MOD}/parseSample"                        
include {
    parse_sample;
    tosamplelist;
    parse_sample_frombam;
    bam;
    bam2                              } from "${params.MOD}/parseSample"
include {
    fq;
    fq1; 
    readLen;
    readLen as readLenPf;
    basecount;
    basecount as basecountPf;
    splitfq;
    samplePfFq;
    sampleStlfrFq;
    fqcheck;
    fqcheck as fqcheckPf;
    fqdist;
    fqdist as fqdistPf;
    eachstat_fastq } from "${params.MOD}/fq"
include {
    qc as qc_pf;
    qc as qc_stlfr;
    qc_stlfr_stats                                 } from "${params.MOD}/qc"
include {fqstats as fqstats_pf                     } from "${params.MOD}/fqstats"
include {fqstats as fqstats_stlfr                  } from "${params.MOD}/fqstats"
include {
    barcode_split;
    barcode_split_stLFRreseq;
    makeLog                                        } from "${params.MOD}/barcodeSplit"
include {
    tofake10xHash;
    tofake10x;
    mergeMaps;
    bwaMegabolt;    //megabolt module: alignmentsortmarkdupbqsr
    bwaMegabolt as bwaMegaboltPf;
    bqsr as bqsrPf;
    bqsr as bqsrStlfrLariat;
    bqsr as bqsrStlfrBwa;
    bqsrMegabolt;
    bwa;
    bwa as bwaPf;
    lariatBC;
    fake10x2lariat;
    mergeFq;
    lariat;
    sortbam;
    markdup as markdupPf;
    markdup as markdupStlfrBwa;
    markdup as markdupStlfrLariat;
    intersect as intersectLariat;
    intersect as intersectBwa;
    mergeBam as mergeBamLariat;
    mergeBam as mergeBamBwa;
    sampleBam as sampleBamPf;
    sampleBam as sampleBamStlfrBwa;
    sampleBam as sampleBamStlfrLariat;
    stLFRQC} from "${params.MOD}/align"
include { mapq;
    mapq as mapq_frombam            } from "${params.MOD}/bam"
include {
    frag1;
    frag2                           } from "${params.MOD}/frag"
include {
    coverage;
    coverage as coveragePf;
    coverageMean;
    coverageMean as coverageMeanPf;
    coverageAvg;
    mosdepth;
    mosdepth as mosdepthPf;
    mosdepthPlot                                  } from "${params.MOD}/genedepth"
include {
    gatk_interval;
    gatherVcfsHc;
    gatherVcfsHc as gatherVcfsHcBwa;
    gatherVcfsVqsr;
    gatherVcfsVqsr as gatherVcfsVqsrBwa;
    deepvariantv16;
    deepvariantv16 as deepvariantv16Bwa;
    deepvariantv16 as deepvariantv16BwaPf;
    dvMegabolt;
    dvMegabolt as dvMegaboltBwa;
    hcMegabolt;
    hcMegabolt as hcMegaboltBwa;
    hcMegabolt as hcMegaboltBwaPf;
    vqsrMegabolt;
    vqsrMegabolt as vqsrMegaboltBwa;
    vqsrMegabolt as vqsrMegaboltBwaPf;
    hc;
    hc as hcBwa;
    hcSplit;
    hcSplit as hcSplitBwa;
    vqsrSnp;
    vqsrSnp as vqsrSnpBwa;
    vqsrIndel;
    vqsrIndel as vqsrIndelBwa } from "${params.MOD}/callvariants"

include {
    splitBam4phasing;
    splitBam4phasing as splitBam4phasingBwa                          } from "${params.MOD}/splitbam"
include {
    vcfeval as vcfevalBwaDv;
    vcfeval as vcfevalBwaGatk;
    vcfeval as vcfevalLariatDv;
    vcfeval as vcfevalLariatGatk;
    vcfeval as vcfevalPf;
    eachstat_vcf;
    variant_fix } from "${params.MOD}/vcfeval"

include {
    variant_stats as varStatsBwaDv;
    variant_stats as varStatsBwaGatk;
    variant_stats as varStatsLariatDv;
    variant_stats as varStatsLariatGatk } from "${params.MOD}/variantstats"

include {
    split_vcf as splitVcfBwaDv;
    split_vcf as splitVcfBwaGatk;
    split_vcf as splitVcfLariatDv;
    split_vcf as splitVcfLariatGatk } from "${params.MOD}/splitvcf"

include {
    getchrs;
    phase as phaseBwaDv;
    phase as phaseBwaGatk;
    phase as phaseLariatDv;
    phase as phaseLariatGatk;
    phase_cat as phaseCatBwaDv;
    phase_cat as phaseCatBwaGatk;
    phase_cat as phaseCatLariatDv;
    phaseCat_cwx;
    phase_cat as phaseCatLariatGatk;
    phaseCatRef;
    eachstat_phase;
    hapKaryotype;
    hapKaryotype_bak;
    circos } from "${params.MOD}/phase"
include {
    cnv;
    smoove;
    svpre;
    lfrsv;
    mergesv } from "${params.MOD}/cnvsv"
include {
    bamdepth;
    bamdepth as bamdepthPf;
    barcodeStat;
    sameChrBCratio;
    samtools_flagstat;
    samtools_flagstat as samtoolsFlagstatPf;
    samtools_stats;
    samtools_stats as samtoolsStatsPf;
    insertsize;
    insertsize as insertsizePf;
    gcbias;
    eachstat_cov;
    eachstat_depth;
    samtools_depth;
    samtools_depth as samtoolsDepthPf;
    samtools_depth0;
    samtools_depth0 as samtoolsDepth0Pf;
    splitDepth as splitDepthPf;
    splitDepth as splitDepthStlfr;
    genomeDepth as genomeDepthPf;
    genomeDepth as genomeDepthStlfr;
    align_cat;
    eachstat_aligncat;
    align_cat as alignCatPf;
    align_catAll } from "${params.MOD}/bamstats"

include {
    report0 as reportBwaGatk;
    report0 as reportBwaDv;
    report0 as reportLariatGatk;
    report0 as reportLariatDv;
    report01 as reportBwaGatk1;
    report01 as reportBwaDv1;
    report01 as reportLariatGatk1;
    report01 as reportLariatDv1;
    report;
    html } from "${params.MOD}/report"


workflow CWGS {
    def ch_libpf = "pf"
    def ch_libstlfr = "stlfr"
    def ch_merge = "merge"

    def ch_bwa = "bwa"
    def ch_lariat = "lariat"

    def ch_dv = "dv"
    def ch_gatk = "gatk"
    def ch_gatk3 = "gatk3"
    def ch_gatk4 = "gatk4"

    ch_bam1 = Channel.empty()
    
    //// stLFRreseq v1.4
    if (params.stLFR_only) {
        // input: sampleid, fqpath, barcode, ref
        parse_sample(tosamplelist(ch_input)).set {ch_fq}
        fq1(ch_fq).set {ch_stlfrfq}

        // split barcode & QC & fqstats
        barcode_split_stLFRreseq(ch_stlfrfq).reads.set{ch_splitfq}
        splitLog = barcode_split_stLFRreseq.out.log 
        qc_stlfr(ch_libstlfr, ch_splitfq).reads.set{ch_stlfrfq}
        qc_stlfr.out.bssq.set{ch_stlfrbssq}
        fqstats_stlfr(ch_stlfrbssq)

        fqcheck(ch_libstlfr, ch_stlfrfq).set {ch_stlfrfqcheck}
        fqdist(ch_libstlfr, ch_stlfrfqcheck) //base, comp png
        

        // mapping (bwa only) & bamstats
        def chrs = (1..22).collect { "chr$it" } + ["chrX", "chrY"]
        if (params.use_megabolt) {
            bwaMegabolt(ch_libstlfr, ch_splitfq).set {ch_stlfrbam}
        } else {
            bwa(ch_libstlfr, ch_splitfq).set {ch_stlfrBwaBam}
            markdupStlfrBwa(ch_libstlfr, ch_bwa, ch_stlfrBwaBam).set {ch_stlfrbam}

            if (params.var_tool.contains("gatk")) { // not use megabolt
                bqsrStlfrBwa(ch_libstlfr, ch_bwa, ch_stlfrbam).set {ch_stlfrbam}
            }
        }
        samtools_flagstat(ch_libstlfr, ch_stlfrbam).set {ch_flagstat}
        samtools_stats(ch_libstlfr, ch_stlfrbam).set {ch_stat}
        insertsize(ch_libstlfr, ch_stlfrbam).insertsize.set {ch_insertsize} 
        samtools_depth(ch_libstlfr, ch_stlfrbam).set {ch_depthreport}
        align_cat(ch_libstlfr, ch_flagstat.join(ch_stat).join(ch_depthreport).join(ch_insertsize)).set {ch_aligncatstlfr} //info
        eachstat_cov(ch_stlfrbam).info1.set {info1}
        eachstat_depth(ch_stlfrbam).info2.set {info2}        //id.Sequencing.depth.pdf, id.Sequencing.depth.accumulation.pdf, id.sorted.bam.info_2.xls
        gcbias(ch_libstlfr, ch_stlfrbam)

        // split bam & LFR stats & stats
        splitBam4phasingBwa(ch_bwa, ch_stlfrbam, chrs).set {ch_eachbambwa}
        frag1(ch_eachbambwa).groupTuple().set {ch_frag1}
        frag2(ch_frag1) // id.frag_*.pdf 
        eachstat_aligncat(ch_flagstat.join(ch_stat).join(info1).join(info2).join(ch_insertsize))
        eachstat_fastq(ch_stlfrbssq.join(splitLog).join(frag2.out.bcinfo))

        // call variants
        if (params.var_tool == "dv") {
            if (params.use_megabolt) {dvMegaboltBwa(ch_bwa, ch_stlfrbam).set {ch_vcf}}
            else if (params.dv_version == "v1.6") {deepvariantv16Bwa(ch_bwa, ch_stlfrbam).set {ch_vcf}}
        } else if (params.var_tool == "gatk") {
            if (params.use_megabolt) {
                vqsrMegaboltBwa(ch_bwa, hcMegaboltBwa(ch_bwa, ch_stlfrbam)).set {ch_vcf}
            } else {
                hc2(ch_bwa, ch_stlfrbam).set {ch_mergevcf4}
                vqsrSnpBwa(ch_bwa, hc(ch_bwa, ch_stlfrbam)).set {ch_vqsrsnp}
                vqsrIndelBwa(ch_bwa, hc(ch_bwa, ch_stlfrbam)).set {ch_vqsrindel}
                gatherVcfsVqsrBwa(ch_bwa, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_vcf}
            }
        }

        // var stats
        vcfevalBwaDv(ch_merge, ch_vcf).set {ch_vcfevalBwaDv}
        eachstat_vcf(ch_vcf).set {ch_vartable}
        varStatsBwaDv(ch_vcf)

        // split vcf for phasing & phasing & phsing stats
        def ch_var = (params.var_tool == "dv") ? ch_dv : ch_gatk
        splitVcfBwaDv(ch_bwa, ch_var, ch_vcf, chrs).eachvcf.set {ch_eachvcf3}
        vcfs3 = splitVcfBwaDv.out.vcf.groupTuple()
        ch_eachbambwa.combine(ch_eachvcf3, by: [0,1]).set {ch_eachchr3}
        pvcfs3 = phaseBwaDv(ch_bwa, ch_var, ch_eachchr3).phasedvcf.groupTuple()  
        lfs3 = phaseBwaDv.out.lf.groupTuple()  
        hbs3 = phaseBwaDv.out.hapblock.groupTuple()  
        stats3 = phaseBwaDv.out.stat.groupTuple()  
        phaseCatBwaDv(ch_bwa, ch_var, vcfs3.join(pvcfs3).join(lfs3).join(hbs3).join(stats3)).report.set {ch_phasereport3}//report
        eachstat_phase(phaseCatBwaDv.out.hapcutstat.join(pvcfs3))

        pvcfs3.join(lfs3).join(hbs3).map { items ->
            def id = items[0]
            def paths = items[1..-1]
            return [id, paths.flatten()]
        } set {ch_phaseallBwaDv}

        // cnvsv
        svpre(phaseBwaDv.out.svpre).groupTuple().map { items ->
            def id = items[0]
            def paths = items[1].flatten()
            return [id, paths]
        }.set {ch_svsplit}   

        //cnvsv
        smoove(ch_stlfrbam).set {ch_smoovevcf}
        ch_stLFRreadLen = barcode_split_stLFRreseq.out.rlen
        lfrsv(ch_stLFRreadLen.join(ch_stlfrbam).join(ch_svsplit)).vcf.set {ch_lfrsvvcf}
        lfrsv.out.bamfinal.set {ch_bamfinal}
        mergesv(ch_smoovevcf.join(ch_lfrsvvcf)).set {ch_sv}
        // stlfrbam.join(ch_vcf).join(phaseall).view()
        cnv(ch_stlfrbam.join(ch_vcf).join(ch_phaseallBwaDv)).set {ch_cnv}  

        hapKaryotype(phaseCatBwaDv.out.hb).pdf.set {ch_haploplot}            // ${id}.haplotype.pdf
        variant_fix(ch_vartable.join(ch_sv).join(ch_cnv))   //${id}.varianttable.xls
        circos(ch_vcf.join(ch_sv).join(ch_cnv).join(ch_bamfinal)).flg.set {ch_circos} // id.circos.png
        html(ch_circos.join(ch_haploplot)) //final html report

    }


    if (params.ref == 'hs37d5') {
        chrs = (1..22).collect { it.toString() } + ['X', 'Y']
    } else if (params.ref == 'hg19' || params.ref == 'hg38') {
        chrs = (1..22).collect { "chr$it" } + ["chrX", "chrY"]
    } else {
        getchrs().set { txt }
        chrs = txt.splitText().map { it.trim() }.collect()
    }

    //// CWGS from fastq
    if (!params.frombam) {
        println("!!! run CWGS from fastq")
        parse_sample (ch_input)
        .reads
        .set { ch_fq }

        fq(ch_fq)
        .stlfr
        .set {ch_stlfrfq}
        
        ch_pffq = fq.out.pf
        qc_pf(ch_libpf, ch_pffq).reads.set {ch_qcpffq} 
        readLenPf(qc_pf.out.bssq).set {ch_PFreadLen} 

        // pffq qc
        qc_pf.out.bssq.set {ch_pfbssq}
        fqcheckPf (ch_libpf, ch_qcpffq).set {ch_pffqcheck} 
        // .groupTuple(sort: sort_filenames('/'))
        // .set {ch_pffqcheck}

        fqdistPf (ch_libpf, ch_pffqcheck) //report 23 25
        fqstats_pf(ch_pfbssq) //report 21

        barcode_split(ch_stlfrfq).reads.set {ch_splitfq}
        splitLog = barcode_split.out.log
        
        if (params.sampleFq) { 
            qc_stlfr_stats(ch_splitfq).basecnt.set {ch_stlfrbasecount}
            qc_stlfr_stats.out.rlen.set {ch_stLFRreadLen}
            sampleStlfrFq(ch_stlfrbasecount.join(ch_stLFRreadLen).join(ch_splitfq).map {id, base, rlen, r1, r2 ->
                tuple(id, base, rlen, [r1, r2])}.transpose())
                .collect()
                .map {id, read1, idd, read2 -> 
                    def r1 = read1.toString().contains("_1.") ? read1 : read2
                    def r2 = read1.toString().contains("_1.") ? read2 : read1
                    tuple(id, r1, r2)}
                // .view()
                .set {ch_stlfrsampledfq}
            basecountPf(ch_pfbssq).set {ch_pfbasecount}
            samplePfFq(ch_pfbasecount.join(ch_PFreadLen).join(ch_qcpffq).map {id, base, rlen, r1, r2 ->
            tuple(id, base, rlen, [r1, r2])}.transpose())
                .collect()
                .map {id, read1, idd, read2 -> 
                    def r1 = read1.toString().contains("_1.") ? read1 : read2
                    def r2 = read1.toString().contains("_1.") ? read2 : read1
                    tuple(id, r1, r2)}
                // .view()
                .set {ch_pfsampledfq}
        } else {
            ch_stlfrsampledfq = ch_splitfq
            ch_pfsampledfq = ch_qcpffq
        }
        
        //// pf align
        if (params.use_megabolt) {
            bwaMegaboltPf(ch_libpf, ch_pfsampledfq).set {ch_pfbam}
        } else {
            bwaPf(ch_libpf, ch_pfsampledfq).set {ch_pfsortbam}
            markdupPf(ch_libpf, ch_bwa, ch_pfsortbam).set {ch_pfbam}
            if (params.var_tool.contains("gatk") && params.run_bqsr) { bqsrPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} }
        } 
        if (params.sampleBam) { sampleBamPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} }
        mapq(ch_pfbam).set {ch_pfbam}

        //pf bam call variant
        deepvariantv16BwaPf(ch_bwa, ch_pfbam).set {ch_pfdvvcf} 
        
        if (!params.ref.startsWith('/')) {
            vcfevalPf(ch_libpf, ch_pfdvvcf).set {ch_vcfevalPf}
            coveragePf(ch_libpf, ch_pfbam).join(coverageMeanPf(ch_libpf, ch_pfbam)).set { ch_PfGeneCov }
        }

        //pf bam stats
        samtoolsFlagstatPf(ch_libpf, ch_pfbam).set {ch_flagstat2}
        samtoolsStatsPf(ch_libpf, ch_pfbam).set {ch_stat2}
        samtoolsDepthPf(ch_libpf, ch_pfbam).set {ch_depthreport2} 

        insertsizePf(ch_libpf, ch_pfbam).insertsize.set {ch_insertsize2} 
        alignCatPf(ch_libpf, ch_flagstat2.join(ch_stat2).join(ch_depthreport2).join(ch_insertsize2)).set {ch_aligncatpf} //info
        bamdepthPf(ch_libpf, ch_pfbam).set {ch_pfbamdepth}

        ////stlfr lariat
        if (params.align_tool.contains("lariat")) {
            if (params.lariatStLFRBC) { //use stLFR bc for lariat
                splitfq(ch_stlfrsampledfq).fq1s.transpose().set {ch_fq1s} //[id, num, fq1, fq2]
                splitfq.out.fq2s.transpose().set {ch_fq2s}
                ch_fq1s.join(ch_fq2s).set {ch_splitstlfrfq}
                lariatBC(ch_splitstlfrfq).groupTuple().set {ch_splitlariatfqs}
                mergeFq(ch_splitlariatfqs).set {ch_lariatfq}
                // ch_splitlariatfqs.map {it -> it[1]}.collect().view()
                // cleanup(ch_splitlariatfqs.map {it -> it[1]}.collect())
            } else {
                tofake10xHash(splitLog).set {ch_hash}
                if (params.lariatSplitFqNum != 1) { 
                    splitfq(ch_stlfrsampledfq).fq1s.transpose().set {ch_fq1s} //[id, num, fq1, fq2]
                    splitfq.out.fq2s.transpose().set {ch_fq2s}
                    ch_fq1s.join(ch_fq2s).set {ch_splitstlfrfq}
                    // ch_splitstlfrfq.view()
                    tofake10x(ch_splitstlfrfq.combine(ch_hash, by:0)).reads.set {ch_splitfake10xfq}
                    fake10x2lariat(ch_splitfake10xfq).groupTuple().set {ch_splitlariatfqs}
                    mergeFq(ch_splitlariatfqs).set {ch_lariatfq}
                } else { fake10x2lariat(tofake10x(ch_stlfrsampledfq.join(ch_hash))).set {ch_lariatfq} }
            }
            
            sortbam(lariat(ch_lariatfq)).set {ch_lariatbam0} 
            markdupStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam0).set {ch_lariatbam}

            if (params.sampleBam) { sampleBamStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam).set {ch_lariatbam} } 
            // sameChrBCratio(barcodeStat(ch_lariatbam)).set {ch_rds}

            if (params.var_tool.contains("gatk") && params.run_bqsr && !params.use_megabolt) { // not use megabolt
                bqsrStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam).set {ch_lariatbam}
            } else if (params.var_tool.contains("gatk") && params.run_bqsr && params.use_megabolt) { // use megabolt
                bqsrMegabolt(ch_lariatbam).set {ch_lariatbam}
            }
            
            //split stLFR bam for phasing
            splitBam4phasing(ch_lariat, ch_lariatbam, chrs).set {ch_eachbamlariat}

            //merge bam
            intersectLariat(ch_lariat, ch_lariatbam.join(ch_pfbam)).set {ch_bed}
            mergeBamLariat(ch_lariat, ch_pfbam.join(ch_lariatbam).join(ch_bed)).set {ch_mergeLariatBam} 
            
            // call variant
            if (params.var_tool.contains("dv")) {
                if (params.use_megabolt && params.dv_version != "v1.6") {dvMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
                else if (params.dv_version == "v1.6") {deepvariantv16(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
                
                if (!params.ref.startsWith('/')) {
                    vcfevalLariatDv(ch_merge, ch_mergevcf).set {ch_vcfevalLariatDv}
                    varStatsLariatDv(ch_mergevcf) 
                }
                

                //phase
                splitVcfLariatDv(ch_lariat, ch_dv, ch_mergevcf, chrs).eachvcf.set {ch_eachvcf}
                vcfs = splitVcfLariatDv.out.vcf.groupTuple()
                ch_eachbamlariat.combine(ch_eachvcf, by: [0,1]).set {ch_eachchr}
                pvcfs = phaseLariatDv(ch_lariat, ch_dv, ch_eachchr).phasedvcf.groupTuple()  
                lfs = phaseLariatDv.out.lf.groupTuple()  
                hbs = phaseLariatDv.out.hapblock.groupTuple()  
                stats = phaseLariatDv.out.stat.groupTuple()  
                if (!params.ref.startsWith('/')) {
                    phaseCatLariatDv(ch_lariat, ch_dv, vcfs.join(pvcfs).join(lfs).join(hbs).join(stats)).report.set {ch_phasereport}//report
                } else {
                    phaseCatRef(ch_lariat, ch_dv, txt, vcfs.join(pvcfs).join(lfs).join(hbs)) //.report.set {ch_phasereport}
                }
                
                // phaseCat_cwx(ch_lariat, ch_dv, vcfs.join(pvcfs).join(lfs).join(hbs).join(stats))

                pvcfs.join(lfs).join(hbs).map { items ->
                    def id = items[0]
                    def paths = items[1..-1]
                    return [id, paths]
                } set {ch_phaseallLariatDv}

            }
            if (params.var_tool.contains("gatk")) {
                if (params.use_megabolt) {
                    if (params.run_vqsr) {
                        vqsrMegabolt(ch_lariat, hcMegabolt(ch_lariat, ch_mergeLariatBam)).set {ch_mergevcf2}
                    } else {
                        hcMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf2}
                    }
                } else {
                    if (params.split_by_intervals) {
                        hcSplit(ch_lariat, ch_mergeLariatBam, intervals).set {ch_mergevcfSplit}
                        gatherVcfsHc(ch_lariat, ch_mergevcfSplit.groupTuple()).set {ch_mergevcf2}
                        if (params.run_vqsr) {
                            vqsrSnp(ch_lariat, ch_mergevcf2).set {ch_vqsrsnp}
                            vqsrIndel(ch_lariat, ch_mergevcf2).set {ch_vqsrindel}
                            gatherVcfsVqsr(ch_lariat, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf2}
                        }
                        
                    } else {
                        hc(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf2}
                        if (params.run_vqsr) {
                            vqsrSnp(ch_lariat, hc(ch_lariat, ch_mergeLariatBam)).set {ch_vqsrsnp}
                            vqsrIndel(ch_lariat, hc(ch_lariat, ch_mergeLariatBam)).set {ch_vqsrindel}
                            gatherVcfsVqsr(ch_lariat, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf2}
                        } 
                    }
                }
                vcfevalLariatGatk(ch_merge, ch_mergevcf2).set {ch_vcfevalLariatGatk}
                varStatsLariatGatk(ch_mergevcf2)

                //phase
                splitVcfLariatGatk(ch_lariat, ch_gatk, ch_mergevcf2, chrs).eachvcf.set {ch_eachvcf2}
                vcfs2 = splitVcfLariatGatk.out.vcf.groupTuple()
                ch_eachbamlariat.combine(ch_eachvcf2, by: [0,1]).set {ch_eachchr2}
                pvcfs2 = phaseLariatGatk(ch_lariat, ch_gatk, ch_eachchr2).phasedvcf.groupTuple()  
                lfs2 = phaseLariatGatk.out.lf.groupTuple()  
                hbs2 = phaseLariatGatk.out.hapblock.groupTuple()  
                stats2 = phaseLariatGatk.out.stat.groupTuple()  
                phaseCatLariatGatk(ch_lariat, ch_gatk, vcfs2.join(pvcfs2).join(lfs2).join(hbs2).join(stats2)).report.set {ch_phasereport2}//report
                pvcfs2.join(lfs2).join(hbs2).set {ch_phaseallLariatGatk}
            }
        }
        ////stlfr bwa
        if (params.align_tool.contains("bwa")) {
            if (params.use_megabolt) {
                bwaMegabolt(ch_libstlfr, ch_splitfq).set {ch_stlfrBwaBam}
            } else {
                bwa(ch_libstlfr, ch_splitfq).set {ch_stlfrBwaBam}
                markdupStlfrBwa(ch_libstlfr, ch_bwa, ch_stlfrBwaBam).set {ch_stlfrBwaBam}
                if (params.var_tool.contains("gatk") && params.run_bqsr) { // not use megabolt
                    bqsrStlfrBwa(ch_libstlfr, ch_bwa, ch_stlfrBwaBam).set {ch_stlfrBwaBam}
                }
            }
            if (params.sampleBam) { sampleBamStlfrBwa(ch_libstlfr, ch_bwa, ch_stlfrBwaBam).set {ch_stlfrBwaBam} }
            splitBam4phasingBwa(ch_bwa, ch_stlfrBwaBam, chrs).set {ch_eachbambwa}

            //merge bam
            intersectBwa(ch_bwa, ch_stlfrBwaBam.join(ch_pfbam)).set {ch_bed2}
            mergeBamBwa(ch_bwa, ch_pfbam.join(ch_stlfrBwaBam).join(ch_bed2)).set {ch_mergebwabam}


            if (params.var_tool.contains("dv")) {
                if (params.use_megabolt) {dvMegaboltBwa(ch_bwa, ch_mergebwabam).set {ch_mergevcf3}}
                else if (params.dv_version == "v1.6") {deepvariantv16Bwa(ch_bwa, ch_mergebwabam).set {ch_mergevcf3}}

                vcfevalBwaDv(ch_merge, ch_mergevcf3).set {ch_vcfevalBwaDv}
                varStatsBwaDv(ch_mergevcf3)

                splitVcfBwaDv(ch_bwa, ch_dv, ch_mergevcf3, chrs).eachvcf.set {ch_eachvcf3}
                vcfs3 = splitVcfBwaDv.out.vcf.groupTuple()
                ch_eachbambwa.combine(ch_eachvcf3, by: [0,1]).set {ch_eachchr3}
                pvcfs3 = phaseBwaDv(ch_bwa, ch_dv, ch_eachchr3).phasedvcf.groupTuple()  
                lfs3 = phaseBwaDv.out.lf.groupTuple()  
                hbs3 = phaseBwaDv.out.hapblock.groupTuple()  
                stats3 = phaseBwaDv.out.stat.groupTuple()  
                phaseCatBwaDv(ch_bwa, ch_dv, vcfs3.join(pvcfs3).join(lfs3).join(hbs3).join(stats3)).report.set {ch_phasereport3}//report
                pvcfs3.join(lfs3).join(hbs3).map { items ->
                    def id = items[0]
                    def paths = items[1..-1]
                    return [id, paths.flatten()]
                } set {ch_phaseallBwaDv}
            }
            if (params.var_tool.contains("gatk")) {
                if (params.use_megabolt) {
                    if (params.run_vqsr) {
                        vqsrMegaboltBwa(ch_bwa, hcMegaboltBwa(ch_bwa, ch_mergebwabam)).set {ch_mergevcf4}
                    } else {
                        hcMegaboltBwa(ch_bwa, ch_mergebwabam).set {ch_mergevcf4}
                    }
                } else {
                    if (params.split_by_intervals) {
                        hcSplitBwa(ch_bwa, ch_mergebwabam, intervals).set {ch_mergevcfSplit}
                        gatherVcfsHcBwa(ch_bwa, ch_mergevcfSplit.groupTuple()).set {ch_mergevcf4}
                        if (params.run_vqsr) {
                            vqsrSnpBwa(ch_bwa, ch_mergevcf2).set {ch_vqsrsnp}
                            vqsrIndelBwa(ch_bwa, ch_mergevcf2).set {ch_vqsrindel}
                            gatherVcfsVqsrBwa(ch_bwa, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf4}
                        }
                        
                    } else {
                        hc2(ch_bwa, ch_mergebwabam).set {ch_mergevcf4}
                        if (params.run_vqsr) {
                            vqsrSnpBwa(ch_bwa, hc(ch_bwa, ch_mergebwabam)).set {ch_vqsrsnp}
                            vqsrIndelBwa(ch_bwa, hc(ch_bwa, ch_mergebwabam)).set {ch_vqsrindel}
                            gatherVcfsVqsrBwa(ch_bwa, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf4}
                        } 
                    }
                }
                vcfevalBwaGatk(ch_merge, ch_mergevcf4).set {ch_vcfevalBwaGatk}
                varStatsBwaGatk(ch_mergevcf4)

                splitVcfBwaGatk(ch_bwa, ch_gatk, ch_mergevcf4, chrs).eachvcf.set {ch_eachvcf4}
                vcfs4 = splitVcfBwaGatk.out.vcf.groupTuple()
                ch_eachbambwa.combine(ch_eachvcf4, by: [0,1]).set {ch_eachchr4}
                pvcfs4 = phaseBwaGatk(ch_bwa, ch_gatk, ch_eachchr4).phasedvcf.groupTuple()  
                lfs4 = phaseBwaGatk.out.lf.groupTuple()  
                hbs4 = phaseBwaGatk.out.hapblock.groupTuple()  
                stats4 = phaseBwaGatk.out.stat.groupTuple()  
                phaseCatBwaGatk(ch_bwa, ch_gatk, vcfs4.join(pvcfs4).join(lfs4).join(hbs4).join(stats4)).report.set {ch_phasereport4}//report
                pvcfs4.join(lfs4).join(hbs4).map { items ->
                    def id = items[0]
                    def paths = items[1..-1]
                    return [id, paths.flatten()]
                } set {ch_phaseallBwaGatk}
            }
        }

        //// stlfr bam stats & report
        if (params.align_tool.contains("lariat")) { 
            stlfrbam = ch_lariatbam 
            ch_mergebam = ch_mergeLariatBam
            ch_eachbam = ch_eachbamlariat
        } else { 
            stlfrbam = ch_stlfrBwaBam
            ch_mergebam = ch_mergebwabam
            ch_eachbam = ch_eachbambwa
        }

        bamdepth(ch_libstlfr, stlfrbam).set {ch_stlfrbamdepth}
        samtools_flagstat(ch_libstlfr, stlfrbam).set {ch_flagstat}
        samtools_stats(ch_libstlfr, stlfrbam).set {ch_stat}
        insertsize(ch_libstlfr, stlfrbam).insertsize.set {ch_insertsize} 
            
        if (!params.ref.startsWith('/')) {
            coverage(ch_merge, ch_mergebam).join(coverageMean(ch_merge, ch_mergebam)).set { ch_MergeGeneCov }
            coverageAvg(ch_PfGeneCov.join(ch_MergeGeneCov)).set {ch_avgCov}
        }

        samtools_depth(ch_libstlfr, stlfrbam).set {ch_depthreport}
        align_cat(ch_libstlfr, ch_flagstat.join(ch_stat).join(ch_depthreport).join(ch_insertsize)).set {ch_aligncatstlfr} //info
        stLFRQC(stlfrbam).report.set {ch_lfr}
        // if (!params.frombam) { fqstats_stlfr(ch_stlfrbssq) }//report 22

        ch_reports = Channel.empty() 
        if (params.align_tool.contains("bwa") && params.var_tool.contains("gatk")) {
            ch_vcf = ch_mergevcf4
            ch_phase = ch_phasereport4  
            phaseall = ch_phaseallBwaGatk
            hapcutstat = phaseCatBwaGatk.out.hapcutstat
            hb = phaseCatBwaGatk.out.hb
            reportBwaGatk(ch_bwa, ch_gatk, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalBwaGatk).join(ch_vcfevalPf)).collect().mix(ch_reports).set {ch_reports}
            
        }
        if (params.align_tool.contains("bwa") && params.var_tool.contains("dv")) {
            ch_vcf = ch_mergevcf3
            ch_phase = ch_phasereport3
            phaseall = ch_phaseallBwaDv
            hapcutstat = phaseCatBwaDv.out.hapcutstat
            hb = phaseCatBwaDv.out.hb
            // ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_PfGeneCov).join(ch_vcfevalBwaDv).join(ch_vcfevalPf).view()
            reportBwaDv(ch_bwa, ch_dv, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalBwaDv).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
        } 
        if (params.align_tool.contains("lariat") && params.var_tool.contains("gatk")) {
            ch_vcf = ch_mergevcf2
            ch_phase = ch_phasereport2
            phaseall = ch_phaseallLariatGatk
            hapcutstat = phaseCatLariatGatk.out.hapcutstat
            hb = phaseCatLariatGatk.out.hb

            reportLariatGatk(ch_lariat, ch_gatk, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalLariatGatk).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
        } 
        if (params.align_tool.contains("lariat") && params.var_tool.contains("dv")) {
            ch_vcf = ch_mergevcf
            phaseall = ch_phaseallLariatDv

            if (!params.ref.startsWith('/')) {
                hapcutstat = phaseCatLariatDv.out.hapcutstat
                ch_phase = ch_phasereport
                hb = phaseCatLariatDv.out.hb
                hapKaryotype(hb) // ${id}.haplotype.pdf

                reportLariatDv(ch_lariat, ch_dv, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalLariatDv).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
                report(ch_reports)
            } else {

            }
            
        } 
        

    } else {
        // from bam (lariat)
        if (!params.fromMergedBam) {
            println("!!! run CWGS from stlfr and pf bams")
            parse_sample_frombam(ch_input).bam.set {ch_bam}
            bam(ch_bam).stlfr.set {ch_lariatbam}
            bam.out.pf.set {ch_pfbam}       

            if (params.sampleBam) { 
                sampleBamPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} 
                sampleBamStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam).set {ch_lariatbam}
            }
            mapq_frombam(ch_pfbam).set {ch_pfbam}
            deepvariantv16BwaPf(ch_bwa, ch_pfbam).set {ch_pfdvvcf} 
            vcfevalPf(ch_libpf, ch_pfdvvcf).set {ch_vcfevalPf}
            coveragePf(ch_libpf, ch_pfbam).join(coverageMeanPf(ch_libpf, ch_pfbam)).set { ch_PfGeneCov }
            // mosdepthPf(ch_libpf, ch_pfbam).set {ch_}
            //depth bed
            // depth_pf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbed}

            //pf bam stats
            samtoolsFlagstatPf(ch_libpf, ch_pfbam).set {ch_flagstat2}
            samtoolsStatsPf(ch_libpf, ch_pfbam).set {ch_stat2}
            samtoolsDepthPf(ch_libpf, ch_pfbam).set {ch_depthreport2} 

            insertsizePf(ch_libpf, ch_pfbam).insertsize.set {ch_insertsize2} 
            alignCatPf(ch_libpf, ch_flagstat2.join(ch_stat2).join(ch_depthreport2).join(ch_insertsize2)).set {ch_aligncatpf} //info
            bamdepthPf(ch_libpf, ch_pfbam).set {ch_pfbamdepth}

            //merge bam
            intersectLariat(ch_lariat, ch_lariatbam.join(ch_pfbam)).set {ch_bed}
            mergeBamLariat(ch_lariat, ch_pfbam.join(ch_lariatbam).join(ch_bed)).set {ch_mergeLariatBam} 
        } else {
            println("!!! run CWGS from merge bam")
            parse_sample_frombam(ch_input).bam.set {ch_bam}
            bam2(ch_bam).stlfr.set {ch_lariatbam}
            bam2.out.merge.set {ch_mergeLariatBam}
        }
        
    

        if (params.ref == 'hs37d5') {chrs = (1..22).collect { it.toString() } + ['X', 'Y']}
        else {chrs = (1..22).collect { "chr$it" } + ["chrX", "chrY"]}

        if (params.var_tool.contains("gatk") && !params.use_megabolt && params.split_by_intervals) { //run gatk-hc with -L and then vqsr on each vcf
            gatk_interval().splitText().map { it.trim() }.collect().set {intervals}
        }


        if (params.var_tool.contains("gatk") && params.run_bqsr && !params.use_megabolt) { // not use megabolt
            bqsrStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam).set {ch_lariatbam}
        } else if (params.var_tool.contains("gatk") && params.run_bqsr && params.use_megabolt) { // use megabolt
            bqsrMegabolt(ch_lariatbam).set {ch_lariatbam}
        }
            

        // split stLFR bam for phasing
        splitBam4phasing(ch_lariat, ch_lariatbam, chrs).set {ch_eachbamlariat}

        if (params.var_tool.contains("dv")) {
            if (params.use_megabolt && params.dv_version != "v1.6") {dvMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
            else if (params.dv_version == "v1.6") {deepvariantv16(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
            
            vcfevalLariatDv(ch_merge, ch_mergevcf).set {ch_vcfevalLariatDv}//report 52 
            varStatsLariatDv(ch_mergevcf) //report 51

            //phase
            splitVcfLariatDv(ch_lariat, ch_dv, ch_mergevcf, chrs).eachvcf.set {ch_eachvcf}
            vcfs = splitVcfLariatDv.out.vcf.groupTuple()
            ch_eachbamlariat.combine(ch_eachvcf, by: [0,1]).set {ch_eachchr}
            pvcfs = phaseLariatDv(ch_lariat, ch_dv, ch_eachchr).phasedvcf.groupTuple()  
            lfs = phaseLariatDv.out.lf.groupTuple()  
            hbs = phaseLariatDv.out.hapblock.groupTuple()  
            stats = phaseLariatDv.out.stat.groupTuple()  
            phaseCatLariatDv(ch_lariat, ch_dv, vcfs.join(pvcfs).join(lfs).join(hbs).join(stats)).report.set {ch_phasereport}//report
            // phaseCat_cwx(ch_lariat, ch_dv, vcfs.join(pvcfs).join(lfs).join(hbs).join(stats))

            pvcfs.join(lfs).join(hbs).map { items ->
                def id = items[0]
                def paths = items[1..-1]
                return [id, paths]
            } set {ch_phaseallLariatDv}

        }
        if (params.var_tool.contains("gatk")) {
            if (params.use_megabolt) {
                if (params.run_vqsr) {
                    vqsrMegabolt(ch_lariat, hcMegabolt(ch_lariat, ch_mergeLariatBam)).set {ch_mergevcf2}
                } else {
                    hcMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf2}
                }
            } else {
                if (params.split_by_intervals) {
                    hcSplit(ch_lariat, ch_mergeLariatBam, intervals).set {ch_mergevcfSplit}
                    gatherVcfsHc(ch_lariat, ch_mergevcfSplit.groupTuple()).set {ch_mergevcf2}
                    if (params.run_vqsr) {
                        vqsrSnp(ch_lariat, ch_mergevcf2).set {ch_vqsrsnp}
                        vqsrIndel(ch_lariat, ch_mergevcf2).set {ch_vqsrindel}
                        gatherVcfsVqsr(ch_lariat, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf2}
                    }
                    
                } else {
                    hc(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf2}
                    if (params.run_vqsr) {
                        vqsrSnp(ch_lariat, hc(ch_lariat, ch_mergeLariatBam)).set {ch_vqsrsnp}
                        vqsrIndel(ch_lariat, hc(ch_lariat, ch_mergeLariatBam)).set {ch_vqsrindel}
                        gatherVcfsVqsr(ch_lariat, ch_vqsrsnp.join(ch_vqsrindel)).set{ch_mergevcf2}
                    } 
                }
            }
            vcfevalLariatGatk(ch_merge, ch_mergevcf2).set {ch_vcfevalLariatGatk}
            varStatsLariatGatk(ch_mergevcf2)

            //phase
            splitVcfLariatGatk(ch_lariat, ch_gatk, ch_mergevcf2, chrs).eachvcf.set {ch_eachvcf2}
            vcfs2 = splitVcfLariatGatk.out.vcf.groupTuple()
            ch_eachbamlariat.combine(ch_eachvcf2, by: [0,1]).set {ch_eachchr2}
            pvcfs2 = phaseLariatGatk(ch_lariat, ch_gatk, ch_eachchr2).phasedvcf.groupTuple()  
            lfs2 = phaseLariatGatk.out.lf.groupTuple()  
            hbs2 = phaseLariatGatk.out.hapblock.groupTuple()  
            stats2 = phaseLariatGatk.out.stat.groupTuple()  
            phaseCatLariatGatk(ch_lariat, ch_gatk, vcfs2.join(pvcfs2).join(lfs2).join(hbs2).join(stats2)).report.set {ch_phasereport2}//report
            pvcfs2.join(lfs2).join(hbs2).set {ch_phaseallLariatGatk}
        }
        
        stlfrbam = ch_lariatbam 
        ch_mergebam = ch_mergeLariatBam
        ch_eachbam = ch_eachbamlariat


        //stlfr bam stats
        bamdepth(ch_libstlfr, stlfrbam).set {ch_stlfrbamdepth}
        samtools_flagstat(ch_libstlfr, stlfrbam).set {ch_flagstat}
        samtools_stats(ch_libstlfr, stlfrbam).set {ch_stat}
        insertsize(ch_libstlfr, stlfrbam).insertsize.set {ch_insertsize} 
            
        coverage(ch_merge, ch_mergebam).join(coverageMean(ch_merge, ch_mergebam)).set { ch_MergeGeneCov }
        if (!params.fromMergedBam) {
            coverageAvg(ch_PfGeneCov.join(ch_MergeGeneCov)).set {ch_avgCov}
        }
        samtools_depth(ch_libstlfr, stlfrbam).set {ch_depthreport}
        align_cat(ch_libstlfr, ch_flagstat.join(ch_stat).join(ch_depthreport).join(ch_insertsize)).set {ch_aligncatstlfr} //info
        stLFRQC(stlfrbam).report.set {ch_lfr}
        // if (!params.frombam) { fqstats_stlfr(ch_stlfrbssq) }//report 22

        ch_reports = Channel.empty() 
        if (params.align_tool.contains("lariat") && params.var_tool.contains("gatk")) {
            ch_vcf = ch_mergevcf2
            ch_phase = ch_phasereport2
            phaseall = ch_phaseallLariatGatk
            hapcutstat = phaseCatLariatGatk.out.hapcutstat
            hb = phaseCatLariatGatk.out.hb
            if (!params.fromMergedBam) {
                reportLariatGatk1(ch_lariat, ch_gatk, ch_vcf.join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalLariatDv).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
            }
            
        } 
        if (params.align_tool.contains("lariat") && params.var_tool.contains("dv")) {
            ch_vcf = ch_mergevcf
            ch_phase = ch_phasereport
            phaseall = ch_phaseallLariatDv
            hapcutstat = phaseCatLariatDv.out.hapcutstat
            hb = phaseCatLariatDv.out.hb

            if (!params.fromMergedBam) {
                reportLariatDv1(ch_lariat, ch_dv, ch_vcf.join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalLariatDv).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
            }
        } 
        report(ch_reports)
        hapKaryotype(hb)            // ${id}.haplotype.pdf
        hapKaryotype_bak(hb)
    }
}
workflow.onComplete {
    println "CWGS started at: $workflow.start"
    println "CWGS completed at: $workflow.complete"
    println "The duration is: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow {
    println "Cmd line: $workflow.commandLine"
    println "CWGS started at: $workflow.start"
    CWGS()
}


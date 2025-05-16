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
    fq_stlfronly;
    fq_pfonly;
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
    kff;
    vg;
    changeid;
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
    combinebam;
    stLFRQC} from "${params.MOD}/align"
include { mapq;
    mapq as mapq_frombam            } from "${params.MOD}/bam"
include {
    frag1;
    frag2;
    fragstats                           } from "${params.MOD}/frag"
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
include { vcfstats } from "${params.MOD}/vcfstats"
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
    phase_cat as phaseCatLariatGatk;
    phaseCatRef;
    eachstat_phase;
    hapKaryotype;
    hapKaryotype_bak;
    ideogram;
    cumu;
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
    bamstats;
    eachstat_aligncat;
    align_cat as alignCatPf;
    align_catAll } from "${params.MOD}/bamstats"

include {
    report0 as reportBwaGatk;
    report0 as reportBwaDv;
    report0 as reportLariatGatk;
    report0 as reportLariatDv;
    reportref;
    report_stlfronly;
    report_stlfronly_ref;
    report01 as reportBwaGatk1;
    report01 as reportBwaDv1;
    report01 as reportLariatGatk1;
    report01 as reportLariatDv1;
    report_frombam_ref;
    report;
    html } from "${params.MOD}/report"


def ch_libpf = "pf"
def ch_libstlfr = "stlfr"
def ch_merge = "merge"

def ch_bwa = "bwa"
def ch_lariat = "lariat"

def ch_dv = "dv"
def ch_gatk = "gatk"
def ch_gatk3 = "gatk3"
def ch_gatk4 = "gatk4"

if (params.ref == 'hs37d5') {
    chrs = (1..22).collect { it.toString() } + ['X', 'Y']
} else if (params.ref == 'hg19' || params.ref == 'hg38') {
    if (params.chr == 'all') {
        chrs = (1..22).collect { "chr$it" } + ["chrX", "chrY"]
    } else {
        chrs = [params.chr]
    }
    
} else {
    getchrs().set { txt }
    chrs = txt.splitText().map { it.trim() }.collect()
}

workflow CWGS {
    println("!!! run CWGS from fastq")
    parse_sample (ch_input)
    .reads
    .set { ch_fq }

    if (params.stLFR_only) {
        fq_stlfronly(ch_fq).set {ch_stlfrfq}
        barcode_split(ch_stlfrfq).reads.set {ch_splitfq}
        splitLog = barcode_split.out.log
    } else if (params.PF_only) {
        fq_pfonly(ch_fq).set {ch_pffq}
    } else {
        fq(ch_fq)
        .stlfr
        .set {ch_stlfrfq}

        barcode_split(ch_stlfrfq).reads.set {ch_splitfq}
        splitLog = barcode_split.out.log

        fq.out.pf.set {ch_pffq}
    }
    
    if (!params.stLFR_only) {
        qc_pf(ch_libpf, ch_pffq).reads.set {ch_qcpffq} 
        readLenPf(qc_pf.out.bssq).set {ch_PFreadLen} 

        // pffq qc
        qc_pf.out.bssq.set {ch_pfbssq}
        fqcheckPf (ch_libpf, ch_qcpffq).set {ch_pffqcheck} 
        // .groupTuple(sort: sort_filenames('/'))
        // .set {ch_pffqcheck}

        fqdistPf (ch_libpf, ch_pffqcheck) //report 23 25
        fqstats_pf(ch_pfbssq) //report 21
    }
    
    def ch_stlfrsampledfq, ch_pfsampledfq
    if (params.sampleFq) { 
        if (!params.PF_only) {
            // sample stlfr fq
            qc_stlfr_stats(ch_splitfq).basecnt.set {ch_stlfrbasecount}
            qc_stlfr_stats.out.rlen.set {ch_stLFRreadLen}
            sampleStlfrFq(ch_stlfrbasecount.join(ch_stLFRreadLen).join(ch_splitfq).map {id, base, rlen, r1, r2 ->
                tuple(id, base, rlen, [r1, r2])}.transpose())
                .collect()
                .map {id, read1, idd, read2 -> 
                    def r1 = read1.toString().contains("_1.") ? read1 : read2
                    def r2 = read1.toString().contains("_1.") ? read2 : read1
                    tuple(id, r1, r2)}
                .set {ch_stlfrsampledfq}
        }
        if ( !params.stLFR_only) {
            // sample pf fq
            basecountPf(ch_pfbssq).set {ch_pfbasecount}
            samplePfFq(ch_pfbasecount.join(ch_PFreadLen).join(ch_qcpffq).map {id, base, rlen, r1, r2 ->
            tuple(id, base, rlen, [r1, r2])}.transpose())
                .collect()
                .map {id, read1, idd, read2 -> 
                    def r1 = read1.toString().contains("_1.") ? read1 : read2
                    def r2 = read1.toString().contains("_1.") ? read2 : read1
                    tuple(id, r1, r2)}
                .set {ch_pfsampledfq}
        }
        
    } else {
        if (!params.PF_only) { ch_stlfrsampledfq = ch_splitfq }
        if (!params.stLFR_only) { ch_pfsampledfq = ch_qcpffq }
    }
    
    // pf align
    if (!params.stLFR_only) {
        if (params.pfAligner == "bwa") {
            if (params.use_megabolt) {
                bwaMegaboltPf(ch_libpf, ch_pfsampledfq).set {ch_pfbam}
            } else {
                bwaPf(ch_libpf, ch_pfsampledfq).set {ch_pfsortbam}
                markdupPf(ch_libpf, ch_bwa, ch_pfsortbam).set {ch_pfbam}
                
            } 
        } else { // pf vg
            kff(ch_pfsampledfq).set {ch_kff}
            changeid(vg(ch_kff.join(ch_pfsampledfq))).set {ch_pfbam}
            //vg(ch_pfsampledfq).set {ch_pfbam}
        }
        if (params.var_tool.contains("gatk") && params.run_bqsr) { bqsrPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} }
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
    }
    
    if (!params.PF_only) {
        // stlfr align
        //// stlfr lariat
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
            if (params.stLFR_only) {
                ch_mergeLariatBam = ch_lariatbam
            } else {
                if (params.just_combine) {
                    combinebam(ch_lariatbam.join(ch_pfbam)).set {ch_mergeLariatBam}
                } else {
                    intersectLariat(ch_lariat, ch_lariatbam.join(ch_pfbam)).set {ch_bed}
                    mergeBamLariat(ch_lariat, ch_pfbam.join(ch_lariatbam).join(ch_bed)).set {ch_mergeLariatBam} 
                }
            }
            
            // call variant
            if (params.var_tool.contains("dv")) {
                if (params.use_megabolt && params.dv_version != "v1.6") {dvMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
                else if (params.dv_version == "v1.6") {deepvariantv16(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
                
                if (!params.ref.startsWith('/')) {
                    vcfevalLariatDv(ch_merge, ch_mergevcf).set {ch_vcfevalLariatDv}
                    varStatsLariatDv(ch_mergevcf) 
                }
                
                if (!params.PF_only) {
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
                        phaseCatRef(ch_lariat, ch_dv, txt, vcfs.join(pvcfs).join(lfs).join(hbs)).report.set {ch_phasereport}
                    }
                    
                    pvcfs.join(lfs).join(hbs).map { items ->
                        def id = items[0]
                        def paths = items[1..-1]
                        return [id, paths]
                    } set {ch_phaseallLariatDv}
                }
                
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

                if (!params.PF_only) {
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
            if (params.stLFR_only) {
                ch_mergebwabam = ch_stlfrBwaBam
            } else {
                intersectBwa(ch_bwa, ch_stlfrBwaBam.join(ch_pfbam)).set {ch_bed2}
                mergeBamBwa(ch_bwa, ch_pfbam.join(ch_stlfrBwaBam).join(ch_bed2)).set {ch_mergebwabam}
            }

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
            if (!params.stLFR_only) {
                coverageAvg(ch_PfGeneCov.join(ch_MergeGeneCov)).set {ch_avgCov}
            } 
        }

        samtools_depth(ch_libstlfr, stlfrbam).set {ch_depthreport}
        align_cat(ch_libstlfr, ch_flagstat.join(ch_stat).join(ch_depthreport).join(ch_insertsize)).set {ch_aligncatstlfr} //info
        stLFRQC(stlfrbam).report.set {ch_lfr}
        fragstats(splitLog.join(ch_lfr))
        // if (!params.frombam) { fqstats_stlfr(ch_stlfrbssq) }//report 22
        ch_reports = Channel.empty() 
        if (!params.stLFR_only) {
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
                    // ideogram(hb) 

                    reportLariatDv(ch_lariat, ch_dv, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalLariatDv).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).set {ch_report}
                    bamstats(ch_aligncatstlfr.join(ch_aligncatpf).join(ch_avgCov).join(ch_stlfrbamdepth).join(ch_pfbamdepth))
                    vcfstats(ch_vcf.join(ch_phase))
                    ch_report.collect().mix(ch_reports).set {ch_reports}
                    report(ch_reports)

                    // html
                    // html(ch_reports)
                } else {
                    ch_phase = ch_phasereport
                    hb = phaseCatRef.out.hb

                    reportref(ch_lariat, ch_dv, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
                    report(ch_reports)
                }
                
            } 
        } else { // stLFR only
            ch_vcf = ch_mergevcf
            phaseall = ch_phaseallLariatDv
            ch_phase = ch_phasereport
            if (!params.ref.startsWith('/')) {    
                report_stlfronly(ch_lariat, ch_dv, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_phase).join(ch_vcfevalLariatDv).join(ch_stlfrbamdepth)).collect().mix(ch_reports).set {ch_reports}
                report(ch_reports)
            } else { 
                report_stlfronly_ref(ch_lariat, ch_dv, ch_vcf.join(splitLog).join(ch_lfr).join(ch_aligncatstlfr).join(ch_phase).join(ch_stlfrbamdepth)).collect().mix(ch_reports).set {ch_reports}
                report(ch_reports)
            }

        }
    }
}

workflow CWGS_alignOnly {
    parse_sample (ch_input)
    .reads
    .set { ch_fq }

    if (params.stLFR_only) {
        fq_stlfronly(ch_fq).set {ch_stlfrfq}
        barcode_split(ch_stlfrfq).reads.set {ch_splitfq}
        splitLog = barcode_split.out.log
    } else if (params.PF_only) {
        fq_pfonly(ch_fq).set {ch_pffq}
    } else {
        fq(ch_fq)
        .stlfr
        .set {ch_stlfrfq}

        barcode_split(ch_stlfrfq).reads.set {ch_splitfq}
        splitLog = barcode_split.out.log

        fq.out.pf.set {ch_pffq}
    }
    
    if (!params.stLFR_only) {
        qc_pf(ch_libpf, ch_pffq).reads.set {ch_qcpffq} 
        readLenPf(qc_pf.out.bssq).set {ch_PFreadLen} 

        // pffq qc
        qc_pf.out.bssq.set {ch_pfbssq}
        fqcheckPf (ch_libpf, ch_qcpffq).set {ch_pffqcheck} 
        // .groupTuple(sort: sort_filenames('/'))
        // .set {ch_pffqcheck}

        fqdistPf (ch_libpf, ch_pffqcheck) //report 23 25
        fqstats_pf(ch_pfbssq) //report 21
    }
    
    def ch_stlfrsampledfq, ch_pfsampledfq
    if (params.sampleFq) { 
        if (!params.PF_only) {
            // sample stlfr fq
            qc_stlfr_stats(ch_splitfq).basecnt.set {ch_stlfrbasecount}
            qc_stlfr_stats.out.rlen.set {ch_stLFRreadLen}
            sampleStlfrFq(ch_stlfrbasecount.join(ch_stLFRreadLen).join(ch_splitfq).map {id, base, rlen, r1, r2 ->
                tuple(id, base, rlen, [r1, r2])}.transpose())
                .collect()
                .map {id, read1, idd, read2 -> 
                    def r1 = read1.toString().contains("_1.") ? read1 : read2
                    def r2 = read1.toString().contains("_1.") ? read2 : read1
                    tuple(id, r1, r2)}
                .set {ch_stlfrsampledfq}
        }
        if ( !params.stLFR_only) {
            // sample pf fq
            basecountPf(ch_pfbssq).set {ch_pfbasecount}
            samplePfFq(ch_pfbasecount.join(ch_PFreadLen).join(ch_qcpffq).map {id, base, rlen, r1, r2 ->
            tuple(id, base, rlen, [r1, r2])}.transpose())
                .collect()
                .map {id, read1, idd, read2 -> 
                    def r1 = read1.toString().contains("_1.") ? read1 : read2
                    def r2 = read1.toString().contains("_1.") ? read2 : read1
                    tuple(id, r1, r2)}
                .set {ch_pfsampledfq}
        }
        
    } else {
        if (!params.PF_only) { ch_stlfrsampledfq = ch_splitfq }
        if (!params.stLFR_only) { ch_pfsampledfq = ch_qcpffq }
    }
    
    // pf align
    if (!params.stLFR_only) {
        if (params.use_megabolt) {
            bwaMegaboltPf(ch_libpf, ch_pfsampledfq).set {ch_pfbam}
        } else {
            bwaPf(ch_libpf, ch_pfsampledfq).set {ch_pfsortbam}
            markdupPf(ch_libpf, ch_bwa, ch_pfsortbam).set {ch_pfbam}
        } 

        if (params.sampleBam) { sampleBamPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} }
        mapq(ch_pfbam).set {ch_pfbam}
    }
    
    if (!params.PF_only) {
        //// stlfr lariat
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
    }
    
}

workflow CWGS_frombam {
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
        if (!params.ref.startsWith('/')) {
            vcfevalPf(ch_libpf, ch_pfdvvcf).set {ch_vcfevalPf}
            coveragePf(ch_libpf, ch_pfbam).join(coverageMeanPf(ch_libpf, ch_pfbam)).set { ch_PfGeneCov }
        }
        

        
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
        
        if (!params.ref.startsWith('/')) {
            vcfevalLariatDv(ch_merge, ch_mergevcf).set {ch_vcfevalLariatDv}//report 52 
            varStatsLariatDv(ch_mergevcf) //report 51
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
            phaseCatRef(ch_lariat, ch_dv, txt, vcfs.join(pvcfs).join(lfs).join(hbs)).report.set {ch_phasereport}
        }

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
    
    if (!params.ref.startsWith('/')) {
        coverage(ch_merge, ch_mergebam).join(coverageMean(ch_merge, ch_mergebam)).set { ch_MergeGeneCov }
        if (!params.fromMergedBam) {
            coverageAvg(ch_PfGeneCov.join(ch_MergeGeneCov)).set {ch_avgCov}
        }
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
        

        if (!params.fromMergedBam) {
            if (!params.ref.startsWith('/')) {
                reportLariatDv1(ch_lariat, ch_dv, ch_vcf.join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_avgCov).join(ch_vcfevalLariatDv).join(ch_vcfevalPf).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
            } else {
                report_frombam_ref(ch_lariat, ch_dv, ch_vcf.join(ch_aligncatstlfr).join(ch_aligncatpf).join(ch_phase).join(ch_stlfrbamdepth).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
            }
            
        }
    } 
    report(ch_reports)
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
    if (params.alignOnly) {
        CWGS_alignOnly()
    } else if (params.frombam) {
        CWGS_frombam()
    } else {
        CWGS()
    }
    
}


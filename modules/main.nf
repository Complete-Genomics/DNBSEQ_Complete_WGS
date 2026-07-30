nextflow.enable.dsl=2

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.gatk_version != "v3" && params.gatk_version != "v4") { exit 1, 'wrong gatk version!'}
if (!params.use_megabolt && params.dv_version == "v0.7" ) { exit 1, 'dv v0.7 can only be used in megabolt!'}

include {
    parse_sample;
    parse_sample_frombam;
    bam;
    bam2                  } from "${params.MOD}/parseSample"
include { WF_align_pf;
        WF_align_stlfr;
        WF_align_stlfr2;
        bqsr as bqsrStlfrLariat;
        bqsrMegabolt;
        sampleBam as sampleBamPf;
        sampleBam as sampleBamStlfrLariat;
        stLFRQC             } from "${params.MOD}/align"
include {
    WF_mergebam;
    intersect as intersectLariat;
    mergeBam as mergeBamLariat } from "${params.MOD}/mergebam"
include {
    WF_callvariants;
    gatk_interval;
    gatherVcfsHc;
    gatherVcfsVqsr;
    deepvariant;
    deepvariant as dvBwaPf;
    dvMegabolt;
    hcMegabolt;
    vqsrMegabolt;
    hc;
    hcSplit;
    vqsrSnp;
    vqsrIndel } from "${params.MOD}/callvariants"
include { WF_phase;
        getchrs;
        splitBam4phasing;
        splitvcf as splitVcfLariatDv;
        splitvcf as splitVcfLariatGatk;
        phase as phaseLariatDv;
        phase as phaseLariatGatk;
        phaseCat as phaseCatLariatDv;
        phaseCat as phaseCatLariatGatk;
        phaseCat as phaseCatRef;
        ideogram;
        cumuplot;
        hapblock2bed        } from "${params.MOD}/phase"
include {
    coverage;
    coverage as coveragePf;
    coverageMean;
    coverageMean as coverageMeanPf;
    coverageAvg;
    cmrg_cnt                              } from "${params.MOD}/genedepth"
include { vcfstats } from "${params.MOD}/vcfstats"

include {
    vcfeval;
    vcfeval as vcfevalLariatDv;
    vcfeval as vcfevalLariatGatk;
    vcfeval as vcfevalPf;
    variant_fix } from "${params.MOD}/vcfeval"
include {
    variant_stats as varStatsLariatDv;
    variant_stats as varStatsLariatGatk } from "${params.MOD}/variantstats"
include { mapq as mapq_frombam } from "${params.MOD}/bam"

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
    samtools_depth as samtoolsDepthMerge;
    samtools_depth0;
    samtools_depth0 as samtoolsDepth0Pf;
    samtools_depth0 as samtoolsDepth0Merge;
    splitDepth as splitDepthPf;
    splitDepth as splitDepthStlfr;
    genomeDepth as genomeDepthPf;
    genomeDepth as genomeDepthStlfr;
    align_cat;
    eachstat_aligncat;
    align_cat as alignCatPf } from "${params.MOD}/bamstats"

include { vep;
    vep as vep_frombam;
    vep_data } from "${params.MOD}/annot"
include {gangstr} from "${params.MOD}/gangstr"
include {frag1; frag2} from "${params.MOD}/frag"
include {pangenie;
    pangenie_plot;
    pangenie_var_plot;
    pangenie_frombam } from "${params.MOD}/pangenie"
include {hlala} from "${params.MOD}/hlala"
include {WF_haplodenovo} from "${params.MOD}/haplodenovo"

include {report0;
    report0 as reportLariatGatk1;
    report0 as reportLariatDv;
    report;
    reportref;
    report_stlfronly;
    report_frombam_PFonly;
    FQC;
    html } from "${params.MOD}/report"

def ch_libpf = "pf"
def ch_libstlfr = "stlfr"
def ch_merge = "merge"

def ch_bwa = "bwa"
def ch_lariat = "lariat"
def ch_dv = "dv"
def ch_gatk = "gatk"


workflow CWGS {
    parse_sample(ch_input).set {ch_datapath}
    // ch_datapath.view()


    ch_datapath.branch { meta, path ->
        fq_stlfr    : meta.type == 'fq'     && meta.lib == 'stlfr'
            [meta.id, path]
        fq_stlfr2   : meta.type == 'fq'     && meta.lib == 'stlfr2'
            [meta.id, path]
        fq_pf       : meta.type == 'fq'     && meta.lib == 'pf'
            [meta.id, path]
        bam_stlfr   : meta.type == 'bam'    && meta.lib == 'stlfr'
            [meta.id, path]
        bam_stlfr2  : meta.type == 'bam'    && meta.lib == 'stlfr2'
            [meta.id, path]
        bam_pf      : meta.type == 'bam'    && meta.lib == 'pf'
            [meta.id, path]
    }.set {ch_data}

    // regular stLFR (PE, lariat/bwa)
    ch_data.bam_stlfr.mix(WF_align_stlfr(ch_data.fq_stlfr)).set {ch_stlfrbam_only}

    // stLFR2 SE 600/700bp (vg, barcodes at read end, header reformatted inside vg process)
    ch_data.bam_stlfr2.mix(WF_align_stlfr2(ch_data.fq_stlfr2)).set {ch_stlfr2bam}

    // combined stlfr-side bam channel fed into mergebam
    ch_stlfrbam_only.mix(ch_stlfr2bam).set {ch_stlfrbam}

    ch_data.bam_pf.mix(WF_align_pf(ch_data.fq_pf)).set {ch_pfbam}

    // merge (or just combine)
    WF_mergebam(ch_stlfrbam.join(ch_pfbam)).set {ch_mergebam}

    // bamstats – stLFRQC only for regular stLFR (not stlfr2 SE reads)
    stLFRQC(ch_stlfrbam_only).report.set {ch_lfr}
    samtoolsDepthMerge('merge', ch_mergebam).set {ch_depthreport}


    // call variants
    // Route per sample composition:
    //   both stlfr + pf → merge.bam
    //   stlfr-only      → stlfr.bam
    //   pf-only         → pf.bam
    ch_stlfrbam.map { id, bam -> tuple(id, 'stlfr') }
        .mix(ch_pfbam.map { id, bam -> tuple(id, 'pf') })
        .groupTuple(by: 0)
        .branch { id, srcs ->
            both       : srcs.size() == 2
            stlfr_only : srcs == ['stlfr']
            pf_only    : srcs == ['pf']
        }.set { ch_sampleType }

    ch_sampleType.both      .join(ch_mergebam).map { id, srcs, bam -> [id, bam] }
        .mix(
            ch_sampleType.stlfr_only.join(ch_stlfrbam).map { id, srcs, bam -> [id, bam] },
            ch_sampleType.pf_only   .join(ch_pfbam)   .map { id, srcs, bam -> [id, bam] }
        ).set { ch_callbam }

    WF_callvariants(ch_callbam).set {ch_mergevcf}
    
    // vcfeval (hg38)
    vcfeval('merge', ch_mergevcf).set {ch_vcfeval}

    // cmrg (hg38)
    coverage('merge', ch_mergebam)      .set {ch_cmrgMergebamhistbed}
    coverageMean('merge', ch_mergebam)  .set {ch_cmrgMergebammeanbed}
    
    // phase
    WF_phase(ch_stlfrbam.join(ch_mergevcf)).vcf.set {ch_phasedvcf}
    WF_phase.out.hb.set {hb}
    WF_phase.out.report.set {ch_phasereport}

    ideogram(hb)
    cumuplot(hb)
    hapblock2bed(hb)

    // per-haplotype overlapping-window de novo assembly
    if (params.run_haplodenovo) {
        WF_haplodenovo(ch_mergebam.join(ch_phasedvcf))
    }

    // sv
    pangenie(ch_data.fq_pf).set {ch_pangenie}
    pangenie_plot(ch_pangenie.join(hb))
    pangenie_var_plot(ch_pangenie)
    gangstr(ch_mergebam)
    hlala(ch_mergebam)

    // annot
    vep_data(vep_frombam(ch_phasedvcf).html)

    // WF_report()
    if (params.ref == 'hg38' || params.ref.contains('GRCh38')) {
        report0(ch_phasedvcf.join(ch_lfr).join(ch_cmrgMergebamhistbed).join(ch_cmrgMergebammeanbed).join(ch_depthreport).join(ch_phasereport)).collect().set {ch_reports}
        html(report(ch_reports))
    }
}

workflow CWGS_frombam {
    if (params.ref.startsWith('/')) {
        getchrs().set { txt }
        chrs = txt.splitText().map { it.trim() }.collect()
    }
    if (!params.fromMergedBam) {
        println("!!! run CWGS from stlfr and pf bams")
        parse_sample_frombam(ch_input).set { ch_parsed_frombam }
        ch_parsed_frombam.bam.set    { ch_bam }
        ch_parsed_frombam.fq_pf.set  { ch_fq_pf_frombam }   // pf1/pf2 from samplesheet for pangenie
        bam(ch_bam).stlfr.set {ch_lariatbam}
        bam.out.pf.set {ch_pfbam}       

        if (params.sampleBam) { 
            sampleBamPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} 
            sampleBamStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam).set {ch_lariatbam}
        }
        mapq_frombam(ch_pfbam).set {ch_pfbam}
        dvBwaPf(ch_bwa, ch_pfbam).set {ch_pfdvvcf} 
        if (!params.ref.startsWith('/')) {
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
        

        // //pf bam stats
        // samtoolsFlagstatPf(ch_libpf, ch_pfbam).set {ch_flagstat2}
        // samtoolsStatsPf(ch_libpf, ch_pfbam).set {ch_stat2}
        // samtoolsDepthPf(ch_libpf, ch_pfbam).set {ch_depthreport2} 

        // insertsizePf(ch_libpf, ch_pfbam).insertsize.set {ch_insertsize2} 
        // bamdepthPf(ch_libpf, ch_pfbam).set {ch_pfbamdepth}

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
        frag2(frag1(ch_eachbamlariat).groupTuple())

    if (params.var_tool.contains("dv")) {
        if (params.use_megabolt && params.dv_version == 'v0.6') {dvMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
        else {deepvariant(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
        
        if (params.ref == 'hg38' || params.ref.contains('GRCh38')) {
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
            vep_frombam(phaseCatLariatDv.out.phasedvcf)
            phaseCatLariatDv.out.phasedvcf.set {ch_phasedvcf}
        } else {
            phaseCatRef(ch_lariat, ch_dv, txt, vcfs.join(pvcfs).join(lfs).join(hbs)).report.set {ch_phasereport}
            phaseCatRef.out.phasedvcf.set {ch_phasedvcf}
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
    // bamdepth(ch_libstlfr, stlfrbam).set {ch_stlfrbamdepth}
    // samtools_flagstat(ch_libstlfr, stlfrbam).set {ch_flagstat}
    // samtools_stats(ch_libstlfr, stlfrbam).set {ch_stat}
    // insertsize(ch_libstlfr, stlfrbam).insertsize.set {ch_insertsize} 

    
    samtools_depth(ch_libstlfr, stlfrbam).set {ch_depthreport}
    stLFRQC(stlfrbam).report.set {ch_lfr}

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
        ch_vcf = ch_phasedvcf
        ch_phase = ch_phasereport
        hb = phaseCatLariatDv.out.hb
        ideogram(hb)
        cumuplot(hb)
        hapblock2bed(hb)

        if (!params.fromMergedBam) {
            if (params.ref == 'hg38' || params.ref.contains('GRCh38')) {
                coverage(ch_merge, ch_mergebam).set {ch_cmrgMergebamhistbed}
                coverageMean(ch_merge, ch_mergebam).set {ch_cmrgMergebammeanbed}

                vep_data(vep_frombam(ch_phasedvcf).html)

                gangstr(ch_mergebam)
                hlala(ch_mergebam)

                pangenie(ch_fq_pf_frombam).set {ch_pangenie}
                pangenie_plot(ch_pangenie.join(hb))
                pangenie_var_plot(ch_pangenie)

                if (params.run_haplodenovo) {
                    WF_haplodenovo(ch_mergebam.join(ch_phasedvcf))
                }

                reportLariatDv(ch_lariat, ch_dv, ch_vcf.join(ch_lfr).join(ch_cmrgMergebamhistbed).join(ch_cmrgMergebammeanbed).join(ch_depthreport).join(ch_phase)).set {ch_report}

                ch_report.collect().mix(ch_reports).set {ch_reports}
                report(ch_reports)

                ch_reports.mix(vep_data.out, hlala.out, cumuplot.out, pangenie_var_plot.out).collect().set {ch_flg}
                html(ch_flg)

            } else {
                reportref(ch_lariat, ch_dv, ch_vcf.join(ch_lfr).join(ch_depthreport).join(ch_phase)).set {ch_report}

                ch_report.collect().mix(ch_reports).set {ch_reports}
                report(ch_reports)
            }
        }
    } 
}
workflow CWGS_frombam_stLFRonly {
    if (params.ref.startsWith('/')) {
        getchrs().set { txt }
        chrs = txt.splitText().map { it.trim() }.collect()
    }
    if (!params.fromMergedBam) {
        println("!!! run CWGS from stlfr bams")
        parse_sample_frombam(ch_input).bam.set {ch_bam}
        bam(ch_bam).stlfr.set {ch_lariatbam}

        if (params.sampleBam) { 
            sampleBamStlfrLariat(ch_libstlfr, ch_lariat, ch_lariatbam).set {ch_lariatbam}
        }

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
    frag2(frag1(ch_eachbamlariat).groupTuple())

    if (params.var_tool.contains("dv")) {
        if (params.use_megabolt && params.dv_version == "v0.6" ) {dvMegabolt(ch_lariat, ch_mergeLariatBam).set {ch_mergevcf}}
        else {deepvariant(ch_lariat, ch_lariatbam).set {ch_vcf}}
        

        //phase
        splitVcfLariatDv(ch_lariat, ch_dv, ch_vcf, chrs).eachvcf.set {ch_eachvcf}
        vcfs = splitVcfLariatDv.out.vcf.groupTuple()
        ch_eachbamlariat.combine(ch_eachvcf, by: [0,1]).set {ch_eachchr}
        pvcfs = phaseLariatDv(ch_lariat, ch_dv, ch_eachchr).phasedvcf.groupTuple()  
        lfs = phaseLariatDv.out.lf.groupTuple()  
        hbs = phaseLariatDv.out.hapblock.groupTuple()  
        stats = phaseLariatDv.out.stat.groupTuple()  
        if (!params.ref.startsWith('/')) {
            phaseCatLariatDv(ch_lariat, ch_dv, vcfs.join(pvcfs).join(lfs).join(hbs).join(stats)).report.set {ch_phasereport}//report
            phaseCatLariatDv.out.phasedvcf.set {ch_phasedvcf}
        } else {
            phaseCatRef(ch_lariat, ch_dv, txt, vcfs.join(pvcfs).join(lfs).join(hbs)).report.set {ch_phasereport}
            phaseCatRef.out.phasedvcf.set {ch_phasedvcf}
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
    ch_eachbam = ch_eachbamlariat


    //stlfr bam stats
    bamdepth(ch_libstlfr, stlfrbam).set {ch_stlfrbamdepth}
    samtools_flagstat(ch_libstlfr, stlfrbam).set {ch_flagstat}
    // samtools_stats(ch_libstlfr, stlfrbam).set {ch_stat}
    // insertsize(ch_libstlfr, stlfrbam).insertsize.set {ch_insertsize} 

    samtools_depth(ch_libstlfr, stlfrbam).set {ch_stlfrbamdepthreport}
    // align_cat(ch_libstlfr, ch_flagstat.join(ch_stat).join(ch_depthreport).join(ch_insertsize)).set {ch_aligncatstlfr} //info
    stLFRQC(stlfrbam).report.set {ch_lfr}

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
        ch_vcf = ch_phasedvcf
        ch_phase = ch_phasereport
        

        if (!params.fromMergedBam) {
            report_stlfronly(ch_lariat, ch_dv, ch_vcf.join(ch_lfr).join(ch_flagstat).join(ch_phase).join(ch_stlfrbamdepth).join(ch_stlfrbamdepthreport)).set {ch_report}
            ch_report.collect().mix(ch_reports).set {ch_reports}
            FQC(report(ch_reports))
        }
    } 
}
workflow CWGS_frombam_PFonly {
    println("!!! run CWGS from pf bams")
    
    parse_sample_frombam(ch_input).bam.set {ch_bam}
    bam(ch_bam).stlfr.set {ch_lariatbam}
    bam.out.pf.set {ch_pfbam}       

    if (params.sampleBam) {  sampleBamPf(ch_libpf, ch_bwa, ch_pfbam).set {ch_pfbam} }
    mapq_frombam(ch_pfbam).set {ch_pfbam}
    dvBwaPf(ch_bwa, ch_pfbam).set {ch_pfdvvcf} 
    vcfevalPf(ch_libpf, ch_pfdvvcf).set {ch_vcfevalPf}
    coveragePf(ch_libpf, ch_pfbam).join(coverageMeanPf(ch_libpf, ch_pfbam)).set { ch_PfGeneCov }

    //pf bam stats
    samtoolsFlagstatPf(ch_libpf, ch_pfbam).set {ch_flagstat2}
    samtoolsStatsPf(ch_libpf, ch_pfbam).set {ch_stat2}
    samtoolsDepthPf(ch_libpf, ch_pfbam).set {ch_depthreport2} 

    insertsizePf(ch_libpf, ch_pfbam).insertsize.set {ch_insertsize2} 
    alignCatPf(ch_libpf, ch_flagstat2.join(ch_stat2).join(ch_depthreport2).join(ch_insertsize2)).set {ch_aligncatpf} //info
    bamdepthPf(ch_libpf, ch_pfbam).set {ch_pfbamdepth}

    ch_vcf = ch_pfdvvcf
    ch_reports = Channel.empty() 
    report_frombam_PFonly(ch_vcf.join(ch_aligncatpf).join(ch_pfbamdepth)).collect().mix(ch_reports).set {ch_reports}
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
    if (params.frombam) {
        if (params.PF_only) {
            CWGS_frombam_PFonly()
        } else if (params.stLFR_only) {
            CWGS_frombam_stLFRonly()
        } else {
            CWGS_frombam()
        }
    } else {
        CWGS()
    }
}

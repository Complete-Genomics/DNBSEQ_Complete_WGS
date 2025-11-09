nextflow.enable.dsl=2

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.gatk_version != "v3" && params.gatk_version != "v4") { exit 1, 'wrong gatk version!'}
if (!params.use_megabolt && params.dv_version == "v0.7" ) { exit 1, 'dv v0.7 can only be used in megabolt!'}

include { parse_sample      } from "${params.MOD}/parseSample"                   
include { WF_align_pf;
        WF_align_stlfr;
        stLFRQC             } from "${params.MOD}/align"    
include { WF_mergebam       } from "${params.MOD}/mergebam"
include { WF_callvariants   } from "${params.MOD}/callvariants"
include { WF_phase;
        ideogram;
        cumuplot            } from "${params.MOD}/phase"
include {
    coverage;
    coverage as coveragePf;
    coverageMean;
    coverageMean as coverageMeanPf;
    coverageAvg;
    cmrg_cnt                              } from "${params.MOD}/genedepth"
include { vcfstats } from "${params.MOD}/vcfstats"

include {
    vcfeval
    variant_fix } from "${params.MOD}/vcfeval"

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
    vep_data } from "${params.MOD}/annot"
include {gangstr} from "${params.MOD}/gangstr"
include {pangenie;
    pangenie_plot;
    pangenie_var_plot } from "${params.MOD}/pangenie"
include {hlala} from "${params.MOD}/hlala"

include {report0;
    report;
    reportref;
    html } from "${params.MOD}/report"


workflow CWGS {
    parse_sample(ch_input).set {ch_datapath}
    // ch_datapath.view()


    ch_datapath.branch { meta, path ->
        fq_stlfr    : meta.type == 'fq'     && meta.lib == 'stlfr'  
            [meta.id, path] 
        fq_pf       : meta.type == 'fq'     && meta.lib == 'pf'     
            [meta.id, path]
        bam_stlfr   : meta.type == 'bam'    && meta.lib == 'stlfr'  
            [meta.id, path]
        bam_pf      : meta.type == 'bam'    && meta.lib == 'pf'  
            [meta.id, path]
    }.set {ch_data}

    // ch_data.fq_stlfr.view()
    // ch_data.bam_pf.view()

    ch_data.bam_stlfr   .mix(WF_align_stlfr(ch_data.fq_stlfr))  .set {ch_stlfrbam}
    ch_data.bam_pf      .mix(WF_align_pf(ch_data.fq_pf))        .set {ch_pfbam}

    // merge (or just combine)
    WF_mergebam(ch_stlfrbam.join(ch_pfbam)).set {ch_mergebam}

    // bamstats
    stLFRQC(ch_stlfrbam).report.set {ch_lfr}
    samtoolsDepthMerge('merge', ch_mergebam).set {ch_depthreport}


    // call variants
    // ch_mergebam.mix(ch_stlfrbam).first().view() //set {ch_bam}
    // ch_bam.view()
    WF_callvariants(ch_mergebam).set {ch_mergevcf}
    
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

    // sv
    pangenie(ch_data.fq_pf).set {ch_pangenie}
    pangenie_plot(ch_pangenie.join(hb))
    pangenie_var_plot(ch_pangenie)
    gangstr(ch_mergebam)
    hlala(ch_mergebam)

    // annot
    vep_data(vep(ch_phasedvcf).html)

    // WF_report()
    if (params.ref == 'hg38' || params.ref.contains('GRCh38')) {
        report0(ch_phasedvcf.join(ch_lfr).join(ch_cmrgMergebamhistbed).join(ch_cmrgMergebammeanbed).join(ch_depthreport).join(ch_phasereport)).collect().set {ch_reports}
        html(report(ch_reports))
    } else {
        reportref(ch_phasedvcf.join(ch_lfr).join(ch_depthreport).join(ch_phasereport)).collect().set {ch_reports}
        report(ch_reports)
    }
    
}

workflow {
    CWGS()
}

#!/usr/bin/env nextflow
params.outdir_bams="bam_int_files"
params.outdir_pangenome="pangenome_mapped_bams"
params.fastq_dir = "/home/luisluna/links/scratch/genotype_GA4K_strs_strkit/fastq_out"
params.gbz_file="pangenomes_resources/hprc-v1.1-mc-chm13.gbz"
params.index_files = "pangenomes_resources/hprc-v1.1-mc-chm13.{gbz,dist,min}"
params.vg_container  = "vg_tools/vg.sif"
params.pbtk_container = "dockers/pbtk.sif"
process merge_bams {
    container params.pbtk_container
    publishDir "${params.outdir_bams}/${family_id}", mode: 'symlink'
    scratch true

    tag "$sample_id"

    input:
        tuple val(family_id), val(sample_id), path(bam_files)
    output:
        //path "${tag}_merged.bam", emit: merge_bam
        tuple val(family_id), val(sample_id), path("${sample_id}_merged_2.bam"), path("${sample_id}_merged_2.bam.pbi"), emit: merge_bam
        //path "${tag}_merged.bam.bai", emit: merge_bam_index

    script:
    """
     echo "Using sample_id: ${sample_id}"
     echo "Merging BAM files: ${bam_files.join(', ')}"

     pbmerge -o ${sample_id}_merged_2.bam ${bam_files} -j 8

    """
}

process map_to_pangenome {
    container params.vg_container
    publishDir "${params.outdir_pangenome}/${family_id}", mode: 'symlink'
    scratch true
    label "big_mem"
    tag "$family_id"

    input:
    tuple val(family_id), path(fastq_files)   // all 3 staged, processed sequentially
    path index_files

    output:
    tuple val(family_id), path("${family_id}_mapped.gam"), emit: mapped_gam

    script:
    def gbz = index_files.find { it.name.endsWith('.gbz') }
    // Sort to ensure consistent order: -01, -02, -03
    def sorted_fqs = fastq_files.sort { it.name }
    """
    # Stream FASTQs one after another into vg giraffe
    # No temp file needed — cat streams directly
    cat ${sorted_fqs.join(' ')} | vg giraffe \
        -b hifi \
        -Z ${gbz} \
        -f - \
        -p \
        --threads ${task.cpus} \
        > ${family_id}_mapped.gam
    """
} 

workflow {
    index_files_ch = Channel.fromPath(params.index_files).collect()

    merge_bams_files = Channel.fromPath("/project/6007512/shared/C3G/projects/ga4k_from_cedar/pacbio/aligned2/UNMC-000034*/*.bam")
        .map { bam_path ->
            def parent_dir = bam_path.getParent().getName()
            def sample_id = parent_dir.tokenize('.')[0]           // e.g. UNMC-0034-01
            def family_id = sample_id.tokenize('-')[0..1].join('-') // e.g. UNMC-0034
            tuple(family_id, sample_id, bam_path)
        }
        .groupTuple(by:1)
        .map { family_id_list, sample_id, bam_list ->
        tuple(family_id_list[0], sample_id, bam_list)  // take first element of family_id list
    }
        .view { it -> "Grouped BAM files by family: ${it}" }
        .set { grouped_bams }

   

    

    fastq_ch = Channel
        .fromPath("${params.fastq_dir}/**/*_merged_2.fastq.gz")
        .map { fq ->
            def sample_id = fq.getName().tokenize('_')[0..2].join('-')  // e.g. UNMC-000034-01
            def family_id = sample_id.tokenize('-')[0..1].join('-')     // e.g. UNMC-000034
            tuple(family_id, sample_id, fq)
        }
        .view { "FASTQ file: ${it}" }

    ////Run processes
    //Bam merging
    merged = merge_bams(grouped_bams)
    merged.merge_bam.view()
    //Align to pangenome
    mapped = map_to_pangenome(fastq_ch, index_files_ch)

}  
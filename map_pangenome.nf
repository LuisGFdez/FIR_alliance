#!/usr/bin/env nextflow
params.outdir_bams="bam_int_files"
process merge_bams {

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

     dockers/apptainer exec pbtk.sif pbmerge ${bam_files.join(' ')} -j ${task.cpus} > ${sample_id}_merged_2.bam 

    """
}

workflow {
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

    merged = merge_bams(grouped_bams)
    
    merged.merge_bam.view()

}  
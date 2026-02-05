#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.bed_file = "/home/luisluna/gentotype_ga4k_strkit/strkit_files/hg38_ver6_modified.sorted.bed"
params.bed_file_trgt = "/home/luisluna/scratch/hg38_ver6_modified.sorted.bed"

params.snps_vcf="/home/luisluna/scratch/genotype_ga4k_strkit/strkit_int_files/00-common_all.vcf.gz"
params.snps_vcf_index="/home/luisluna/scratch/genotype_ga4k_strkit/strkit_int_files/00-common_all.vcf.gz.tbi"
params.reference_genome = "/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa"
params.outdir    = "results_vcf"
params.outdir_fastas="fasta_int_file"
params.outdir_bams="bam_int_files"
params.outdir_trgt="trgt_genotypes"

process bgzip_index_fasta {
    publishDir params.outdir_fastas, mode: 'symlink' 
    
    input:
        path reference_genome
    output:
        path "${reference_genome}.gz", emit: fasta_gz
        path "${reference_genome}.fai", emit: fasta_fai
        path "${reference_genome}.gz.gzi", emit: fasta_gzi
    script:
    """
    
    echo "Compressing reference genome with bgzip ${reference_genome}"
    samtools faidx ${reference_genome}
    ls
    bgzip -@ $task.cpus ${reference_genome} > ${reference_genome}.gz
    ls
    samtools faidx --gzi-idx ${reference_genome}.gz
    """
}

process merge_bams {

    publishDir params.outdir_bams, mode: 'symlink'
    //scratch true

    tag { tag }

    input:
        tuple val(tag), path(bam_files)
    output:
        //path "${tag}_merged.bam", emit: merge_bam
        tuple path("${tag}_merged.bam"), path("${tag}_merged.bam.bai"), emit: merge_bam
        //path "${tag}_merged.bam.bai", emit: merge_bam_index

    script:
    """
    echo "Using tag: ${tag}"
    echo "Merging BAM files: ${bam_files.join(', ')}"
    samtools merge -@ $task.cpus ${tag}_merged.bam ${bam_files.join(' ')}
    samtools index ${tag}_merged.bam
    """
}

process genotype_strkit {

    tag {input_bam.simpleName}
    
    publishDir params.outdir, mode: 'copy'

    input:
        tuple path (input_bam), path (input_bam_index)
        path bed_tr_file
        path snp_vcf_file
        path snp_vcf_index_file
        path reference_genome

    output:
        path "${input_bam.simpleName}_strkit_genotypes.vcf" , emit: vcf_output
        path "${input_bam.simpleName}_strkit_genotypes_sorted.vcf.gz" , emit: vcf_compressed
        path "${input_bam.simpleName}_strkit_genotypes_sorted.vcf.gz.csi", emit: vcf_index
        

    script:
    """
    echo "Genotyping with StrKit for BAM: ${input_bam}"
    echo "Bam index_file: ${input_bam_index}"

    strkit call ${input_bam} --hq --realign\
    --ref ${reference_genome}\
    --loci ${bed_tr_file} \
    --vcf ${input_bam.simpleName}_strkit_genotypes.vcf\
    --incorporate-snvs ${snp_vcf_file} \
    --seed 183 \
    --processes 10 \
    --no-tsv

    bcftools sort ${input_bam.simpleName}_strkit_genotypes.vcf | bcftools view -O z -o ${input_bam.simpleName}_strkit_genotypes_sorted.vcf.gz
    bcftools index ${input_bam.simpleName}_strkit_genotypes_sorted.vcf.gz

    """
}
process genotype_TRGT {
    tag {input_bam.simpleName}

    publishDir params.outdir_trgt, mode: 'copy'

    input:
        tuple path (input_bam), path (input_bam_index)
        path bed_tr_file
        path reference_genome
        path reference_genome_gz
        path reference_genome_index
        path reference_genome_gzi_index
    output:
        tuple path("${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz"), path("${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz.csi"), emit: vcf_file_trgt
        tuple path("${input_bam.simpleName}.spanning.bam"), path("${input_bam.simpleName}.spanning.sorted.bam"), path("${input_bam.simpleName}.spanning.sorted.bam.bai"), emit: spanning_bam
    
    script:
    """
    echo "Genotyping with TRGT for BAM: ${input_bam}"       
    /home/luisluna/scratch/trgt-v5.0.0-x86_64-unknown-linux-gnu/trgt genotype --genome ${reference_genome} \
    --repeats ${bed_tr_file} \
    --reads ${input_bam} \
    --output-prefix ${input_bam.simpleName}_trgt_genotypes \
    --threads $task.cpus

    bcftools sort -Oz -o ${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz ${input_bam.simpleName}_trgt_genotypes.vcf.gz
    bcftools index ${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz

    samtools sort -o ${input_bam.simpleName}.spanning.sorted.bam ${input_bam.simpleName}.spanning.bam
    samtools index ${input_bam.simpleName}.spanning.sorted.bam
    """

}

process mendelian_inheritance {

    publishDir params.outdir, mode: 'copy'

    input:
         path genotype_str_vcf 
         path genotype_str_vcf_gz
         path genotype_str_vcf_csi

    output:
         path "mi_report.json", emit: mendelian_report

     script:

     """
     echo "Performing Mendelian Inheritance Check on VCF files: ${genotype_str_vcf.join(', ')}"
     ###echo ${genotype_str_vcf_gz.join(', ')}  
     ###echo ${genotype_str_vcf_csi.join(', ')}    
     echo "Vcf files: ${genotype_str_vcf_gz[0]}"
     echo "Vcf files: ${genotype_str_vcf_gz[1]}"
     echo "Vcf files: ${genotype_str_vcf_gz[2]}"

     strkit mi --caller strkit-vcf \
     --json mi_report.json \
     --caller strkit-vcf \
     --mismatch-out-mi strict \
     --hist \
     --test x2 \
     --sig-level 0.05 \
     ${genotype_str_vcf_gz[0]} \
     ${genotype_str_vcf_gz[1]} \
     ${genotype_str_vcf_gz[2]}

     """
}

process targt_denovo {
    publishDir params.outdir_trgt, mode: 'copy'

    input:
         path reference_genome
         path bed_tr_file
         tuple path (genotype_TRGT_vcfs)
         tuple path (genotype_TRGT_bams)

    output:
         path "trgt_denovo_report.tsv", emit: trgt_denovo_report

     script:

     """
     echo "Performing De Novo Mutation Detection on VCF files: ${genotype_trgt_vcf.join(', ')}"

     TRGT-denovo trio --reference ${reference_genome}\
     --bed repeat.bed
     --father ${genotype_str_vcf_gz[1]} \
     --mother ${genotype_str_vcf_gz[2]} \
     --child ${genotype_str_vcf_gz[0]} \
     --out trgt_denovo_report.tsv

     """
}
///project/6007512/wcheung/pacbio
/*
 * Workflow definition
 */
workflow {
    merge_bams_files = Channel.fromPath('/project/6007512/wcheung/pacbio/aligned2/cmh002019*/*.bam')
        .map { bam_path ->
            def parent_dir = bam_path.getParent().getName()
            def tag = parent_dir.tokenize('.')[0]
            tuple(tag, bam_path)
        }
        .groupTuple()
        .map { tag, bam_list -> tuple(tag, bam_list) }
        //.view { it -> "tag: ${it[0]}, bam files: ${it[1].collect{ it.getName() }}" }

    // To use in merge_bams process:
    merge_bams_files.set {grouped_bams}
    grouped_bams.view()
    bed_tr_file = Channel.value(file(params.bed_file))
    snp_files   = Channel.value(file(params.snps_vcf))
    snps_index  = Channel.value(file(params.snps_vcf_index))
    reference_genome = Channel.fromPath(params.reference_genome)

    //bed_tr_file = Channel.fromPath(params.bed_file)
    bed_tr_file_trgt = Channel.value(file(params.bed_file_trgt))             
                        // .view { it -> "bed_tr_file: ${it}" } 
    //reference_genome=Channel.fromPath(params.reference_genome)
                           // .view { it -> "ref_genome: ${it}" }
    //snp_files=Channel.fromPath(params.snps_vcf)
    //snps_index=Channel.fromPath(params.snps_vcf_index)
                           // .view { it -> "snp_files: ${it}" }        
    ////define processes              
    bgzip_index_fasta(reference_genome)  
    merged = merge_bams(grouped_bams)

    //merged.merge_bam.view { it -> "Merged BAM: ${it}" }
   
    //merged.merge_bam.view()
    //merged.merge_bam_index.view { it -> "Merged BAM Index: ${it}" }

    genotype_strkit(merged.merge_bam,bed_tr_file,snp_files,snps_index,bgzip_index_fasta.out.fasta_gz)
    sorted_genotypes=genotype_strkit.out.vcf_compressed.toSortedList { a, b -> a[0] <=> b[0] }.view()
    genotype_TRGT(merged.merge_bam, bed_tr_file_trgt,reference_genome,bgzip_index_fasta.out.fasta_gz,bgzip_index_fasta.out.fasta_fai,bgzip_index_fasta.out.fasta_gzi)
    
    genotype_str_vcf=genotype_strkit.out.vcf_output.collect()//.view { it -> "Genotyped VCF files: ${it}" }  
    genotype_str_vcf_gz=genotype_strkit.out.vcf_compressed.collect()
    genotype_str_vcf_csi=genotype_strkit.out.vcf_index.collect()

    genotype_TRGT_vcfs = genotype_TRGT.out.vcf_file_trgt.collect()
    genotype_TRGT_vcfs.view { it -> "Genotyped TRGT VCF files: ${it}" }
    //genotype_str_vcf_trgt = genotype_TRGT_vcfs.map { it[0] }.view { it -> "Genotyped TRGT VCF files: ${it}" }
    //genotype_str_vcf_gz_trgt = genotype_TRGT_vcfs.map { it[1] }.view { it -> "Genotyped TRGT VCF GZ files: ${it}" }
    //genotype_str_vcf_csi_trgt = genotype_TRGT_vcfs.map { it[2] }.view { it -> "Genotyped TRGT VCF CSI files: ${it}" }

   
    //genotype_TRGT_bams = genotype_TRGT.out.spanning_bam.collect()
   
    //genotype_STR_bams_trgt = genotype_TRGT_bams.map { it[0] }.view { it -> "Spanning BAM files: ${it}" }
    //genotype_STR_bams_sorted_trgt = genotype_TRGT_bams.map { it[1] }.view { it -> "Spanning Sorted BAM files: ${it}" }
    //genotype_STR_bams_index_trgt = genotype_TRGT_bams.map { it[2] }.view { it -> "Spanning BAM Index files: ${it}" }

    mendelian_inheritance(genotype_str_vcf,genotype_str_vcf_gz,genotype_str_vcf_csi)

    //targt_denovo(bgzip_index_fasta.out[0].first(), bed_tr_file.first(),genotype_TRGT_vcfs,genotype_TRGT_bams)

}    
//nextflow clean $(nextflow log -q) -f
//nextflow run genotype_strkit.nf -bg -with-trace -resume
//source ./GA4K/bin/activate
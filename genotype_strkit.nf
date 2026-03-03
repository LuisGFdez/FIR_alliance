#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
params.bed_file = "tr_bed_files/hg38_ver6_strkit.sorted.bed"
params.bed_file_trgt = "tr_bed_files/hg38_ver6_trgt.sorted.bed"

params.snps_vcf="strkit_int_files/00-All.vcf.gz"
params.snps_vcf_index="strkit_int_files/00-All.vcf.gz.tbi"
params.reference_genome = "/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa"
params.reference_genome_index= "/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa.fai"
params.outdir    = "results_vcf"
params.outdir_fastas="fasta_int_file"
params.outdir_bams="bam_int_files"
params.outdir_trgt="trgt_genotypes"
params.bed_files="tr_bed_files"
process bgzip_index_fasta {
    publishDir params.outdir_fastas, mode: 'symlink' 
    
    input:
        path reference_genome
    output:
        path "${reference_genome}.gz", emit: fasta_gz
        path "${reference_genome}.gz.fai", emit: fasta_fai
        path "${reference_genome}.gz.gzi", emit: fasta_gzi
    script:
    """
    echo "Compressing reference genome with bgzip ${reference_genome}"
    bgzip -@ ${task.cpus} ${reference_genome}
    echo "Compressing reference genome with faidx --gzi ${reference_genome}"
    samtools faidx ${reference_genome}.gz
    """
}
process create_tr_catlog {
    publishDir "${params.bed_files}", mode: 'symlink'
    //container params.snakemake_container
    scratch true
    input:
        path reference_genome
        path reference_genome_index    // needs the .fai index

    output:
        path "chr*.fa", emit: chr_fastas

    script:
    """
    # Split into one FA per chromosome
    cut -f1 ${reference_genome_index} \
        | grep -E '^chr([0-9]+|X|Y)\$' \
        | while read chr; do
            samtools faidx ${reference_genome} \$chr > \${chr}.fa
          done
    """
}
process merge_bams {

    publishDir "${params.outdir_bams}/${family_id}", mode: 'symlink'
    scratch true

    tag "$sample_id" 

    input:
        tuple val(family_id), val(sample_id), path(bam_files)
    output:
        //path "${tag}_merged.bam", emit: merge_bam
        tuple val(family_id), val(sample_id), path("${sample_id}_merged.bam"), path("${sample_id}_merged.bam.bai"), emit: merge_bam
        //path "${tag}_merged.bam.bai", emit: merge_bam_index

    script:
    """
    echo "Using sample_id: ${sample_id}"
    echo "Merging BAM files: ${bam_files.join(', ')}"

    if [ \$(echo "${bam_files}" | wc -w) -eq 1 ]; then
        samtools view -o ${sample_id}_merged.bam ${bam_files}
    else
        samtools merge -@ ${task.cpus} ${sample_id}_merged.bam ${bam_files.join(' ')}
    fi
        samtools index -@ ${task.cpus} ${sample_id}_merged.bam
    """
}

process genotype_strkit {

    tag "$sample_id"
    scratch true
    publishDir "${params.outdir}/${family_id}", mode: 'copy'

    input:
        tuple val(family_id), val(sample_id), path (input_bam), path (input_bam_index)
        path bed_tr_file
        path snp_vcf_file
        path snp_vcf_index_file
        path reference_genome

    output:
        path "${sample_id}_strkit_genotypes.vcf" , emit: vcf_output
        path "${sample_id}_strkit_genotypes_sorted.vcf.gz" , emit: vcf_compressed
        path "${sample_id}_strkit_genotypes_sorted.vcf.gz.csi", emit: vcf_index
        

    script:
    """
    echo "Genotyping with StrKit for BAM: ${input_bam}"
    echo "Bam index_file: ${input_bam_index}"

    strkit call ${input_bam} --hq --realign\
    --ref ${reference_genome}\
    --loci ${bed_tr_file} \
    --vcf ${sample_id}_strkit_genotypes.vcf\
    --incorporate-snvs ${snp_vcf_file} \
    --seed 183 \
    --processes 10 \
    --no-tsv

    bcftools sort ${sample_id}_strkit_genotypes.vcf | bcftools view -O z -o ${sample_id}_strkit_genotypes_sorted.vcf.gz
    bcftools index ${sample_id}_strkit_genotypes_sorted.vcf.gz

    """
}
process genotype_TRGT {

    tag "$input_bam.simpleName"
    scratch true
    label 'big_mem'
    publishDir "${params.outdir_trgt}/${family_id}", mode: 'copy'

    input:
        tuple val(family_id), val(sample_id), path (input_bam), path (input_bam_index)
        path bed_tr_file
        path reference_genome
        path reference_genome_index
        path reference_genome_gz
        path reference_genome_gzi_index
    output:
        tuple val(sample_id), path("${sample_id}_trgt_genotypes.vcf.gz") , path("${sample_id}_trgt_genotypes_sorted.vcf.gz"), path("${sample_id}_trgt_genotypes_sorted.vcf.gz.csi"), emit: vcf_file_trgt
        tuple val(sample_id), path("${sample_id}_trgt_genotypes.spanning.bam"), path("${sample_id}_trgt_genotypes.spanning.sorted.bam"), path("${sample_id}_trgt_genotypes.spanning.sorted.bam.bai"), emit: spanning_bam
    
    script:
    """
    echo "Genotyping with TRGT for BAM: ${input_bam}"       
    /home/luisluna/links/scratch/genotype_GA4K_strs_strkit/trgt-v5.0.0-x86_64-unknown-linux-gnu/trgt genotype --genome ${reference_genome} \
    --repeats ${bed_tr_file} \
    --reads ${input_bam} \
    --output-prefix ${sample_id}_trgt_genotypes \
    --threads $task.cpus

    bcftools sort -Oz -o ${sample_id}_trgt_genotypes_sorted.vcf.gz ${sample_id}_trgt_genotypes.vcf.gz
    bcftools index ${sample_id}_trgt_genotypes_sorted.vcf.gz

    samtools sort -o ${sample_id}_trgt_genotypes.spanning.sorted.bam ${sample_id}_trgt_genotypes.spanning.bam
    samtools index ${sample_id}_trgt_genotypes.spanning.sorted.bam
    """

}

process mendelian_inheritance {

    publishDir "${params.outdir}/${family_ids}", mode: 'copy'

    input:
         val family_ids
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
    publishDir "${params.outdir_trgt}/${family_ids}", mode: 'copy'
    scratch true
    input:
         val family_ids
         path reference_genome
         path reference_genome_index
         path bed_tr_file
         tuple path(child_files), path(father_files), path(mother_files)
         tuple path(child_bams), path(father_bams), path(mother_bams)

    output:
         path "trgt_denovo_report.tsv", emit: trgt_denovo_report

     script:

     """
     echo "Reference genome: ${reference_genome}"
     echo "BED file: ${bed_tr_file}"
     echo "child VCF: ${child_files[1]}"
     echo "father VCF: ${father_files[1]}"
     echo "mother VCF: ${mother_files[1]}"
     /home/luisluna/links/scratch/genotype_GA4K_strs_strkit/trgt-denovo-v0.3.0-x86_64-unknown-linux-gnu/trgt-denovo trio --reference ${reference_genome}\
     --bed ${bed_tr_file} \
     --father-vcf ${father_files[1]} \
     --mother-vcf ${mother_files[1]} \
     --child-vcf ${child_files[1]} \
     --father-bam ${father_bams[1]} \
     --mother-bam ${mother_bams[1]} \
     --child-bam ${child_bams[1]} \
     --out trgt_denovo_report.tsv

     """
}
///project/6007512/wcheung/pacbio
/*
 * Workflow definition
 */
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
        

    // To use in merge_bams process:

    //grouped_bams.view()
    bed_tr_file = Channel.value(file(params.bed_file))
    snp_files   = Channel.value(file(params.snps_vcf))
    snps_index  = Channel.value(file(params.snps_vcf_index))
    reference_genome = Channel.value(file(params.reference_genome))
    reference_genome_index = Channel.value(file(params.reference_genome_index))

    bed_tr_file_trgt = Channel.value(file(params.bed_file_trgt))             
                        // .view { it -> "bed_tr_file: ${it}" } 
    //reference_genome=Channel.fromPath(params.reference_genome)
                           // .view { it -> "ref_genome: ${it}" }
    //snp_files=Channel.fromPath(params.snps_vcf)
    //snps_index=Channel.fromPath(params.snps_vcf_index)
                           // .view { it -> "snp_files: ${it}" }        
    ////define processes              
    bgzip_index_fasta(reference_genome) 
    create_tr_catlog(reference_genome,reference_genome_index)  
    merged = merge_bams(grouped_bams)

    //merged.merge_bam.view { it -> "Merged BAM: ${it}" }
   
    merged.merge_bam.view()
    family_ids = merged.merge_bam.map { family_id, sample_id, bam, bai ->
    family_id
     }.unique()
    //merged.merge_bam_index.view { it -> "Merged BAM Index: ${it}" }

    genotype_strkit(merged.merge_bam,bed_tr_file,snp_files,snps_index,bgzip_index_fasta.out.fasta_gz)
    genotype_TRGT(merged.merge_bam, bed_tr_file_trgt,reference_genome,reference_genome_index,bgzip_index_fasta.out.fasta_gz,bgzip_index_fasta.out.fasta_gzi)
    
    genotype_str_vcf=genotype_strkit.out.vcf_output.collect()//.view { it -> "Genotyped VCF files: ${it}" } 
    genotype_str_vcf_gz=genotype_strkit.out.vcf_compressed.collect()
    genotype_str_vcf_csi=genotype_strkit.out.vcf_index.collect()
    flat_vcfs = genotype_str_vcf_gz.flatten()//.view { it -> "flattened VCF GZ files: ${it}" }
    sorted_genotypes = flat_vcfs.toSortedList { a, b ->
              def na = (a.name =~ /-(\d+)_/)[0][1].toInteger()
              def nb = (b.name =~ /-(\d+)_/)[0][1].toInteger()
              na <=> nb }
    sorted_genotypes.view { it -> "Sorted Genotyped VCF files: ${it}" }
    genotype_TRGT.out.vcf_file_trgt
    .toSortedList { a, b -> a[0] <=> b[0] }   // sort by sample_id (val)
    .map { sorted ->
        sorted.collect { sample_id, vcf, sorted_vcf, csi -> [vcf, sorted_vcf, csi] }
    }
    .set { trio_vcfs }

    genotype_TRGT.out.spanning_bam
    .toSortedList { a, b -> a[0] <=> b[0] }   // sort by sample_id (val)
    .map { sorted ->
        sorted.collect { sample_id, bam, sorted_bam, bai -> [bam, sorted_bam, bai] }
    }
    .set { trio_bams }
    

    mendelian_inheritance(family_ids,genotype_str_vcf,sorted_genotypes,genotype_str_vcf_csi)

    targt_denovo(family_ids,reference_genome, reference_genome_index,bed_tr_file_trgt,trio_vcfs,trio_bams)

}    
//nextflow clean $(nextflow log -q) -f
//nextflow run genotype_strkit.nf -bg -with-trace -resume
//source ./GA4K/bin/activate
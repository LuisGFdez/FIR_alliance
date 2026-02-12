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

process merge_bams {

    publishDir params.outdir_bams, mode: 'symlink'
    scratch true

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
    samtools merge -@ ${task.cpus} ${tag}_merged.bam ${bam_files.join(' ')}
    samtools index -@ ${task.cpus} ${tag}_merged.bam
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
        path reference_genome_index
        path reference_genome_gz
        path reference_genome_gzi_index
    output:
        tuple path("${input_bam.simpleName}_trgt_genotypes.vcf.gz") , path("${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz"), path("${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz.csi"), emit: vcf_file_trgt
        tuple path("${input_bam.simpleName}_trgt_genotypes.spanning.bam"), path("${input_bam.simpleName}_trgt_genotypes.spanning.sorted.bam.bai"), path("${input_bam.simpleName}_trgt_genotypes.spanning.sorted.bam.bai"), emit: spanning_bam
    
    script:
    """
    echo "Genotyping with TRGT for BAM: ${input_bam}"       
    /home/luisluna/links/scratch/genotype_GA4K_strs_strkit/trgt-v5.0.0-x86_64-unknown-linux-gnu/trgt genotype --genome ${reference_genome} \
    --repeats ${bed_tr_file} \
    --reads ${input_bam} \
    --output-prefix ${input_bam.simpleName}_trgt_genotypes \
    --threads $task.cpus

    bcftools sort -Oz -o ${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz ${input_bam.simpleName}_trgt_genotypes.vcf.gz
    bcftools index ${input_bam.simpleName}_trgt_genotypes_sorted.vcf.gz

    samtools sort -o ${input_bam.simpleName}_trgt_genotypes.spanning.sorted.bam ${input_bam.simpleName}_trgt_genotypes.spanning.bam
    samtools index ${input_bam.simpleName}_trgt_genotypes.spanning.sorted.bam
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
         tuple val(sample_id_vcf) , path (genotype_TRGT_vcfs)
         tuple val(sample_id_bam) , path (genotype_TRGT_bams)

    output:
         path "trgt_denovo_report.tsv", emit: trgt_denovo_report

     script:

     """
     echo "Performing De Novo Mutation Detection on VCF files: ${genotype_TRGT_vcfs.join(', ')}"
     echo "Reference genome: ${reference_genome}"
     echo "BED file: ${bed_tr_file}"
     echo "Sample IDs: ${sample_id_vcf.join(', ')}"
     echo "child VCF: ${sample_id_vcf[0]}"
     echo "father VCF: ${sample_id_vcf[1]}"
     echo "mother VCF: ${sample_id_vcf[2]}"
     /home/luisluna/links/scratch/genotype_GA4K_strs_strkit/trgt-denovo-v0.3.0-x86_64-unknown-linux-gnu/trgt-denovo trio --reference ${reference_genome}\
     --bed ${bed_tr_file}
     --father ${sample_id_vcf[2]} \
     --mother ${sample_id_vcf[1]} \
     --child ${sample_id_vcf[0]} \
     --out trgt_denovo_report.tsv

     """
}
///project/6007512/wcheung/pacbio
/*
 * Workflow definition
 */
workflow {
    merge_bams_files = Channel.fromPath("/project/6007512/shared/C3G/projects/ga4k_from_cedar/pacbio/aligned2/cmh002019*/*.bam")
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
    merged = merge_bams(grouped_bams)

    //merged.merge_bam.view { it -> "Merged BAM: ${it}" }
   
    //merged.merge_bam.view()
    //merged.merge_bam_index.view { it -> "Merged BAM Index: ${it}" }

    genotype_strkit(merged.merge_bam,bed_tr_file,snp_files,snps_index,bgzip_index_fasta.out.fasta_gz)
    genotype_TRGT(merged.merge_bam, bed_tr_file_trgt,reference_genome,reference_genome_index,bgzip_index_fasta.out.fasta_gz,bgzip_index_fasta.out.fasta_gzi)
    
    genotype_str_vcf=genotype_strkit.out.vcf_output.collect().view { it -> "Genotyped VCF files: ${it}" } 
    genotype_str_vcf_gz=genotype_strkit.out.vcf_compressed.collect()
    genotype_str_vcf_csi=genotype_strkit.out.vcf_index.collect()
    flat_vcfs = genotype_str_vcf_gz.flatten().view { it -> "flattened VCF GZ files: ${it}" }
    sorted_genotypes = flat_vcfs.toSortedList { a, b ->
    def na = (a.name =~ /-(\d+)_/)[0][1].toInteger()
    def nb = (b.name =~ /-(\d+)_/)[0][1].toInteger()
    na <=> nb }
    sorted_genotypes.view { it -> "Sorted Genotyped VCF files: ${it}" }

    genotype_TRGT.out.vcf_file_trgt.flatten().map {file ->
                                                          def sample = (file.name =~ /^(cmh\d+-\d+)/)[0][1]
                                                          tuple(sample, file)
                                                          }
                                              .groupTuple()
                                              .toSortedList { a, b -> a[0] <=> b[0] }
                                              //.flatten()
                                              //.buffer(size:2)
                                              .map { sorted_list ->
                                                      def samples = sorted_list.collect { it[0] }   // [cmh002019-01, cmh002019-02, cmh002019-03]
                                                      def files   = sorted_list.collect { it[1] }   // [[files01], [files02], [files03]]
                                                      tuple(samples, files)
                                                   }
                                                      // group all 3 files per sample
                                              .set{group_trgt_vcfs}
    // Reshape channel to include both sorted vcf AND its index
    group_trgt_vcfs
        .map { sample_ids, file_groups ->
        def child  = file_groups[0]   // [merged.vcf.gz, sorted.vcf.gz, sorted.vcf.gz.csi]
        def father = file_groups[1]
        def mother = file_groups[2]
        tuple(
            sample_ids,
            child[1],    // child sorted vcf
            child[2],    // child index .csi
            father[1],   // father sorted vcf
            father[2],   // father index .csi
            mother[1],   // mother sorted vcf
            mother[2]    // mother index .csi
        )
    }
    .set { trio_vcfs }
   trio_vcfs.view()                                                 //.view { sample, files -> "Sample: $sample\n  Files: $files" }
   // group_trgt_vcfs.view()
   trio_vcfs.view{it->"Group VCF files: ${it[0][1]}"}
   trio_vcfs.view{it->"Group VCF files: ${it[0]}"}
   genotype_TRGT.out.spanning_bam.flatten().map { file ->
                                                          def sample = (file.name =~ /^(cmh\d+-\d+)/)[0][1]
                                                          tuple(sample, file)
                                                          }
                                              .groupTuple()
                                              .toSortedList { a, b -> a[0] <=> b[0] }
                                              .map { sorted_list ->
                                                      def samples = sorted_list.collect { it[0] }   // [cmh002019-01, cmh002019-02, cmh002019-03]
                                                      def files   = sorted_list.collect { it[1] }   // [[files01], [files02], [files03]]
                                                      tuple(samples, files)
                                                   }
                                                      // group all 3 files per sample
                                              .set{group_trgt_bams} 
    group_trgt_bams.map { sample_ids, file_groups ->
        def child  = file_groups[0]   // [merged.vcf.gz, sorted.vcf.gz, sorted.vcf.gz.csi]
        def father = file_groups[1]
        def mother = file_groups[2]
        tuple(
            sample_ids,
            child[0],    // child spanning bam
            child[1],    // child spanning bam index
            child[2],    // child spanning bam index
            father[0],   // father spanning bam
            father[1],   // father spanning bam index
            father[2],   // father spanning bam index
            mother[0],   // mother spanning bam
            mother[1],   // mother spanning bam index
            mother[2]    // mother spanning bam index
        )
    }
    .set { trio_bams }
    trio_bams.view{it->"Group BAM files: ${it[0]}"}
    trio_bams.view{it->"Group BAM files: ${it[1]}"} 
    


    mendelian_inheritance(genotype_str_vcf,sorted_genotypes,genotype_str_vcf_csi)

    targt_denovo(reference_genome, bed_tr_file_trgt,trio_vcfs,trio_bams)

}    
//nextflow clean $(nextflow log -q) -f
//nextflow run genotype_strkit.nf -bg -with-trace -resume
//source ./GA4K/bin/activate
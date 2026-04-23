#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { checkDirs; checkFiles } from './utils/helpers.nf'

user=System.getenv('USER')
params.scratch="/scratch/alpine/$USER"
params.workdir="${params.scratch}/cef_strains/work"
params.datadir="${params.scratch}/cef_strains/20210421_gDNA_ONT_2/fastq_pass/cat_fastq"
params.plasmid_size_select = "10000000" // Passed to seqkit seq -m {}. 
                                        // Karla-Vasco used 3M file size using 
                                        // find but I'm subsampling input fastqs

infile_pat="${params.datadir}/50x*.fastq.gz" // smallest sample
//infile_pat="${params.datadir}/*.*.fastq.gz" // only subsampled by frac
//infile_pat="${params.datadir}/*.fastq.gz" // everything

// DIAMOND ANNOTATIONS
diamond_database_path="${launchDir}/HMDARG/hmd-arg.dmnd"
hmdarg_database_path="${launchDir}/HMDARG/arg_v5_linear.fasta"
hmdarg_annotation_path="${launchDir}/HMDARG/annotations_hmd-arg.csv"

// OUTPUT/PUBLISH DIRS
// use built-in launchDir so that everything 
// is relative to where nextflow is invoked
params.nanoq_stats_outdir="${launchDir}/nanoq_stats"
// assembly flye/medaka
params.assembly_stats_outdir="${launchDir}/assemblies/stats"
params.medaka_consensus_outdir="${launchDir}/assemblies/medaka"
params.medaka_gaps_outdir="${launchDir}/assemblies/medaka_gaps"
params.assembly_outdir="${launchDir}/assemblies" 
params.prokka_out="${launchDir}/assemblies/prokka"
// prodigal
params.prodigal_outdir="${launchDir}/prodigal"
params.prodigal_coords_gbk_outdir="${params.prodigal_outdir}/gbk"
params.prodigal_proteins_outdir="${params.prodigal_outdir}/proteins"
// resistome
params.resistome_outdir="${launchDir}/resistome"

// show and check settings
c_bold = "\033[1m"
c_red = "\033[0;31m"
c_green = "\033[0;32m"
c_yellow = "\033[0;33m"
c_grey = "\033[90m"
c_reset = "\033[0m"

println("${c_bold}Environment and Directory settings${c_reset}");
println("\tuser: $USER")
println("\tparams.scratch: $params.scratch")
println("\tparams.workdir: $params.workdir")
println("\tparams.datadir: $params.datadir")


/*
println("${c_bold}Check for directories${c_reset}");
checkDirs([
    "params.scratch": params.scratch,
    "params.workdir": params.workdir,
    "params.datadir": params.datadir]
)
println("${c_bold}Check for diamond database sources${c_reset}");
checkFiles([
    "diamond_database_path": diamond_database_path,
    "hmdarg_database_path": hmdarg_database_path,
    "hmdarg_annotation_path": hmdarg_annotation_path
], true)
*/

process PORECHOP {
    tag "${sample_name}"
    conda params.nanop_env

    input:
    tuple val(sample_name), path(fastq_gz)

    output:
    tuple val(sample_name), path("${sample_name}_chop.fastq.gz")

    script:
    """
    porechop -t ${task.cpus} -i ${fastq_gz} -o ${sample_name}_chop.fastq.gz
    """

}

process NANOFILT {
    tag "${sample_name}"
    conda params.nanop_env
    publishDir "${params.nanoq_stats_outdir}", pattern: "*_nanoq_stats.txt"

    input: 
    tuple val(sample_name), path(chop_fastq_gz)

    output:
    tuple val(sample_name), path("${sample_name}_filt.fastq.gz"), emit: filtered
    path("${sample_name}_nanoq_stats.txt"), emit: stats

    script:
    """
    pigz -dc -p ${task.cpus} ${chop_fastq_gz} | nanoq -q 10 -t ${task.cpus} -r \
        ${sample_name}_nanoq_stats.txt -o ${sample_name}_filt.fastq.gz
    """
}

process FLYE {
    tag "${sample_name}"
    conda params.nanop_env
    label 'highres' 

    input:
    tuple val(sample_name), path(nanofilt_fastq_gz)

    output:
    tuple val(sample_name), path(nanofilt_fastq_gz), path("30-contigger/contigs.fasta")

    script:
    """
    flye --nano-raw ${nanofilt_fastq_gz} \
        --out-dir . \
        --threads ${task.cpus} \
        --stop-after contigger \
        --meta

    """
}

process MEDAKA {
    tag "${sample_name}"
    container params.medaka_sif
    publishDir "${params.medaka_consensus_outdir}", pattern: "medaka_out/*.consensus.fasta"
    publishDir "${params.medaka_gaps_outdir}", pattern: "medaka_out/*.consensus.fasta.gaps_in_draft_coords.bed"

    input: 
    tuple val(sample_name), path(nanofilt_fastq_gz), path(contig_fasta)

    output:
    tuple val(sample_name), path("medaka_out/*.consensus.fasta"), emit: consensus
    path "medaka_out/*.consensus.fasta.gaps_in_draft_coords.bed", emit: gaps


    script:
    """
    samtools faidx ${contig_fasta}

    medaka_consensus \
    -i ${nanofilt_fastq_gz} \
    -d ${contig_fasta} \
    -o medaka_out \
    -t ${task.cpus} \
    -f \
    -b ${25 * task.cpus} \
    -m r941_prom_hac_g507

    ln -v medaka_out/consensus.fasta medaka_out/${sample_name}.consensus.fasta
    ln -v medaka_out/consensus.fasta.gaps_in_draft_coords.bed medaka_out/${sample_name}.consensus.fasta.gaps_in_draft_coords.bed
    """

}

process REMOVE_PLASMIDS {
    tag "${sample_name}"
    conda params.nanop_env
    publishDir "${params.assembly_outdir}", pattern: "filtered/*.fasta"

    input:
    tuple val(sample_name), path(full_assembly_fasta)

    output:
    tuple val(sample_name), path("filtered/*.fasta")

    script:
    """
    mkdir -v filtered
    seqkit seq -m ${params.plasmid_size_select} < ${full_assembly_fasta} > filtered/${sample_name}.fasta
    """

    stub:
    def outfile = "filtered/${full_assembly_fasta.name}"
    """
    mkdir -v filtered
    touch ${outfile}
    """
}


process PRODIGAL {
    tag "${sample_name}"
    conda params.nanop_env
    publishDir "${params.prodigal_coords_gbk_outdir}", pattern: "*_coords.gbk"
    publishDir "${params.prodigal_proteins_outdir}", pattern: "*_proteins.faa"

    input:
    tuple val(sample_name), path(assembly_fasta)

    output:
    path "*_coords.gbk", emit: prodigal_gbk
    path "*_proteins.faa", emit: prodigal_faa

    script:
    """
    prodigal -i ${assembly_fasta} \
    -o ${sample_name}_coords.gbk \
    -a ${sample_name}_proteins.faa \
    -p meta
    """
}
process PROKKA {
    tag "${sample_name}"
    conda params.nanop_env
    publishDir params.prokka_out, mode: 'copy'

    input:
    tuple val(sample_name), path(assembly_fasta)

    output:
    path "*_prokka"

    script:
    """
    prokka ${assembly_fasta} \
    --outdir ${sample_name}_prokka --force \
    --prefix ${sample_name}_prokka \
    --genus Escherichia \
    --evalue 0.001 \
    --cpus ${task.cpus} \
    --addgenes  
    """
}
process DIAMOND {
    tag "${sample_name}"
    conda params.nanop_env

    input:
    tuple val(sample_name), path(proteins_faa)
    path diamond_db

    output:
    tuple val(sample_name), path("${sample_name}_hmdarg_matches.sam")

    script:
    """
    diamond blastp \
    -d ${diamond_db} \
    -q ${proteins_faa} \
    -o ${sample_name}_hmdarg_matches.sam \
    -f 101 \
    --id 80 \
    -k 1 \
    --max-hsps 1
    """
}

process RESISTOME {
    tag "${sample_name}"
    conda params.nanop_env
    publishDir "${params.resistome_outdir}", pattern: "*_hmdarg_*", mode: 'copy'

    input:
    tuple val(sample_name), path(assembly_fa)

    script:
    """
    resistome \
    -ref_fp ${assembly_fa} \
    -annot_fp ${params.annotation} \
    -sam_fp ${sample_name}_hmdarg_matches.sam \
    -gene_fp ${sample_name}_hmdarg_gene.tsv \
    -group_fp ${sample_name}_hmdarg_mobility.tsv \
    -mech_fp ${sample_name}_hmdarg_mech.tsv \
    -class_fp ${sample_name}_hmdarg_class.tsv \
    -t 80
    """
}

// this will still work with nextflow -preview, whereas the Chanel().view will not
println "Files in datadir: " + file(infile_pat).collect { it.name }

workflow {
    // set up source channels
    ch_cat_fastq = Channel
                    .fromPath(infile_pat, checkIfExists: true)
                    .map { pth -> 
                            tuple(
                              pth.name.replace('.fastq.gz', ''), // isolate sample name to pass through channels 
                              pth
                            ) 

                    } 
    //ch_diamond_db = DIAMOND_DB.broadcast() // download/format databases

    // create an assembly from each input (cat_fastq)
    ch_cat_fastq | PORECHOP | NANOFILT 
    FLYE(NANOFILT.out.filtered) | MEDAKA 
    PROKKA(MEDAKA.out.consensus)
    PRODIGAL(MEDAKA.out.consensus)
    
    //ch_assembly = REMOVE_PLASMIDS(MEDAKA.out.consensus)
    // branch on the assembly
    //ch_assembly | PRODIGAL //| DIAMOND(ch_diamond_db) | RESISTOME
    //ch_assembly | PROKKA
}

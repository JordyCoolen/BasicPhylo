#!/usr/bin/env nextflow

// (Tool) params
params.reads                = "$baseDir/data/raw_input/20300134001-1_MB_R{1,2}.fastq.gz"
params.outDir               = "$baseDir"
params.threads              = 4
// SKA Fastq
params.kmer                 = "15"
params.quality_score        = "20"
params.coverage             = "4"
// SKA Align
params.proportion           = "0.3"
// Iq Tree
params.bootstrap            = 1000
params.contree              = "false" // is not working

// Parsing the input parameters
outDir                      = "$params.outDir"
def samplename                = file("$params.reads").simpleName[0].split('_')[0]
threads                     = "$params.threads"
// SKA Fastq
kmer                        = "$params.kmer"
quality_score               = "$params.quality_score"
coverage                    = "$params.coverage"
// SKA Align
proportion                  = "$params.proportion"
// Iq Tree
bootstrap                   = "$params.bootstrap"
contree                     = "$params.contree"

// Database
database                    = "$baseDir/db/ref_15mer_${proportion}_variants.skf"

// Path prefixes
r_folder                    = "$baseDir/R/"

// Tool paths
cluster_sample              = "$baseDir/R/clustering.py"
create_tree                 = "$baseDir/R/treemaker.R"
to_json                     = "$baseDir/python/finalize.py"

parameters = "${kmer},${quality_score},${coverage},${proportion},${bootstrap},${contree}"

log.info """
BASICPHYLO V0.1

============INPUT===============
reads      : $params.reads
filename   : $samplename

~~~~~~~~~~~~parameters~~~~~~~~~~
threads    : $threads
bootstrap  : $bootstrap
proportion : $proportion
kmer       : $kmer
quality_score   : $quality_score
coverage        : $coverage

~~~~~~~~~~~Databases~~~~~~~~~~~~
database    : $database

~~~~~~~~~~~Authors~~~~~~~~~~~~~~
        Paul Verhoeven
        J.P.M. Coolen
================================
"""

Channel
    .fromFilePairs( params.reads )
    .ifEmpty{ "cannot find read pairs in path"}
    .set{ raw_reads }

// Pairs forward and reverse fastq files and puts them in to pairReads channel.
process pairReads {
	input:
	set pair_ID, file(reads) from raw_reads

	output:
	set file("${reads[0]}"), file("${reads[1]}") into pair_reads
	"""
	"""
}

// Clean reads (adapter and read length filter)
process '1A_clean_reads' {
    tag '1A'
    conda 'bioconda::fastp=0.20.1'
    publishDir outDir, mode: 'copy', pattern: "*.fastp.html"
    input:
        set file(forward_read), file(reverse_read) from pair_reads
    output:
        set file("${forward_read.simpleName}_fastp.fastq.gz"), file("${reverse_read.simpleName}_fastp.fastq.gz") into fastp_2A
        file "${forward_read.simpleName}.fastp.json"
        file ".command.*"
    script:
        """
        fastp -i $forward_read -I $reverse_read -o ${forward_read.simpleName}_fastp.fastq.gz -O ${reverse_read.simpleName}_fastp.fastq.gz \
        --trim_poly_x --length_required 100 --json ${forward_read.simpleName}.fastp.json --html ${forward_read.simpleName}.fastp.html \
        --thread ${threads}
        """
}


// SKA Fastq
process splitKmerReads{
    conda 'bioconda::ska=1.0'
    publishDir outDir + "/${samplename}/ska/", mode: 'copy', pattern: "*.skf"
	input:
	set file(forward_read), file(reverse_read) from fastp_2A

	output:
	file "*.skf" into split_to_align, split_to_compare

	script:
	"""
	ska fastq -k ${kmer} -c ${coverage} -q ${quality_score} -o ${samplename}_${proportion} $forward_read $reverse_read
	"""
}

// SKA Align
process alignSplitFile{
    conda 'bioconda::ska=1.0'
    publishDir outDir + "/${samplename}/ska/", mode: 'copy', pattern: "*.aln"
    input:
    file(split_kmer) from split_to_align

    output:
    file "*.aln" into alignment
    script:
    """
    ska align -v -p ${proportion} ${database} ${split_kmer} -o ${samplename}_${proportion}
    """
}

process iqTree{
    conda 'bioconda::iqtree=2.0.3'
    publishDir outDir + "/${samplename}/iqtree/", mode: 'copy'
    input:
    file(alignment_file) from alignment

    output:
    file "*.contree" into tree_file

    script:
    """
    iqtree -s ${alignment_file} -st DNA -m GTR+G+ASC -nt ${threads} -bb ${bootstrap} -pre ${samplename}_${proportion}
    """
}

process rCode{
    //conda 'python=3.8.5 r::r-base=3.6.1 conda-forge::r-cowplot=1.1.0 //bioconda::bioconductor-ggtree=1.8.2 r::r-ggplot2=3.1.1 bioconda::bioconductor-treeio=1.0.2'

    //r::r-base=3.6.1
    //bioconda::bioconductor-ggtree=1.8.2
    //r::r-ggplot2=3.1.1
    //bioconda::bioconductor-treeio=1.0.2
    //conda-forge::r-cowplot=1.1.0

    conda "${baseDir}/conda/env-R/"
    input:
    file(newick) from tree_file

    script:
    filename = newick.baseName
    pdf_output = outDir + "/${samplename}/${filename}.pdf"
    input = outDir + "/${samplename}/iqtree/${newick}"
    """
    python ${cluster_sample} ${filename} ${r_folder}
    Rscript "${r_folder}treemaker.R" "${r_folder}sample_cluster.txt" ${input} ${pdf_output}
    """
}

process skaCompare{
    conda 'bioconda::ska=1.0'
    publishDir outDir + "/${samplename}/ska/", mode: 'copy', pattern: "*.tsv"
	input:
	file(split_kmer) from split_to_compare

	output:
	file "*.tsv" into split_distances

	script:
	"""
	ska compare ${database} -q ${split_kmer} > ${samplename}_${proportion}.tsv
	"""
}

process jsonify{
    conda 'python=3.8.5 anaconda::pandas=1.1.3'
	input:
	file(distance_dataframe) from split_distances

	script:
	result_folder = outDir + "/${samplename}/"
	"""
	python ${to_json} ${distance_dataframe} ${result_folder} ${parameters}
	"""
}


#!/usr/bin/env nextflow

// (Tool) params
params.reads                = "$baseDir/data/raw_input/20300134001-1_MB_R{1,2}.fastq.gz"
params.outDir               = "$baseDir/output"
params.threads              = 4
// MASH
params.sketchsize                 = "10000"

// Parsing the input parameters
outDir                      = "$params.outDir"
def samplename                = file("$params.reads").simpleName[0].split('_')[0]
threads                     = "$params.threads"
// SKA Fastq
sketchsize                        = "$params.sketchsize"

// Database
//database                    = "$baseDir/db/ref_15mer_${proportion}_variants.skf"
//database                    = "$baseDir/db/Tortolli.msh"
//database                    = "$baseDir/db/Tortolli_10k.msh"
database                    = "$baseDir/db/20230303_database_10K.msh"

// Path prefixes
r_folder                    = "$baseDir/R/"

// Tool paths
cluster_sample              = "$baseDir/python/clustering.py"
create_tree                 = "$baseDir/R/treemaker.R"
to_json                     = "$baseDir/python/finalize.py"
dendropy                    = "$baseDir/python/dendro.py"

parameters = "${sketchsize}"

log.info """
BASICPHYLO V0.2

============INPUT===============
reads      : $params.reads
filename   : $samplename
outDir     : $outDir
~~~~~~~~~~~~parameters~~~~~~~~~~
threads    : $threads
sketchsize : $sketchsize

~~~~~~~~~~~Databases~~~~~~~~~~~~
database    : $database

~~~~~~~~~~~Authors~~~~~~~~~~~~~~
        J.P.M. Coolen
        Pieter Koopman
        Paul Verhoeven (v0.1)
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


// Mash sketch of Fastq reads
process create_mash_sketch_reads{
    conda 'bioconda::mash=2.1'
    publishDir outDir + "/${samplename}/mash/", mode: 'copy', pattern: "*.msh"
	input:
	set file(forward_read), file(reverse_read) from fastp_2A

	output:
	file "${samplename}.msh" into to_create_distance, to_get_stats

	script:
	"""
    #export LD_LIBRARY_PATH="/home/location/gsl/lib"
    cat $forward_read $reverse_read > ${samplename}
    mash sketch -s $sketchsize -r -m 2 ${samplename} -o ${samplename}
	"""
}

// Mash paste and distance
process create_mash_distance_file{
    conda 'bioconda::mash=2.1'
    publishDir outDir + "/${samplename}/mash/", mode: 'copy', pattern: "*.dist"
    input:
    file(sample_mash) from to_create_distance

    output:
    file "*.dist" into distance
    file "total.msh" into total_mash
    script:
    """
    mash paste total ${sample_mash} ${database}
    mash dist -s $sketchsize -t total.msh total.msh > table.dist
    """
}

process create_dendogram_of_distances{
    conda 'bioconda::dendropy=4.5.2'
    publishDir outDir + "/${samplename}/dendropy/", mode: 'copy'
    input:
    file(distance_file) from distance

    output:
    file "*.tre" into tree_file

    script:
    """
    python ${dendropy} --input ${distance_file} --output ${samplename}
    """
}

process rCode{
    //conda 'python=3.8.5 r::r-base=3.6.1 conda-forge::r-cowplot=1.1.0 bioconda::bioconductor-ggtree=1.8.2 r::r-ggplot2=3.1.1 bioconda::biocond$
    //conda 'python r::r-ggplot2 bioconda::bioconductor-ggtree conda-forge::r-cowplot bioconda::bioconductor-treeio'
    //conda install -c r r-base=3.6.1
    //conda install -c bioconda bioconductor-ggtree=1.8.2
    //conda install -c r r-ggplot2=3.1.1
    //conda install -c bioconda bioconductor-treeio=1.0.2-0
    //conda install -c conda-forge r-cowplot=1.1.0

    //conda "${baseDir}/conda/env-41038bd246722f75f387d1ef4a449043/"
    conda "${baseDir}/conda/env-R"
    input:
    file(newick) from tree_file

    script:
    filename = newick.baseName
    pdf_output = outDir + "/${samplename}/${filename}.pdf"
    input = outDir + "/${samplename}/dendropy/${newick}"
    """
    python ${cluster_sample} ${filename} ${r_folder}
    Rscript "${r_folder}treemaker.R" "${r_folder}sample_cluster.txt" ${input} ${pdf_output}
    """
}

process mash_dist_stats{
    conda 'bioconda::mash=2.1'
    publishDir outDir + "/${samplename}/mash/", mode: 'copy', pattern: "*.txt"
	input:
	file(total) from total_mash
    file(sample) from to_get_stats

	output:
	file "*.txt" into split_distances

	script:
	"""
    mash dist -s $sketchsize ${sample} ${total} > stats.txt
	"""
}

//TODO: need to fix and check headers etc.
//process jsonify{
//    conda 'python=3.8.5 anaconda::pandas=1.1.3'
//	input:
//	file(distance_dataframe) from split_distances
//
//	script:
//	result_folder = outDir + "/${samplename}/"
//	"""
//	python ${to_json} ${distance_dataframe} ${result_folder} ${parameters}
//	"""
//}



//input_ch = Channel.fromPath("${params.inputFolder}/*.bam")
covidCsv = file("${params.covid}") 
pancCsv = file("${params.panc}") 
unwrappedFasta = file("${params.unwrappedFasta}")
//hisat2_ch.view()

//hisat2_ch = file("${params.hiSat2Index}/*.ht2")
// input_read_ch = Channel
//     .fromFilePairs("${params.inputFolder}*_R{1,2}*.gz")
//     .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
//     .map { it -> [it[0], it[1][0], it[1][1]]}

//input_read_ch.view()

process RSubset{ 
    container 'quay.io/vpeddu/kim-lab-rotation:latest'
    containerOptions = "--user root"
    publishDir("${params.output}/CsvSubsets")
    beforeScript 'chmod o+rw .'
    cpus 1 

    input:
    file covidCsv
    file pancCsv
    file unwrappedFasta
    output: 
    file "*.fasta" into subsetoutCh
    tuple file('covid_filtered.csv'), file('panc_filtered.csv') into filteredCSVCh
    """
    #!/bin/bash

    Rscript --vanilla ${baseDir}/bin/subset.r ${covidCsv} ${pancCsv}

    echo "grabbing covid fasta"
    for line in `cat covidUnique.txt` 
    do
        grep -A1 \$line ${unwrappedFasta} >> covid_unique_alu_subset.fasta
    done

    echo "grabbing panc fasta"
    for line in `cat pancUnique.txt` 
    do
        grep -A1 \$line ${unwrappedFasta} >> panc_unique_alu_subset.fasta
    done

    """
}

process mafft{ 
    container 'staphb/mafft:7.475'
    containerOptions = "--user root"
    publishDir("${params.output}/mafftOut")

    cpus 4
    memory: 

    input:
    file fasta from subsetoutCh.flatten()
    output: 
    tuple env(base), file("*.aligned.fasta") into mafftOutCh
    """
    #!/bin/bash

    ls -lah

    echo "aligning" 
    
    base=`basename -s "_unique_alu_subset.fasta" ${fasta} `
    echo \$base
    mafft --thread ${task.cpus} ${fasta} > \$base.aligned.fasta
    """
}

process raxml{ 
    container 'staphb/raxml'
    containerOptions = "--user root"
    publishDir("${params.output}/raxml")

    cpus 4
    memory: 

    input:
    tuple val(base), file(alignedFasta) from mafftOutCh
    output: 
    tuple val("${base}"), file("*bestTree*") into raxmlOutCh
    tuple val("${base}"), file("${alignedFasta}") into fastaCollectCh
    """
    #!/bin/bash

    ls -lah

    sed 's/\\:/__/g' ${alignedFasta} > ugh.fasta
    sed "s/'/___/g" ugh.fasta > ughh.fasta
    echo 'tree building'
    #FastTree -gtr -nt < ugh.fasta > ${base}.unique.alu.tree.newick
    raxmlHPC-PTHREADS -m "GTRCATX" -p 420 -T ${task.cpus} -s ughh.fasta -n ${base}.unique.alu.tree.newick
    """
}

process treeFigure{ 
    container 'quay.io/vpeddu/kim-lab-rotation:latest'
    containerOptions = "--user root"
    publishDir("${params.output}/treeFigures")

    cpus 1
    memory: 

    input:
    tuple val(base), file(treeFile) from raxmlOutCh
    tuple file(covidfiltered), file(pancfiltered) from filteredCSVCh
    file unwrappedFasta

    output: 
    tuple val("${base}"), file("*.pdf"), file("${base}.topgenes.fasta") into figureOutCh
    file ('*.rds') into rDataCh
    """
    #!/bin/bash

    Rscript --vanilla ${baseDir}/bin/graph_tree.r ${treeFile} ${base}

    sed 's/__/\\:/g' ${base}.top_genes.txt > ugh.txt

    for line in `cat ugh.txt | cut -f1`
    do
    grep -A1 \$line ${unwrappedFasta} >> ${base}.topgenes.fasta
    done
    """
}

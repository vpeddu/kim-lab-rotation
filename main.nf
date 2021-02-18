
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

    cpus 1 

    input:
    file covidCsv
    file pancCsv
    file unwrappedFasta
    output: 
    //tuple val("${base}"), file("${base}.trimmed_val_1.fq.gz"), file("${base}.trimmed_val_2.fq.gz") into trimOut_ch
    file "*.fasta" into subsetoutCh
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




// process runOpossum {
//     container 'quay.io/vpeddu/u2af1-allele-bias:latest'
//     containerOptions = "--user root"
//     publishDir "${params.output}/Opossum" //, mode: 'symlink'
//     //cpus 1
//     //memory 16.GB
//     input:
//     tuple val(base), file(bam),file(bamIndex) from hiSat2Out_ch
//     file index from fastaIndex_ch
//     file referenceFasta
//     //tuple val(base), file(bam) from samtoolsOut_ch

//     output: 
//     tuple val("${base}"), file("${base}.opossum.bam")  into opossumOut_ch
//     file "${index}" into opossumFastaIndex_ch

//     """
//     #!/bin/bash
//     # logging
//     ls -lah

//     outfilename="${base}"".opossum.bam"
//     echo "output file name is \$outfilename"
    
//     # not sure why this is even necessary but fixed the file not found error
//     touch \$outfilename

//     python2.7 /Opossum-master/Opossum.py --BamFile=${bam} --OutFile="\$outfilename"
//     """
// }


// process samtoolsAddMD {
//     container 'jweinstk/samtools:latest'
//     containerOptions = "--user root"

//     input:
//     tuple val(base), file(bam)  from opossumOut_ch

//     output: 
//     //tuple env(outfilename), file("*.MD.bam"), file("*.fai"), file("*.bai") into samtoolsOut_ch
//     tuple val("${base}"), file ("${bam}"), file("${base}.opossum.bam.bai")  into fileIndex_ch

//     """
//     #!/bin/bash
//     # logging
//     ls -lah

//     touch ${base}.opposum.bam.bai

//     # needed for platypus to run  
//     samtools index -@ ${task.cpus} ${bam} ${base}.opossum.bam.bai
//     """
// }

// process runPlatypus {
//     container 'iarcbioinfo/platypus-nf'
//     containerOptions = "--user root"
//     publishDir "${params.output}/Platypus" , mode: 'copy', overwrite: true
//     //cpus 1
//     //memory 16.GB
//     input:
//     tuple val(base), file(bam), file(bamIndex) from fileIndex_ch
//     //tuple val(base), file(bam) from opossumOut_ch.groupTuple().join(fileIndex_ch )
//     //tuple val(base) file(fastaIndex), file(bamindex) from fileIndex_ch
//     file index from opossumFastaIndex_ch
//     file referenceFasta
 
//     output: 
//     file "${base}.platypus.vcf" 
//     file "${bamIndex}"
//     file "${bam}"

//     """
//     #!/bin/bash
//     # logging
//     ls -lah



//     platypus callVariants \
//         --bamFiles ${bam} \
//         --refFile ${referenceFasta} \
//         --filterDuplicates 0 \
//         --minMapQual 0 \
//         --minFlank 0 \
//         --maxReadLength 500 \
//         --minGoodQualBases 10 \
//         --minBaseQual 20 \
//         -o ${base}.platypus.vcf
//     """
// }


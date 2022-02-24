#!/usr/bin/env nextflow
/*
========================================================================================
                         LAVA
========================================================================================
 Longitudinal Analysis of Viral Alleles
 #### Homepage / Documentation
https://github.com/greninger-lab/lava
----------------------------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
nextflow.preview.dsl=2


def helpMessage() {
    log.info"""
    LAVA: Longitudinal Analysis of Viral Alleles
    Usage:
    An example command for running the pipeline is as follows:
    nextflow run greninger-lab/lava \\
        --RAVA          Run in reference-based mode instead of longitudinal. This mode 
                        compares all samples to the given reference fasta.

        --METADATA      Required argument: A two column csv - the first column is the
                        path to all the fastqs you wish to include in your analysis.
                        All fastqs that you want to include need to be specified in
                        this file AND be located in the folder from which you are
                        running lava. The second column is the temporal seperation
                        between the samples. This is unitless so you can input
                        passage number, days, or whatever condition your experiment
                        happens to have. [REQUIRED]
        
        --OUTDIR        Output directory [REQUIRED]
        
        --FASTA         Specify a reference fasta with the majority consensus of the
                        control fastq. This option must be used with the -g flag to
                        specify the protein annotations relative to the start of this
                        fasta. [REQUIRED IF NOT --GENBANK]
        --GFF           Specify a reference gff file with the protein annotations for
                        the reference fasta supplied with the -f flag. This option
                        must be paired with the -f flag. [REQUIRED IF NOT GENBANK]
        --GENBANK       Provide a Genbank accession number. This record will be used
                        to generate a majority consensus from the control fastq, and
                        this consensus will be annotated from the downloaded genbank
                        record as well. [REQUIRED IF NOT --FASTA + --GFF]
        
        --CONTROL_FASTQ The fastq reads for the first sample in
                        your longitudinal analysis [REQUIRED IF NOT --RAVA]

        --TITLE         Title of your plot. This will appear on the tab when opened with
                        an Internet browser.

               
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

params.OUTDIR = false
params.GENBANK = 'False'
params.GFF = 'False'
params.FASTA = 'NO_FILE'
params.CONTROL_FASTQ = 'NO_FILE'
params.RAVA = 'false'
params.TITLE = 'LAVA Plot'

METADATA_FILE = file(params.METADATA)

// Throws exception if CONTROL_FASTQ doesn't exist 
if(params.RAVA=='false') {
    CONTROL_FASTQ = file(params.CONTROL_FASTQ, checkIfExists:true)
}

/*
 * Import the processes used in this workflow
 */

include { CreateGFF_Genbank_RAVA }from './Modules.nf'
include { CreateGFF } from './Modules.nf'
include { Alignment_prep } from './Modules.nf'
include { Align_samples } from './Modules.nf' 
include { Process_samples } from './Modules.nf' 
include { Pipeline_prep } from './Modules.nf'
// include { Create_VCF } from './Modules.nf'
include { Extract_variants } from './Modules.nf'
// include { Annotate_complex } from './Modules.nf'
include { Generate_output } from './Modules.nf'

/*
 * Staging python scripts
 */
PULL_ENTREZ = file("$workflow.projectDir/bin/pull_entrez.py")
WRITE_GFF = file("$workflow.projectDir/bin/write_gff.py")
INITIALIZE_PROTEINS_CSV = file("$workflow.projectDir/bin/initialize_proteins_csv.py")
ANNOTATE_COMPLEX_MUTATIONS = file("$workflow.projectDir/bin/Annotate_complex_mutations.py")
MAT_PEPTIDE_ADDITION = file("$workflow.projectDir/bin/mat_peptide_addition.py")
// RIBOSOMAL_SLIPPAGE = file("$workflow.projectDir/bin/ribosomal_slippage.py")
GENOME_PROTEIN_PLOTS = file("$workflow.projectDir/bin/genome_protein_plots.py")
PALETTE = file("$workflow.projectDir/bin/palette.py")
VCFUTILS = file("$workflow.projectDir/bin/vcfutils.pl")
FIX_VARIANTS_FILE = file("$workflow.projectDir/bin/fix_variants_file.py")

/*
 * Error handling for input flags
 */

//if OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}

// Make sure control fastq is set if not --RAVA
if (!params.CONTROL_FASTQ && params.RAVA=='false'){
    println("Must provide control FASTQ with --CONTROL_FASTQ") 
    exit(1)
}

// Make sure one of the two reference options is correctly set
if(params.GFF == "False" && params.FASTA == 'false' && params.GENBANK == "False"){ 
    println('Either --GENBANK or --FASTA + --GFF are required flags')
    exit(1)
}

// Parsing metadata file
input_read_ch = Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{ row-> tuple(file(row.Sample), (row.Passage)) }

// Throws exception if paths in METADATA are not valid
Channel
    .fromPath(METADATA_FILE)
    .splitCsv(header:true)
    .map{row-> (file(row.Sample, checkIfExists:true))}.ifEmpty{error "Check metadata file"}
    //.map{row-> (file(row.Sample).isEmpty())}
    //.filter{ it == false}.subscribe{println it}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    log.info nfcoreHeader()

    // If we are running with --RAVA
    if(params.RAVA != 'false') {
        if(params.FASTA == "NO_FILE") {
        CreateGff_Genbank_RAVA ( 
            params.GENBANK,
            PULL_ENTREZ,
            WRITE_GFF
        )
        
        Alignment_prep ( 
            CreateGff_Genbank_RAVA.out[0],
            CreateGff_Genbank_RAVA.out[1],
            CreateGff_Genbank_RAVA.out[2],
            CreateGff_Genbank_RAVA.out[3],
            CreateGff_Genbank_RAVA.out[4]
        )   
    }
    else {
        CreateGFF ( 
            params.GENBANK,
            PULL_ENTREZ,
            WRITE_GFF,
            file(params.FASTA),
            file(params.GFF)
        )
        
        Alignment_prep ( 
            CreateGFF.out[0],
            CreateGFF.out[1],
            CreateGFF.out[2],
            CreateGFF.out[3],
            CreateGFF.out[4]
        )   
    }

        Align_samples ( 
            input_read_ch,
            Alignment_prep.out[0],            
        )

        Process_samples (
            Align_samples.out[0],
            VCFUTILS,
            Alignment_prep.out[1],
            Alignment_prep.out[2],
            Alignment_prep.out[0]
        )

        Pipeline_prep ( 
            Align_samples.out[0].collect(),
            Alignment_prep.out[3],
            Alignment_prep.out[0],
            INITIALIZE_PROTEINS_CSV
        )

        Extract_variants ( 
            Process_samples.out[3],
            METADATA_FILE,
            Alignment_prep.out[2],
            Alignment_prep.out[4],
            Pipeline_prep.out[1],
            FIX_VARIANTS_FILE,
            Alignment_prep.out[5],
            MAT_PEPTIDE_ADDITION
        )

        Generate_output( 
            Extract_variants.out[2].collect(),
            Extract_variants.out[3].collect(),
            Extract_variants.out[4].collect(),
            Pipeline_prep.out[0],
            Pipeline_prep.out[1],
            Process_samples.out[1].collect(),
            Process_samples.out[4].collect(),
            GENOME_PROTEIN_PLOTS,
            PALETTE,
            Alignment_prep.out[5],
            MAT_PEPTIDE_ADDITION,
            METADATA_FILE,
            params.TITLE
        )
    }
    // If we are running on LAVA mode
    else {
        if(params.FASTA == 'false') {
            CreateGFF_Genbank ( 
                params.GENBANK, 
                CONTROL_FASTQ,
                PULL_ENTREZ,
                WRITE_GFF
            )
            Alignment_prep ( 
                CreateGFF_Genbank.out[0],
                CreateGFF_Genbank.out[1],
                CreateGFF_Genbank.out[2],
                CreateGFF_Genbank.out[3],
                CreateGFF_Genbank.out[4],
                CreateGFF_Genbank.out[5]
            )

        }
        else {
            if(params.FASTA == 'false') {
                CreateGFF ( 
                params.GENBANK, 
                CONTROL_FASTQ,
                PULL_ENTREZ,
                WRITE_GFF,
                CreateGFF_Genbank.out[6],
                CreateGFF_Genbank.out[7]
                )
            } else {
                CreateGFF ( 
                params.GENBANK, 
                CONTROL_FASTQ,
                PULL_ENTREZ,
                WRITE_GFF,
                file(params.FASTA),
                file(params.GFF)
                )
            }
        }
    }
}

def nfcoreHeader() {

    return """
                       ooO
                     ooOOOo
                   oOOOOOOoooo
                 ooOOOooo  oooo
                /vvv\\
               /V V V\\ 
              /V  V  V\\          
             /         \\            oh wow  look at these alleles
            /           \\          /         /    
          /               \\   	  o          o
__       /                 \\     /-   o     /-
/\\     /                     \\  /\\  -/-    /\\
                                    /\\
 ___      _______  __   __  _______ 
|   |    |   _   ||  | |  ||   _   |
|   |    |  |_|  ||  |_|  ||  |_|  |
|   |    |       ||       ||       |
|   |___ |       ||       ||       |
|       ||   _   | |     | |   _   |
|_______||__| |__|  |___|  |__| |__| Version 2    
    """.stripIndent()
}

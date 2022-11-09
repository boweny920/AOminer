# The gff to be used must be GFF3, format explaination: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

import argparse
from Annotation import annotation_prodigal
import os
import sys
from Bio import SeqIO
from AOminer_process import AO_Find_process

def is_non_zero_file(fpath):
    #Check if file is empty or does not exsit
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def make_output_folder(GCF_Name):
    os.makedirs(GCF_Name)
    return os.path.dirname(GCF_Name)

parser = argparse.ArgumentParser(description="AOminer")
parser.add_argument('-n','--FNA_file',nargs='?',default=None,help="Genome FNA file")
parser.add_argument('-p','--FAA_file', nargs='?',default=None ,help="Genome annotated FAA file")
parser.add_argument('-g','--GFF_file', nargs='?',default=None, help="Genome annotated GFF file")
parser.add_argument('-j','--just_operon',action='store_true',help='Provide option -j/--just_operon if input data is the faa file of "ONE" operon')
parser.add_argument('-m','--mode_prodiagal', nargs='?', default="meta",choices=["single","meta"],help="mode prodigal will be run")
parser.add_argument('-o','--outputFolder',type=str,default="AOminer_Output",help="folder containing all output results")
parser.add_argument('-l','--all_protein_length_in_AcrAca_operon', nargs='?', default=600,type=int,help="max proten length in Acr-Aca operon when length of Acr homolog < 200aa")
parser.add_argument('-i','--intergenic_dist_in_AcrAca_operon',nargs='?',default=250,type=int,help="Maximum Intergenic distance in Acr-Aca operon")
parser.add_argument('-w', '--Prok', action='store_true',help="Provide option -w/--Prok if input data is a Prokaryotic assembled genome/contig")
parser.add_argument('-d','--threads',nargs='?', default="1", help="Number of cpu cores to run the program with")
parser.add_argument('-z','--phamDir', nargs='?', default="all_pFam_hmm", help="Directory of all pfam hmm files with .dat files and other binaries")
parser.add_argument('-r','--Acr_protein_database',nargs='?',default="publishedACRs.faa",help="The Acr proteins published")
parser.add_argument('-a','--HTH_HMM',nargs='?',default="HTH_hmms/HTH_HMM_strict",help="Curated HTH database")
parser.add_argument('-e','--PfamScan_evalue',nargs='?', default="1e-1", help="Evalue cut-off for Pfam annotation")
parser.add_argument('-x','--execute_level', nargs='?', choices=["strict","medium","relaxed"],default="strict",help="Selection levels to predicted AOs, recommeded strict")
args=parser.parse_args()
HMMdb=os.path.join("dbPFhmm","TrainPresent_AOPF_nonAOPF_HMMs_ERROR_CORRECT.hmm")

if os.path.isdir(args.outputFolder) is not True:
    os.makedirs(args.outputFolder)

if args.just_operon:
    from AOminer_process_OPERON import AO_Find_process
    if args.FAA_file != None:
        if len([v.id for v in SeqIO.parse(args.FAA_file,"fasta")]) >1:
            AO_Find_process(args.FAA_file,
                            args.outputFolder,
                            args.all_protein_length_in_AcrAca_operon,
                            args.intergenic_dist_in_AcrAca_operon,
                            args.Acr_protein_database,
                            args.Prok,
                            args.threads, args.phamDir,
                            HMMdb,
                            args.execute_level,
                            args.PfamScan_evalue,
                            args.HTH_HMM).run_process()
        else: print("Input operon must contain at least 2 proteins...")
    else: print("Please provide the input operon protein sequences in fasta format...")

else:
    if (args.FAA_file == None or args.GFF_file == None) and args.FNA_file != None:
        print("No annotations provided, using prodigal to annotate genome :)")
        prodigal_outDir=os.path.join(args.outputFolder,os.path.basename(args.FNA_file)+".prodigalOUT")
        gff,faa,fna = annotation_prodigal(args.FNA_file, args.mode_prodiagal,prodigal_outDir).run_prodigal()
        Prodigal=True

    elif args.FAA_file != None and args.GFF_file != None and args.FNA_file != None:
        gff=args.GFF_file
        faa=args.FAA_file
        Prodigal = False
        with open(os.path.join(args.outputFolder, os.path.basename(args.FNA_file)), "w") as newfile:
            for record in SeqIO.parse(args.FNA_file, "fasta"):
                record.description = ""
                SeqIO.write(record, newfile, "fasta")
        fna=os.path.join(args.outputFolder, os.path.basename(args.FNA_file))


    elif (args.FAA_file is None or args.GFF_file is None) and args.FNA_file == None:
        sys.exit("##No FNA, GFF and FAA detected##\nPlease provide FNA file or FNA with annotated protein FAA file and associated GFF file.")

    # Need the model running script

    AO_Find_process(fna,gff,faa,
                    args.outputFolder,
                    args.all_protein_length_in_AcrAca_operon,
                    args.intergenic_dist_in_AcrAca_operon,
                    args.Acr_protein_database,
                    args.Prok,
                    args.threads,args.phamDir,
                    HMMdb,
                    args.execute_level,
                    args.PfamScan_evalue,
                    args.HTH_HMM,
                    Prodigal).run_process()


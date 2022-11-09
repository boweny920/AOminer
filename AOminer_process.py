import os
import pandas as pd
from functions import *
from Bio import SeqIO
import numpy as np

class AO_Find_process:
    def __init__(self,fna_file,gff_file,faa_file,outputFolder,all_protein_length_in_AcrAca_operon,intergenic_dist_in_AcrAca_operon,KnowAcrFaa,Prok,threads,phamDir,hmm_allPF,execute_level,PfamScan_evalue,HTH_HMM,Prodigal):
        self.fna=fna_file
        self.faa=faa_file
        self.gff=gff_file
        self.outputFolder=outputFolder
        self.all_protein_length_in_AcrAca_operon=all_protein_length_in_AcrAca_operon
        self.intergenic_dist_in_AcrAca_operon=intergenic_dist_in_AcrAca_operon
        self.KnowAcrFaa= KnowAcrFaa
        self.Prok=Prok
        self.threads=threads
        self.phamDir=phamDir
        self.hmm_allPF=hmm_allPF
        self.execute_level=execute_level
        self.PfamScan_evalue=PfamScan_evalue
        self.HTH_HMM=HTH_HMM
        self.Prodigal=Prodigal

    def run_process(self):
        sub_outputfolder_path=self.outputFolder
        file_list=[]
        file=subprocess.Popen(['grep "CDS" %s|grep -v "#"'%(self.gff)],shell=True, stdout=subprocess.PIPE)
        for line in file.stdout:
            line = line.decode('utf-8').rstrip().split('\t')
            file_list.append(line)
        file_dic={}
        for i in enumerate(file_list):
            file_dic.setdefault(i[0],i[1])

        ##Search for Acr homologs##
        diamond_outfile_Acr=run_diamond(self.faa,sub_outputfolder_path,self.KnowAcrFaa,self.threads)
        if is_non_zero_file(diamond_outfile_Acr) is not False:
            print("Known Acr homologs found in input data","...")
            protein_NP_list, dic_Acr = parse_diamond_get_proteinID(diamond_outfile_Acr)
            faa_file_wrote(protein_NP_list, self.faa, os.path.join(sub_outputfolder_path, "Acr_homologs.faa"))
        else:
            protein_NP_list=[]
            dic_Acr={}
        ##Search for HTH homologs##
        HMM_outfile_HTH=run_hmmscan(self.faa,'1e-3',self.HTH_HMM,self.threads, sub_outputfolder_path,"HTH_search")
        if is_non_zero_file(HMM_outfile_HTH) is not False:
            dic_hth = parse_hmmOutfile(HMM_outfile_HTH,"0.8")
        else:
            dic_hth={}

        loci_list = loci_select(file_dic, dic_Acr, self.KnowAcrFaa, self.Prodigal ,self.all_protein_length_in_AcrAca_operon, self.intergenic_dist_in_AcrAca_operon)
        if len(loci_list) > 0:
            print("Short gene operons found in input sequences", "...")
            # Put the new faa files in a directory called "Short_gene_operons"
            output_1_DirPath = os.path.join(sub_outputfolder_path, "Short_Gene_Operons")
            os.makedirs(output_1_DirPath)
            # Write each locus into individual fasta files
            AcrOperons_predicted = []
            for order, locus in enumerate(loci_list):
                SGO_filename = os.path.join(output_1_DirPath,"SGO_OperonNumber-" + str(order) + ".faa")
                faa_file_wrote([v[1] for v in locus["operon"]], self.faa, SGO_filename)
                ##run hmmscan on this SGO to label the operon with AOPF + nonAOPFs
                dic_PF_proID = parse_hmmOutfile(run_hmmscan(SGO_filename, "1e-3", self.hmm_allPF, self.threads, output_1_DirPath),"0.5") #evalue=1e-3, coverage=0.5
                locus_PF=data_list_gen(dic_PF_proID,SGO_filename)
                locus_PF_AOprob,is_AO=HMM_output(locus_PF,self.execute_level)
                if is_AO:
                    locus["HMM_Prob"]=locus_PF_AOprob
                    locus["Operon_num"]=order
                    locus["SGO_filename"]=SGO_filename
                    locus["AOPF_list"]=locus_PF
                    AcrOperons_predicted.append(locus)
                    print("Predicted anti-CRISPR operon protein sequences saved in %s" % SGO_filename)
            ##Check to see if CCtyper and VIBRANT need to run
            if len(AcrOperons_predicted) >0:
                # Search for prophage and CRISPR-Cas
                if self.Prok == True:
                    ##check if any complete CRISPR-Cas found, and check for self-targeting locations if any
                    complete_CRISPR_Cas_systems = find_complete_CRISPR_Cas_and_SelfTargeting(self.fna, os.path.join(sub_outputfolder_path, "CRISPR_Cas_Found"), self.threads,sub_outputfolder_path)  # Name of cctyper output directory is "CRISPR_Cas_Found"
                    # [Contig|CasTyper|Position|self-targeting regions,...]/None
                    prophage_regions = find_prophage(self.fna, sub_outputfolder_path, self.threads)
                    # [Contig:startPos-endPos,...]/None
                else:
                    complete_CRISPR_Cas_systems = None
                    prophage_regions = None

                ###Final output dataframe make###
                protein_faa_dic = SeqIO.to_dict(SeqIO.parse(self.faa, "fasta"))
                fna_dic = SeqIO.to_dict(SeqIO.parse(self.fna, "fasta"))
                df_allResult = pd.DataFrame(
                            columns=["Operon Number", "AO Score","Protein ID", "Contig ID", "Strand", "Protein Length (nt)",
                                     "Start", "End", "Acr Homolog", "With HTH", "Pfam Annotation",
                                     "Complete CRISPR-Cas and STSS", "Operon in Prophage", "Protein Sequence"])
                for AO in AcrOperons_predicted:
                    operon_number_list = ["SGO#"+str(AO["Operon_num"])]*len(AO["check_result"])
                    AO_scor_list= [AO["HMM_Prob"]]*len(AO["check_result"])
                    AO_npArray=np.array(AO["check_result"],dtype=object)
                    ##Acr/Aca Homolog Info##
                    AO_acr_homologs_list=[]
                    for v in AO["operon"]:
                        if v[1] in dic_Acr:
                            AO_acr_homologs_list.append(dic_Acr[v[1]])
                        else: AO_acr_homologs_list.append(np.nan)
                    AO_hth_homologs_list = []
                    for v in AO["operon"]:
                        if v[1] in dic_hth:
                            AO_hth_homologs_list.append(dic_hth[v[1]])
                        else: AO_hth_homologs_list.append(np.nan)
                    ###
                    df_AO=pd.DataFrame(AO_npArray,columns=["Protein ID", "Contig ID", "Strand", "Protein Length (nt)","Start", "End"])
                    ###
                    df_AO["AO Score"]=AO_scor_list
                    df_AO["Operon Number"]=operon_number_list
                    df_AO["Acr Homolog"]=AO_acr_homologs_list
                    df_AO["With HTH"]=AO_hth_homologs_list
                    ##Add prophage & CRISPR into df##
                    if prophage_regions is not None:
                        prophage_containing_AcaOperon_list = [prophage_harboring_operonFind(df_outTable, prophage_regions)] * len(AO["check_result"])  # Contig:startPos-endPos/None
                    else: prophage_containing_AcaOperon_list=[np.nan] * len(AO["check_result"])
                    df_AO["Operon in Prophage"]=prophage_containing_AcaOperon_list
                    if complete_CRISPR_Cas_systems is not None:
                        Complete_CRISPR_Cas_inContig_list=[";".join(complete_CRISPR_Cas_systems)] * len(AO["check_result"])
                    else: Complete_CRISPR_Cas_inContig_list=[np.nan] * len(AO["check_result"])
                    df_AO["Complete CRISPR-Cas and STSS"]=Complete_CRISPR_Cas_inContig_list
                    ###
                    operon_faa_list=[]
                    for protein in AO["operon"]:
                        if "pseudo" in protein[1]:
                            operon_faa_list.append("pseudo")
                        else:
                            operon_faa_list.append(str(protein_faa_dic[protein[1]].seq))
                    df_AO["Protein Sequence"]=operon_faa_list
                    ### Pfam annotation
                    dic_pfam = pfamScan_run(AO["SGO_filename"], self.phamDir, self.threads, self.PfamScan_evalue)
                    #dic_pfam={}
                    pfam_list=[]
                    for protein in AO["operon"]:
                        if protein[1] in dic_pfam:
                            pfam_list.append(dic_pfam[protein[1]])
                        else: pfam_list.append(np.nan)
                    df_AO["Pfam Annotation"]=pfam_list
                    ###
                    df_allResult=pd.concat([df_allResult,df_AO],ignore_index=True)
                df_allResult.to_csv(os.path.join(sub_outputfolder_path, "ALL_AOs_predicted.csv"), index=False)

            else:
                print("No AOs were predicted from input contig/genome :(")
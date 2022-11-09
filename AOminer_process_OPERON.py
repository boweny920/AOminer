import os
import pandas as pd
from functions import *
from Bio import SeqIO
import numpy as np

class AO_Find_process:
    def __init__(self, faa_file, outputFolder, all_protein_length_in_AcrAca_operon,
                 intergenic_dist_in_AcrAca_operon, KnowAcrFaa, Prok, threads, phamDir, hmm_allPF, execute_level,
                 PfamScan_evalue, HTH_HMM):
        self.faa = faa_file
        self.outputFolder = outputFolder
        self.all_protein_length_in_AcrAca_operon = all_protein_length_in_AcrAca_operon
        self.intergenic_dist_in_AcrAca_operon = intergenic_dist_in_AcrAca_operon
        self.KnowAcrFaa = KnowAcrFaa
        self.Prok = Prok
        self.threads = threads
        self.phamDir = phamDir
        self.hmm_allPF = hmm_allPF
        self.execute_level = execute_level
        self.PfamScan_evalue = PfamScan_evalue
        self.HTH_HMM=HTH_HMM

    def run_process(self):
        sub_outputfolder_path = self.outputFolder

        ##Search for Acr homologs##
        diamond_outfile_Acr = run_diamond(self.faa, sub_outputfolder_path, self.KnowAcrFaa, self.threads)
        if is_non_zero_file(diamond_outfile_Acr) is not False:
            print("Known Acr homologs found in input data", "...")
            protein_NP_list, dic_Acr = parse_diamond_get_proteinID(diamond_outfile_Acr)
            faa_file_wrote(protein_NP_list, self.faa, os.path.join(sub_outputfolder_path, "Acr_homologs.faa"))
        else:
            protein_NP_list = []
            dic_Acr = {}
            ##Search for HTH homologs##
        HMM_outfile_HTH = run_hmmscan(self.faa, '1e-3', self.HTH_HMM, self.threads, sub_outputfolder_path,
                                          "HTH_search")
        if is_non_zero_file(HMM_outfile_HTH) is not False:
            dic_hth = parse_hmmOutfile(HMM_outfile_HTH, "0.8")
        else:
            dic_hth = {}

        output_1_DirPath=os.path.join(sub_outputfolder_path,"hmmer_results")
        if os.path.isdir(output_1_DirPath) == False: os.makedirs(output_1_DirPath)
        dic_PF_proID = parse_hmmOutfile(run_hmmscan(self.faa, "1e-3", self.hmm_allPF, self.threads, output_1_DirPath),"0.5")  # evalue=1e-3, coverage=0.5
        locus_PF = data_list_gen(dic_PF_proID, self.faa)
        locus_PF_AOprob, is_AO = HMM_output(locus_PF, self.execute_level)
        if is_AO:
            print("Input operon is predicted to be anti-CRISPR operon")
            dic_faa=SeqIO.to_dict(SeqIO.parse(self.faa,"fasta"))
            list_faa=[record.id for record in SeqIO.parse(self.faa,"fasta")]
            df_AO=pd.DataFrame()
            df_AO["AO Score"]=[locus_PF_AOprob] * len(locus_PF)
            df_AO["Protein ID"]=list_faa
            ##Acr/Aca Homolog Info##
            AO_acr_homologs_list = []
            for v in list_faa:
                if v in dic_Acr:
                    AO_acr_homologs_list.append(dic_Acr[v])
                else:
                    AO_acr_homologs_list.append(np.nan)
            AO_hth_homologs_list = []
            for v in list_faa:
                if v in dic_hth:
                    AO_hth_homologs_list.append(dic_hth[v])
                else:
                    AO_hth_homologs_list.append(np.nan)
            df_AO["Acr Homolog"] = AO_acr_homologs_list
            df_AO["With HTH"] = AO_hth_homologs_list
            ###
            df_AO["Protein Length (nt)"]=[len(record.seq)*3 for record in SeqIO.parse(self.faa,"fasta")]
            df_AO["Protein Sequence"]=[str(record.seq) for record in SeqIO.parse(self.faa,"fasta")]

            ### Pfam annotation
            dic_pfam = pfamScan_run(self.faa, self.phamDir, self.threads, self.PfamScan_evalue)
            #dic_pfam={}
            pfam_list = []
            for protein in list_faa:
                if protein in dic_pfam:
                    pfam_list.append(dic_pfam[protein])
                else:
                    pfam_list.append(np.nan)
            df_AO["Pfam Annotation"] = pfam_list
            df_AO=df_AO[["AO Score","Protein ID", "Protein Length (nt)",
                                    "Acr Homolog", "With HTH", "Pfam Annotation",
                                     "Protein Sequence"]]
            df_AO.to_csv(os.path.join(sub_outputfolder_path, "AO_predicted.csv"), index=False)
        else:
            print("Input operon was not predicted as AO :(")





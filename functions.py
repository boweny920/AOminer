import subprocess
import sys
import pandas as pd
from Bio import SeqIO
import os
import pprint as pp
import numpy as np
import math
import json

def make_output_folder(GCF_Name):
    os.makedirs(GCF_Name)
    return os.path.dirname(GCF_Name)

def get_file_directory_path(file_path):
    directory_path=str(file_path).rsplit("/",1)
    return directory_path[0]

def length(list):
    gene_length=int(list[4])-int(list[3]) + 1
    return gene_length

def contig(list):
    contig=list[0]
    return contig

def strand(list):
    gene_strand=list[6]
    return gene_strand

def proteinInfo(list,Prodigal):
    if Prodigal:
        return  [list[0], list[-1]]
    else:
        return list[-1]

def proteinID(IDinfo,Prodigal):
    if Prodigal:
        IDinfo_P = IDinfo[-1].split(";")
        for v in IDinfo_P:
            if "ID=" in v:
                ID = str(v).split("_", 1)
                ProID = str(IDinfo[0]) + "_" + ID[-1]
                return ProID
    else:
        IDinfo = IDinfo.split(";")
        for v in IDinfo:
            if "protein_id=" in v:
                ProID = str(v).strip('protein_id=')
                return ProID

def protein_pos_start(list):
    protein_start=int(list[3])
    return protein_start

def protein_pos_end(list):
    protein_end = int(list[4])
    return protein_end

def proteinInfo_List_process(proteinInfo,Prodigal):
    loci_list_proteins_WithPseudo = []
    for number, value in enumerate(proteinInfo):
        if "pseudo=true" not in value:
            loci_list_proteins_WithPseudo.append([int(number), proteinID(value,Prodigal)])
        elif "pseudo=true" in value:
            loci_list_proteins_WithPseudo.append([int(number), "pseudo"])
    pseudo_gene_number = 0
    for x in loci_list_proteins_WithPseudo: pseudo_gene_number = pseudo_gene_number + x.count("pseudo")
    if len(loci_list_proteins_WithPseudo) - pseudo_gene_number >= 2:
        ##see how many genes left without the pseudo gene
        return loci_list_proteins_WithPseudo

def result_check_list(list,Prodigal):
    gene_strand=strand(list)
    gene_length=length(list)
    if "pseudo=true" not in list[-1]:ProID=proteinID(proteinInfo(list,Prodigal),Prodigal)
    elif "pseudo=true" in list[-1]:ProID="pseudo"
    # output result in format ProteinID; contig;gene strand; gene length；gene start;gene end, in a list
    check_result = [ProID, contig(list), gene_strand, gene_length, list[3], list[4]]
    return check_result

def loci_select(file_dic, dic_Acr, knwonAcrfaa, Prodigal ,all_protein_lenbp=600, intergenic_dist=250):
    ###################################################################
    ###   Get the regions that contain at least two genes, each gene must
    ###   be on the same strand, each gene must be less than 200 aa or 600 nucleotides
    key = 0
    lst = []
    loci_list = []
    while key in file_dic.keys():
        # First get operons without considering protein length
        if key - 1 in file_dic.keys() and strand(file_dic[key]) == strand(file_dic[key - 1]) and contig(file_dic[key]) == contig(file_dic[key - 1]) and int(protein_pos_start(file_dic[key]) - protein_pos_end(file_dic[key - 1])) < intergenic_dist:
            lst.append(file_dic[key])
            key = key + 1
        elif key + 1 in file_dic.keys() and strand(file_dic[key]) == strand(file_dic[key + 1]) and contig(
                file_dic[key]) == contig(file_dic[key + 1]) and int(protein_pos_start(file_dic[key + 1]) - protein_pos_end(file_dic[key])) < intergenic_dist:
            if len(lst) >= 2:
                # if previous lst has 2 or more genes, append the lst to the loci_list
                loci_list.append(lst)
            lst = []
            lst.append(file_dic[key])
            lst.append(file_dic[key + 1])
            key = key + 2
        else:
            key = key + 1
            if len(lst) >= 2:
                loci_list.append(lst)
            lst = []
    else:
        if len(lst) >= 2:
            loci_list.append(lst)
    knownAcrFaa_dic=SeqIO.to_dict(SeqIO.parse(knwonAcrfaa,"fasta"))
    loci_list_length_adapat=[]
    for operon in loci_list:
        operon_numbered_proIDs=proteinInfo_List_process([proteinInfo(v,Prodigal) for v in operon],Prodigal)
        Acrs_in_operon = [dic_Acr[proteinID(proteinInfo(v,Prodigal),Prodigal)] for v in operon if proteinID(proteinInfo(v,Prodigal),Prodigal) in dic_Acr]
        if len(Acrs_in_operon) > 0:
            Acrs_in_operon_maxLength=max([len(knownAcrFaa_dic[v].seq) for v in Acrs_in_operon])
            if Acrs_in_operon_maxLength * 3 < 600:
                maxProtein_length = all_protein_lenbp
            elif Acrs_in_operon_maxLength * 3 > 600:
                maxProtein_length = Acrs_in_operon_maxLength * 3
        else: maxProtein_length=all_protein_lenbp
        if any(v for v in operon if length(v) > maxProtein_length) == False and operon_numbered_proIDs != None:
            # A list of dictionaries, with operon and check_result as keys
            loci_list_length_adapat.append({"operon":operon_numbered_proIDs,"check_result":[result_check_list(v,Prodigal) for v in operon]})
    return loci_list_length_adapat

def is_non_zero_file(fpath):
    #Check if file is empty or does not exsit
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def make_and_check_output_directory(dpath):
    if os.path.isdir(dpath) is False:
        os.makedirs(dpath)
    return str(dpath)

def run_diamond(database_file,outputdir,query,threads,evalue_cutoff="1e-3",coverage_cut_off="0.7"): #add coverage cutoff
    ##Use this as dmond_out
    outputdir=os.path.join(outputdir,"Homology_Search_Results")
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)
    subprocess.Popen(['diamond','makedb','--in',database_file,'-d',os.path.join(outputdir,"diamond_db")]).wait()
    subprocess.Popen(['diamond','blastp','-d',os.path.join(outputdir,"diamond_db"),'-f','6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen' ,'qstart' ,'qend' ,'sstart' ,'send' ,'evalue' ,'bitscore', 'qlen', 'slen','--quiet','--more-sensitive','-q',query,'-p',threads,'-o',os.path.join(outputdir,"diamond_blastp_%s.txt"%(os.path.basename(query.replace(".faa","")))),'-e',evalue_cutoff]).wait()
    ## Output format 6 columns: Default: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    subprocess.Popen(["awk -F '\t' '($8-$7)/$(NF-1) > %s {print}' %s > %s"%(str(coverage_cut_off),os.path.join(outputdir,"diamond_blastp_%s.txt"%(os.path.basename(query.replace(".faa","")))),os.path.join(outputdir,"diamond_blastp_%s.coverageParsed.txt"%(os.path.basename(query.replace(".faa","")))))],shell=True).wait()
    return os.path.join(outputdir,"diamond_blastp_%s.coverageParsed.txt"%(os.path.basename(query.replace(".faa",""))))

def parse_diamond_get_proteinID(dmond_out):
    ##Use with run_diamond; Use on NON-EMPTY diamond output file
    ##Return also a dict of Protein ID to Acr type
    dic_acr={}
    protein_ID_list=[]
    Pro_ID = subprocess.Popen("awk -F '\t' '{print $1,$2}' %s" % (dmond_out), shell=True, stdout=subprocess.PIPE)
    for line in Pro_ID.stdout:
        line = line.rstrip().decode('utf-8').split()
        proID = str(line[1])
        protein_ID_list.append(proID)
        dic_acr.setdefault(proID,str(line[0]))
    return protein_ID_list,dic_acr

def faa_file_wrote(one_acr_aca_locus,faa_file,newfile_name_dirctory):
    ##The one_acr_aca_locus is a list of proteinID.
    ##The function takes the IDs and the faa file with protein sequences
    ##And make a fasta newfile containing ProteinIDs with corresponding amino acid sequences
    record_dict = SeqIO.to_dict(SeqIO.parse(faa_file, "fasta"))
    with open(newfile_name_dirctory,"w") as newfile:
        for proID in one_acr_aca_locus:
            SeqIO.write(record_dict[proID],newfile,"fasta")

def run_hmmscan(faafile,evalue_cut_off,hmmfile,threads,output_dir,hth=''):
    out_dir=os.path.join(output_dir,"Homology_Search_Results")
    if os.path.isdir(out_dir) == False: os.makedirs(out_dir)
    hmm_outfile=os.path.join(out_dir,os.path.basename(faafile)+".hmmout"+"."+hth)
    subprocess.Popen(['hmmscan','--domtblout',hmm_outfile,'--noali',"-o","log.hmm","--cpu",threads,'-E', evalue_cut_off,hmmfile,faafile]).wait()
    return hmm_outfile

def parse_hmmOutfile(hmm_outfile,hmm_coverage_cutoff): # redo the coverage calculation part!
    parsed_hmm_outfile = str(hmm_outfile) + ".Coverage_parsed"
    #Parse the hmm output in terms of subject coverage: subjust aligned length/ Total subject length
    subprocess.Popen(["grep -v '#' %s | awk '($17-$16)/$3 > %s {print $1,$2,$4, ($17-$16)/$3}'  > %s" % (hmm_outfile, str(hmm_coverage_cutoff),parsed_hmm_outfile)],shell=True).wait()
    dic_AOPF={}
    for line in open(parsed_hmm_outfile).readlines():
        AOPF, pfamID, ID, coverage = line.rstrip().split()
        if ID in dic_AOPF:
            dic_AOPF[ID][float(coverage)]=AOPF
        else:
            dic_AOPF[ID]={float(coverage):AOPF}
    dic_AOPF_bestSelected={}
    for id in dic_AOPF:
        dic_AOPF_bestSelected.setdefault(id,dic_AOPF[id][max(v for v in dic_AOPF[id])])
    return dic_AOPF_bestSelected

def data_list_gen(dic_PF_proID,SGO_FAA):
    print_list=[]
    for record in SeqIO.parse(SGO_FAA, "fasta"):
        if record.id in dic_PF_proID:
            print_list.append(dic_PF_proID[record.id])  # AO-PF annotation
        else:
            print_list.append("NoAO")
    return print_list

def final_result_check_output_generation(Acr_homolog_candidate_list_resultCheck,newfile_name,Aca_Pro_lst_dic,publishedAcaHMM_hits_dic):
    #generation of final result checking file
    #Table will be tsv format
    #### GCF_number    ProteinID   Contig  Strand  Length  Start   End     Species     Acr     Aca  AcaHMM_HIT ####
    Acr_Aca_loci_by_GBA = []
    for i in Acr_homolog_candidate_list_resultCheck:
        if any(v for v in i if v[0] in Aca_Pro_lst_dic):
            for v in i:
                # Check which ProteinID is the Aca homolog
                if v[0] in Aca_Pro_lst_dic:
                    v.append("Aca_protein|"+";".join(Aca_Pro_lst_dic[v[0]]))
                else:
                    v.append("NA")
                # Check which ProteinID has AcaHMM hits
                if v[0] in publishedAcaHMM_hits_dic:
                    v.append(";".join(publishedAcaHMM_hits_dic[v[0]]))
                else:
                    v.append("NA")
            Acr_Aca_loci_by_GBA.append(i)

    with open(newfile_name, "w") as newfile:
        for pro in Acr_Aca_loci_by_GBA:
            for value in pro:
                # for info in value: newfile.write("%s\t" % str(info))
                newfile.write("\t".join([str(v) for v in value]))
                newfile.write("\n")

def pos_rearrange(pos_list):
    if pos_list[0] < pos_list[1]: return pos_list
    else:
        tmp=pos_list[0]
        pos_list[0]=pos_list[1]
        pos_list[1]=tmp
        return pos_list

def distance_cal(pos1,pos2):
    #map function, very useful
    pos1=pos_rearrange([int(v) for v in pos1.split("-")])
    pos2=pos_rearrange([int(v) for v in pos2.split("-")])
    # Handy function to calculate distances
    if pos2[0] >pos1[1]: return pos2[0]-pos1[1]
    if pos1[0]<pos2[0]<pos1[1] or pos1[0]<pos2[1]<pos1[1] : return 0
    if pos2[1]<pos1[0]:return pos1[0]-pos2[1]
    else: return 0

def find_complete_CRISPR_Cas_and_SelfTargeting(fna,outputfile,threads,general_folder):
    subprocess.Popen(["cctyper",fna,outputfile,"--no_plot","-t",threads,"--prodigal","meta"]).wait()
    CRISPR_CasTable = subprocess.run(["find", outputfile, "-name", "CRISPR_Cas.tab"], stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
    if CRISPR_CasTable == "":
        print("No complete CRISPR-Cas found input")
        return None
    else:
        print("Complete CRISPR-Cas found and can be found in %s"%CRISPR_CasTable)
        CasTable = subprocess.run(["find", outputfile, "-name", "cas_operons.tab"],
                                  stdout=subprocess.PIPE).stdout.decode(
            'utf-8').rstrip()
        Crispr_Table = subprocess.run(["find", outputfile, "-name", "crisprs_near_cas.tab"],
                                      stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
        df_CasTable = pd.read_csv(CasTable, sep="\t")
        df_Crispr_Table = pd.read_csv(Crispr_Table, sep="\t")
        fna_dic = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

        df_CRISPRcas = pd.read_csv(CRISPR_CasTable, sep="\t")
        CC_list = []
        fna_blastdb = os.path.join(outputfile, os.path.basename(fna) + ".blastDB")
        subprocess.Popen(["makeblastdb", "-dbtype", "nucl", "-in", fna, "-out", fna_blastdb]).wait()
        df_CCtable = pd.DataFrame()
        CC_contig = []
        CC_location = []
        C_operon_withLocation = []
        CC_type = []
        Cas_operon_withLocation = []
        STSS_region = []
        contig_length = []
        for index, CRISPR_Cas in df_CRISPRcas.iterrows():
            location = "-".join(CRISPR_Cas["Operon_Pos"].lstrip("[").rstrip("]").replace(" ", "").split(","))
            CRCasType = CRISPR_Cas["Prediction"]
            spacer_names = CRISPR_Cas["CRISPRs"].lstrip("[").rstrip("]").replace("'", "").replace(" ", "").split(
                ",")  # Some CC have multiple spacers, thisis to get all the spacers
            # print("CRISPR-Cas spacers: ")
            ## Look for self-targeting regions in the genome
            self_targeting_regions = []
            CC_contig.append(CRISPR_Cas["Contig"])
            CC_type.append(CRCasType)
            CC_location.append(location)
            tmp_list = []
            for spacer in spacer_names:
                tmp_list.append(spacer + "|" + "-".join(
                    [df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["Start"].to_string(index=False),
                     df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["End"].to_string(index=False)]))
            C_operon_withLocation.append(";".join(tmp_list))
            Cas_operon_withLocation.append(CRISPR_Cas["Operon"] + "|" + "-".join(
                [df_CasTable[df_CasTable["Operon"] == CRISPR_Cas["Operon"]]["Start"].to_string(index=False),
                 df_CasTable[df_CasTable["Operon"] == CRISPR_Cas["Operon"]]["End"].to_string(index=False)]))
            contig_length.append(len(fna_dic[CRISPR_Cas["Contig"]].seq))
            for spacer in spacer_names:
                spacer_fna = subprocess.run(["find", outputfile, "-name", spacer + ".fa"],
                                            stdout=subprocess.PIPE).stdout.decode('utf-8').rstrip()
                spacer_fna_blastOutfile = os.path.join(outputfile,
                                                       os.path.basename(fna) + "_VS_" + spacer + ".blastnOUT")
                spacer_contig = CRISPR_Cas["Contig"]
                spacer_location = "-".join(
                    [df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["Start"].to_string(index=False),
                     df_Crispr_Table[df_Crispr_Table["CRISPR"] == spacer]["End"].to_string(index=False)])
                CRISPR_spacer_with_cas_locations=str(min([int(q) for q in spacer_location.split("-") + location.split("-")]))+"-"+str(max([int(p) for p in spacer_location.split("-") + location.split("-")]))

                subprocess.Popen(
                    ["blastn", "-query", spacer_fna, "-db", fna_blastdb, "-out", spacer_fna_blastOutfile,
                     "-num_threads",
                     str(threads), "-outfmt", "6"]).wait()
                df_blastOUT = pd.read_csv(spacer_fna_blastOutfile, header=None, sep="\t")
                for index, row in df_blastOUT.iterrows():
                    blast_location = str(row[8]) + "-" + str(row[9])
                    target_contig = str(row[1])
                    if distance_cal(CRISPR_spacer_with_cas_locations, blast_location) > 5000 and target_contig == spacer_contig:
                        # print(spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
                        self_targeting_regions.append(
                            spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
                    elif target_contig != spacer_contig:
                        # print(spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
                        self_targeting_regions.append(spacer + ":" + spacer_location + "=>" + target_contig + ":" + blast_location)
            if len(self_targeting_regions) > 0:
                CC_list.append(CRISPR_Cas["Contig"] + "|" + CRCasType + "|" + location + "|" + "STSS=" + "+".join(
                    self_targeting_regions)) #This location indicate the Cas operon
                STSS_region.append("+".join(self_targeting_regions))
            else:
                CC_list.append(CRISPR_Cas["Contig"] + "|" + CRCasType + "|" + location + "|" + "No_STSS")
                STSS_region.append("No_STSS")
            # Contig|CasTyper|CasPosition|STSS_info
        df_CCtable["Contig"] = CC_contig
        df_CCtable["CRISPR-Cas Type"] = CC_type
        # df_CCtable["CRISPR-Cas Location"] = CC_location
        df_CCtable["CRISPR operon and location"] = C_operon_withLocation
        df_CCtable["Cas operon and location"] = Cas_operon_withLocation
        df_CCtable["STSS"] = STSS_region
        df_CCtable["Contig Length"] = contig_length
        df_CCtable.to_csv(os.path.join(general_folder, "CRISPR-Cas_found.csv"), index=False)
        return CC_list

def find_prophage(fna,outputdir,threads):
    subprocess.Popen(["VIBRANT_run.py","-i",fna,"-folder",outputdir,"-t",threads,"-no_plot"]).wait()
    phage_combined=subprocess.run(["find %s -name *phages_combined.faa"%outputdir],stdout=subprocess.PIPE,shell=True).stdout.decode('utf-8').rstrip()
    if os.path.isfile(phage_combined) and os.path.getsize(phage_combined) > 0:
        print("All prophage found stored in file %s"%phage_combined)
        phage_pos_start_dic = {}
        phage_pos_end_dic = {}
        for line in subprocess.Popen(["grep", ">", phage_combined], stdout=subprocess.PIPE).stdout:
            line = line.decode('utf-8').rstrip().split()
            contig = line[0].lstrip(">") # this is actually proID
            if "fragment" in contig:
                fragment=[contig.split("fragment_")[-1].split("_")[0]]
            else: fragment=[]
            positions_start_end = [int(s) for s in
                                   [v.lstrip("(").rstrip(")") for v in line if ".." in v][0].split("..")]
            if len(fragment) > 0:
                contig_fragment_key = contig.split("_fragment")[0] + "|" + fragment[0]
                if contig_fragment_key not in phage_pos_start_dic:
                    phage_pos_start_dic.setdefault(contig_fragment_key, [positions_start_end[0]])
                else:
                    phage_pos_start_dic[contig_fragment_key].append(positions_start_end[0])
                if contig_fragment_key not in phage_pos_end_dic:
                    phage_pos_end_dic.setdefault(contig_fragment_key, [positions_start_end[1]])
                else:
                    phage_pos_end_dic[contig_fragment_key].append(positions_start_end[1])
            else:
                contig_key = contig.rsplit("_",1)[0]
                if contig_key not in phage_pos_start_dic:
                    phage_pos_start_dic.setdefault(contig_key, [positions_start_end[0]])
                else:
                    phage_pos_start_dic[contig_key].append(positions_start_end[0])
                if contig_key not in phage_pos_end_dic:
                    phage_pos_end_dic.setdefault(contig_key, [positions_start_end[1]])
                else:
                    phage_pos_end_dic[contig_key].append(positions_start_end[1])
        phage_locations=[]
        prophage_out_table = os.path.join(outputdir, "prophage_locations.csv")
        newfile = open(prophage_out_table, "w")
        newfile.write(",".join(["Contig","Start","End","Contig Length"])+"\n")
        fna_dic=SeqIO.to_dict(SeqIO.parse(fna,"fasta"))
        for key in phage_pos_start_dic:
            pos_end = max(phage_pos_start_dic[key])
            pos_start = min(phage_pos_end_dic[key])
            phage_locations.append(key.split("|")[0] + ":" + str(pos_start) + "-" + str(pos_end)) # Contig:startPos-endPos
            newfile.write(",".join([key.split("|")[0],str(pos_start),str(pos_end),str(len(fna_dic[key.split("|")[0]].seq))])+"\n")
        return phage_locations
    else:
        print("No prophage found in sequence")
        return  None

def prophage_harboring_operonFind(df_outTable,prophage_regions):
    operon_location = df_outTable.head(1)["Start"].to_string(index=False) + "-" + df_outTable.tail(1)["End"].to_string(index=False)
    operon_contig = df_outTable.head(1)["Contig ID"].to_string(index=False)
    prophage_withOperon=[v for v in prophage_regions if distance_cal(v.split(":")[-1],operon_location) == 0 and operon_contig == v.split(":")[0]]
    if len(prophage_withOperon)>0:
        return ";".join(prophage_withOperon) # Contig:startPos-endPos
    else:
        return np.nan

def pfamScan_run(operon_faa_file, pfam_hmm_dir,threads,HTH_alignment_evalue):
    pfamOut_file=operon_faa_file+".pfamScanOut"
    subprocess.Popen(["pfam_scan.pl","-fasta",operon_faa_file,"-dir",pfam_hmm_dir,"-outfile",pfamOut_file,"-e_seq",HTH_alignment_evalue,"-cpu",threads]).wait()
    dic_pfam={}
    for record in SeqIO.parse(operon_faa_file,"fasta"):
        info_list=[]
        for line in subprocess.Popen(["grep", record.id,pfamOut_file],stdout=subprocess.PIPE).stdout:
            line=line.decode('utf-8').rstrip().split()
            if len(line)>1:
                info_list.append("|".join([line[5],line[6],line[12],line[14]])) #hmmacc|hmmname|evalue|clan
        if len(info_list)>0:
            dic_pfam.setdefault(record.id,";".join(info_list))
    return dic_pfam

def IR_find(df_allResult,fna_file,sub_outputfolder_path):
    fna_dic=SeqIO.to_dict(SeqIO.parse(fna_file,"fasta"))
    IR_dir=os.path.join(sub_outputfolder_path,"Inverted_Repeats")
    if os.path.isdir(IR_dir) == False:
        os.makedirs(IR_dir)
    operonID=""
    start_pos_list=[]
    operon_withPAL=[]
    for index, row in df_allResult.iterrows():
        if (row["Operon Number"] != operonID and len(start_pos_list)>0) or (index==df_allResult.index[-1]):
            promoter_end=min(start_pos_list)
            if promoter_end <= 400:
                promoter_start=0
            else:
                promoter_start=promoter_end-400
            ### Process IRs ###
            operonContig=df_allResult.loc[[index-1]]["Contig ID"].to_string(index=False)
            # print(promoter_end,promoter_start,operonID,operonContig)
            promoter_seq=str(fna_dic[operonContig].seq[promoter_start:promoter_end])
            # print(promoter_seq)
            operon_promoter_fna=os.path.join(IR_dir,operonID+"_promoter-region.fna")
            with open(operon_promoter_fna,"w") as newfile:
                newfile.write(">"+operonID+"_promoter-region"+"\n")
                newfile.write(promoter_seq)
            operon_promoter_fna_outfile=operon_promoter_fna+".pal"
            subprocess.Popen(["palindrome",operon_promoter_fna,"-minpallen","10","-maxpallen", "100", "-gaplimit", "100", "-nummismatches", "0", "-overlap", "Y", "-outfile",operon_promoter_fna_outfile]).wait()
            pal_found=[]
            for line in subprocess.Popen(['grep -A 100000 "Palindromes:" %s|grep -v "Palindromes:"'%operon_promoter_fna_outfile],shell=True,stdout=subprocess.PIPE).stdout:
                line=line.decode('utf-8').rstrip()
                if line != "": pal_found.append(line)
            if len(pal_found)>0: operon_withPAL.append(operonID)
            ###
            operonID=row["Operon Number"]
            start_pos_list=[row["Start"]]
        else:
            start_pos_list.append(row["Start"])
            operonID = row["Operon Number"]
    operon_withPAL_column=[]
    for index, row in df_allResult.iterrows():
        if row["Operon Number"] in operon_withPAL:
            operon_withPAL_column.append("IR Present")
        else: operon_withPAL_column.append(np.nan)
    df_allResult["IR in Promoter"] = operon_withPAL_column
    return df_allResult


# get list of AO protein and NoAO proteins

def get_protein_labels(AOPF_IDFile, NonAOPF_IDFile):
    # get protein labels for AOPF's
    AO_proteins = []
    with open(AOPF_IDFile) as f:
        for readline in f:
            line_strip = readline.rstrip('\n')
            AO_proteins.append(line_strip)
    # get protein labels for non-AOPF's
    NoAO_proteins = []
    with open(NonAOPF_IDFile) as f:
        for readline in f:
            line_strip = readline.rstrip('\n')
            NoAO_proteins.append(line_strip)

    NoAO_proteins.append("NoAO")
    return AO_proteins, NoAO_proteins


# get state sequence corresponding to protein sequence
def get_state_seq(protein_seq, AO_proteins, NoAO_proteins):
    state_seq = []
    for i in protein_seq:
        if i in AO_proteins:
            state_seq.append("AO")
        else:
            state_seq.append("NoAO")
    return state_seq


# get transition, emission and initial probabilities
# I: initial probabilities for AO and NoAO states
# A: transition probabilities
# B_AO: emission probabilities for AO state
# B_NoAO: emission probabilities for NoAO state

def get_train_HMM_prob(I_file, A_file, B_AO_file, B_NoAO_file):
    with open(I_file) as f:
        data_I = f.read()
    I = json.loads(data_I)

    with open(A_file) as f:
        data_A = f.read()
    A = json.loads(data_A)

    with open(B_AO_file) as f:
        data_B_AO = f.read()
    B_AO = json.loads(data_B_AO)

    with open(B_NoAO_file) as f:
        data_B_NoAO = f.read()
    B_NoAO = json.loads(data_B_NoAO)

    return I, A, B_AO, B_NoAO


# get a list of threshold : [strict_threshold, moderate_threshold, relaxed_threshold]
def get_threshold(T_file):
    with open(T_file, "r") as fp:
        t = json.load(fp)
    return t


# calculate probability for the given protein_seq
def calculate_prob(protein_seq, state_seq, I, A, B_AO, B_NoAO):
    if state_seq[0] == "AO":
        prob = math.log10(I["AO"])
    else:
        prob = math.log10(I["NoAO"])
    for i in range(1, len(protein_seq)):
        if state_seq[i - 1] == "AO" and state_seq[i] == "AO":
            prob = prob + math.log10(A["AA"]) + math.log10(B_AO[str(protein_seq[i])])
        elif state_seq[i - 1] == "AO" and state_seq[i] == "NoAO":
            prob = prob + math.log10(A["AN"]) + math.log10(B_NoAO[str(protein_seq[i])])
        elif state_seq[i - 1] == "NoAO" and state_seq[i] == "AO":
            prob = prob + math.log10(A["NA"]) + math.log10(B_AO[str(protein_seq[i])])
        else:
            prob = prob + math.log10(A["NN"]) + math.log10(B_NoAO[str(protein_seq[i])])
    return prob
# compare with threshold
# strict
# medium
# realxed
def compare_threshold(prob, t, choice):
    if choice == "strict":
        if prob < t[0]:
            output = [prob, True]
        else:
            output = [prob, False]
    elif choice == "medium":
        if prob < t[1]:
            output = [prob, True]
        else:
            output = [prob, False]
    else:
        if prob < t[2]:
            output = [prob, True]
        else:
            output = [prob, False]
    return output

def HMM_output(protein_seq, choice):
    # get protein labels
    AO_proteins, NoAO_proteins = get_protein_labels('./AOPF_IDs.txt', './nonAOPF_IDs.txt')
    # get state sequence corresponding to protein sequence
    state_seq = get_state_seq(protein_seq, AO_proteins, NoAO_proteins)
    # get HMM trained parameters
    I, A, B_AO, B_NoAO = get_train_HMM_prob('./Initial.txt', 'Transition.txt', 'Emission_AO.txt', 'Emission_NoAO.txt')
    # get threshold strict, medium, relaxed
    t = get_threshold('./threshold.txt')
    # get probability of protein sequence
    prob = calculate_prob(protein_seq, state_seq, I, A, B_AO, B_NoAO)
    # return output : [score,True/False], True if AO and False otherwise
    output = compare_threshold(prob, t, choice)
    return output
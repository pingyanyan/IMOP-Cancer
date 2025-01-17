# -*- coding: utf-8 -*-

import subprocess
import os

def process_cancer_data(cancer_type):
    input_file_path = "E:/comp_allcancer/comp_"+cancer_type+"/NDTinput/purity_"+cancer_type+".txt"
    with open(input_file_path, 'r') as input_file:
        for line in input_file:
            columns = line.strip().split()
            if len(columns) >= 2:
                sample_id = columns[0]
                maf_fn = "E:/comp_allcancer/comp_"+cancer_type+"/NDTinput/" + sample_id + "_input.txt"
                purity = columns[1]
                timepoint = 0
                sif_filename = "{}.sif".format(sample_id)
                output_folder = "E:/comp_allcancer/comp_" + cancer_type + "/NDT/" + sample_id
                sif_filepath = os.path.join(output_folder, sif_filename)
                with open(sif_filepath, 'w') as sif_file:
                    sif_file.write("sample_id\tmaf_fn\tseg_fn\tpurity\ttimepoint\n")
                    sif_file.write("{}\t{}\t\t{}\t{}\n".format(sample_id, maf_fn, purity, timepoint))
                m_filename = "E:/comp_allcancer/comp_"+cancer_type+"/NDT/" + sample_id + "/" + sample_id + ".mut_ccfs.txt"
                c_filename = "E:/comp_allcancer/comp_"+cancer_type+"/NDT/" + sample_id + "/" + sample_id + ".cluster_ccfs.txt"
                command = ['python', 'D:/PhylogicNDT/PhylogicNDT.py','BuildTree','-i',sample_id,'-sif',sif_filename,'-m',m_filename,'-c',c_filename]
                return_code = subprocess.call(command,cwd=output_folder)
if __name__ == "__main__":
    cancer_types = ["LUAD","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
          "KICH","KIRC","KIRP","LAML","LGG","LIHC","LUSC","MESO","OV","PAAD",
          "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC",
          "UCS","UVM"]
    for cancer_type in cancer_types:
        process_cancer_data(cancer_type)
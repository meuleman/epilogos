datasets = ["ROADMAP","ADSERA","GORKIN"]

assemblyR = ['hg19','hg38']
assemblyA = ['hg19','hg38']
assemblyG = ['mm10']

statesR = ['15','18','25']
statesA = ['18']
statesG = ['15']

groupsRoadmap15 = ['Adult_Blood_Reference', 'Brain', 'Digestive', 'ESC_derived', 'Immune_and_neurosphere', 'Muscle', 'Primary_Cell',
            'Adult_Blood_Sample', 'Cell_Line', 'ENCODE_2012', 'Female_donors', 'iPSC', 'Neurospheres', 'Primary_Tissue',
            'All_127_Roadmap_epigenomes', 'Cord_Blood_Reference', 'Epithelial', 'Heart', 'Male_donors', 'Non-T-cells', 'Smooth_Muscle',
            'Blood_and_T-cells', 'Cord_Blood_Sample', 'ESC', 'HSC_and_B-cells', 'Mesenchymal', 'Other', 'Thymus']

groupsRoadmap18 = ['Adult_Blood_Reference', 'Brain', 'Digestive', 'ESC_derived', 'Muscle', 'Primary_Cell',
            'Adult_Blood_Sample', 'Cell_Line', 'ENCODE_2012', 'Female_donors', 'iPSC', 'Primary_Tissue',
            'All_127_Roadmap_epigenomes', 'Cord_Blood_Reference', 'Epithelial', 'Heart', 'Male_donors', 'Non-T-cells', 'Smooth_Muscle',
            'Blood_and_T-cells', 'Cord_Blood_Sample', 'ESC', 'HSC_and_B-cells', 'Mesenchymal', 'Other', 'Thymus']

groupsRoadmap25 = ['Adult_Blood_Reference', 'Brain', 'Digestive', 'ESC_derived', 'Muscle', 'Primary_Cell',
            'Adult_Blood_Sample', 'Cell_Line', 'ENCODE_2012', 'Female_donors', 'iPSC', 'Primary_Tissue',
            'All_127_Roadmap_epigenomes', 'Cord_Blood_Reference', 'Epithelial', 'Heart', 'Male_donors', 'Non-T-cells', 'Smooth_Muscle',
            'Blood_and_T-cells', 'Cord_Blood_Sample', 'ESC', 'HSC_and_B-cells', 'Mesenchymal', 'Other', 'Thymus', 'Neurospheres']

groupsAdsera18 = ['All_833_biosamples']

groupsGorkin15 = ['Day-of-birth', 'Embryonic_day_12.5', 'Embryonic_day_15.5', 'Forebrain', 'Intestine', 'Liver', 'Stomach',
        'Digestive_System', 'Embryonic_day_13.5', 'Embryonic_day_16.5', 'Heart', 'Kidney', 'Lung',
        'Embryonic_day_11.5', 'Embryonic_day_14.5', 'Facial_Prominence', 'Hindbrain',  'Limb', 'Neural_Tube','All_65_epigenomes']

combinationFile = open("combinations.txt", "x")

for dataset in datasets:
    if dataset == "ROADMAP":
        for assembly in assemblyR:
            for state in statesR:
                if state == "15":
                    for group in groupsRoadmap15:
                        for saliency in ["S1", "S2", "S3"]:
                            combinationFile.write("{}\t{}\t{}\t{}\t{}\n".format(dataset, assembly, state, group, saliency))
                elif state == "18":
                    for group in groupsRoadmap18:
                        for saliency in ["S1", "S2", "S3"]:
                            combinationFile.write("{}\t{}\t{}\t{}\t{}\n".format(dataset, assembly, state, group, saliency))
                elif state == "25":
                    for group in groupsRoadmap25:
                        for saliency in ["S1", "S2", "S3"]:
                            combinationFile.write("{}\t{}\t{}\t{}\t{}\n".format(dataset, assembly, state, group, saliency))
                else:
                    print("ERROR: ROADMAP State not supported")
    elif dataset == "ADSERA":
        for assembly in assemblyA:
            for state in statesA:
                if state == "18":
                    for group in groupsAdsera18:
                        for saliency in ["S1", "S2"]:
                            combinationFile.write("{}\t{}\t{}\t{}\t{}\n".format(dataset, assembly, state, group, saliency))
                else:
                    print("ERROR: ADSERA State not supported")
    elif dataset == "GORKIN":
        for assembly in assemblyG:
            for state in statesG:
                if state == "15":
                    for group in groupsGorkin15:
                        for saliency in ["S1", "S2", "S3"]:
                            combinationFile.write("{}\t{}\t{}\t{}\t{}\n".format(dataset, assembly, state, group, saliency))
                else:
                    print("ERROR: GORKIN State not supported")

combinationFile.close()
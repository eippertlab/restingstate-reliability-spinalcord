#!/bin/bash

# 2.3.3. Denoising
# 2.3.3.2. Thermal noise

# create template FEAT and PNM files using Feat and PNM GUI and save them under ...derivatives/sub-FeatExample
# copy these for each subject and change the name 

# see the following preprint: " ", and following dataset: https://openneuro.org/datasets/ds004386
# Merve Kaptan, mkaptan@cbs.mpg.de

datapath=/data/pt_02098/RELIABILITY_FC/Reliability_Spinal_RestingStatefMRI/derivatives/

for zMode in auto_denoised manual_denoised
   do
      
    for pnmFiles in max
        do
        
        templateDesignFilePnm=${datapath}/sub-FeatExample/physio/${zMode}_PNM/evlist_${pnmFiles}.txt;
        templateDesignFileFeat=${datapath}/sub-FeatExample/func/FeatFolders/${zMode}_${pnmFiles}.fsf;

        for subj in sub-ZS001 sub-ZS002 sub-ZS003 sub-ZS004 sub-ZS005 sub-ZS006 sub-ZS007 sub-ZS008 sub-ZS010 sub-ZS011 sub-ZS012 sub-ZS013 sub-ZS014 sub-ZS015 sub-ZS016 sub-ZS017 sub-ZS019 sub-ZS020 sub-ZS021 sub-ZS022 sub-ZS023 sub-ZS024 sub-ZS025 sub-ZS026 sub-ZS027 sub-ZS028 sub-ZS029 sub-ZS031 sub-ZS032 sub-ZS033 sub-ZS034 sub-ZS035 sub-ZS036 sub-ZS037 sub-ZS038 sub-ZS039 sub-ZS040 sub-ZS041 sub-ZS042 sub-ZS043 sub-ZS044 sub-ZS045 sub-ZS046 sub-ZS047 sub-ZS048  

            do
            echo $subj

            #sed 's/\/sub-FeatExample/\/'${subj}'/g' $templateDesignFilePnm > ${datapath}/${subj}/physio/${zMode}_PNM/evlist_${pnmFiles}.txt;
            
            
                if [[ ! -f ${datapath}/${subj}/func/FeatFolders ]]
                
                 then
                
                   mkdir ${datapath}/${subj}/func/FeatFolders
                    
                 fi

            sed 's/\/sub-FeatExample/\/'${subj}'/g' $templateDesignFileFeat > ${datapath}/${subj}/func/FeatFolders/${zMode}_${pnmFiles}.fsf;
            
done
done
done

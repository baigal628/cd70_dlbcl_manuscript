#!/bin/bash
#SBATCH --job-name=cp
#SBATCH --mem=30G        # total memory need
#SBATCH -n 8 #Number of cores
#SBATCH -N 1 # on one node
#SBATCH --error=/liulab/galib/hodgkin_lymphoma_mtang/EGA_data/lsf_%j_%x.err
#SBATCH --output=/liulab/galib/hodgkin_lymphoma_mtang/EGA_data/lsf_%j_%x.err
#SBATCH --mail-type=END,FAIL # email notification when job ends/fails
#SBATCH --mail-user=galib@ds.dfci.harvard.edu # email to notify

cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_1/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_1/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_10/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_10/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_11/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_11/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_12/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_12/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_13/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_13/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_14/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_14/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_15/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_15/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_16/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_16/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_17/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_17/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_18/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_18/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_19/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_19/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_2/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_2/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_20/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_20/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_21/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_21/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_22/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_22/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_23/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_23/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_24/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_24/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_25/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_25/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_26/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_26/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_27/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_27/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_28/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_28/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_29/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_29/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_3/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_3/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_30/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_30/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_31/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_31/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_32/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_32/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_33/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_33/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_34/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_34/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_35/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_35/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_36/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_36/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_38/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_38/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_39/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_39/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_4/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_4/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_40/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_40/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_41/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_41/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_42/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_42/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_43/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_43/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_44/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_44/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_45/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_45/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_46/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_46/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_47/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_47/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_48/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_48/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_49/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_49/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_5/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_5/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_50/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_50/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_51/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_51/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_52/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_52/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_53/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_53/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_54/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_54/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_55/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_55/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_56/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_56/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_57/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_57/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_58/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_58/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_59/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_59/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_6/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_6/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_60/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_60/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_61/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_61/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_62/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_62/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_63/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_63/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_64/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_64/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_65/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_65/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_66/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_66/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_67/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_67/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_68/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_68/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_69/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_69/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_7/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_7/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_70/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_70/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_71/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_71/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_72/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_72/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_73/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_73/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_74/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_74/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_75/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_75/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_76/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_76/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_78/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_78/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_79/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_79/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_8/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_8/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_80/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_80/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_81/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_81/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_82/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_82/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_83/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_83/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_84/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_84/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_85/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_85/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_86/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_86/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_87/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_87/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_88/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_88/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_89/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_89/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_9/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_9/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_90/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_90/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool105_91/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool105_91/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool106-1_52/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool106-1_52/
cp -r /liulab/galib/mouse_scRNAseq_margaret/data/GEXdata/Pool106-1_53/outs/filtered_feature_bc_matrix/ /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool106-1_53/


echo "done copying files"

mkdir /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/all_analysis
mv  /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/Pool10*  /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/all_analysis/
tar -czvf /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/all_analysis.tar.gz /liulab/galib/mouse_scRNAseq_margaret/data/EGA_data/GEXdata/all_analysis/

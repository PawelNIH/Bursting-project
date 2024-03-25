Repository for "3D chromatin architecture, BRD4 and Mediator have distinct roles in regulating genome-wide transcriptional bursting and gene network"

Authors
Pawel Trzaskoma1†*, SeolKyoung Jung1†, Aleksandra Pękowska5†, Christopher H. Bohrer2, Xiang Wang1‡, Faiza Naz1, Stefania Dell'Orso1, Wendy D. Dubois2, Ana Olivera3, Supriya V. Vartak1, Yongbing Zhao1, Subhashree Nayak1, Andrew Overmiller1, Maria I. Morasso1, Vittorio Sartorelli1, Daniel R. Larson2, Carson C. Chow4*, Rafael Casellas1, John J. O'Shea1*

Affiliations 
1.	National Institute of Arthritis and Musculoskeletal and Skin Diseases, National Institutes of Health, Bethesda, MD, USA.
2.	National Cancer Institute, National Institutes of Health, Bethesda, MD, USA.
3.	National Institute of Allergy and Infectious Diseases, National Institutes of Health, Bethesda, MD, USA.
4.	National Institute of Diabetes and Digestive and Kidney Diseases, National Institutes of Health, Bethesda, MD, USA.
5.	Dioscuri Center of Chromatin Biology and Epigenomics, Nencki Institute of Experimental Biology, Polish Academy of Sciences/ National Institute of Arthritis and Musculoskeletal and Skin Diseases, National Institutes of Health, Bethesda, MD, USA.

†These authors contributed equally to this work: Pawel Trzaskoma, SeolKyoung Jung, Aleksandra Pękowska.

‡Present address: Children's Hospital of Philadelphia, Philadelphia, PA, USA.

*Correspondence: Pawel Trzaskoma (pawel.trzaskoma@nih.gov), Carson C. Chow (StochasticGene) (carsonc@niddk.nih.gov) and John J. O'Shea (osheaj@arb.niams.nih.gov).



This repository includes two folders:
1)	Summary_files, where we included transcription parameters for all experiments included in the paper: rates along with summary statistics.
2)	Codes, where we included:
•	Code used for preprocessing of raw scRNA-Seq data (Preprocessing.r), which was used for filtering and generating count matrices used for the fitting. The raw sequencing and processed data (count matrices for all samples) are available on GEO (GSE241338),
•	Code used to fit data (StochasticGene_fit.rtf),
•	Code used for QC filtering of inferred rates (QC_filtering.R),
•	Codes used to generate all plots.

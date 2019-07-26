Readme file for codes and data used in the following publication:		

# Lachmuth S, Molofsky J, Suda J, Milbrath L, Keller SR (accepted for AoB Plants on 17 June 2019) 
# Associations between genomic ancestry, genome size and capitula morphology in the invasive meadow knapweed hybrid complex 
# (Centaurea ×moncktonii C.E. Britton) in eastern North America


For further inquiry, contact Susanne Lachmuth (susanne.lachmuth@botanik.uni-halle.de)				
and/or Stephen R. Keller (corresponding author, srkeller@uvm.edu)	


R-codes (.R):
See file name and comments in the files 

Data files (in folder: data)
File: capitula_data.txt	contains morphometric capitula data																		
File: capitula_data_PCA.txt	contains selected capitula data and taxon based on admixture and NewHybrids analyses																		
File: ind_HybClasses_Cjacea_Panel.txt	contains individual assignment probabilities to the parental taxa and various hybrid classes according to the NewHybrids analysis based on the C. jacea marker panel																		
File: ind_HybClasses_Cnigra_Panel.txt	contains individual assignment probabilities to the parental taxa and various hybrid classes according to the NewHybrids analysis based on the C. nigra marker panel																		
File: ind_metadata.txt	contains individual metadata (pop, taxon)																		
File: ind_SNP_data.vcf	contains individual SNP data																		
File: PCA_data_CjaceaPanel.txt	contains individual scores from morphometric and genomic principal component analyses, genome size and hybrid class according to the NewHybrids analysis based on the C. nigra marker panel																		
File: PCA_data_CnigraPanel.txt	contains individual scores from morphometric and genomic principal component analyses, genome size and hybrid class according to the NewHybrids analysis based on the C. jacea marker panel																		
File: pop_genetic_data.txt	contains individual admixture assignment probabilities, genome sizes, genomic PCA scores and NewHybrids hybrid classes based on both marker panels																		
 
 
Description of variables across all data files:
Variable		Description																		
AnK2_K1			individual assignment probabilities to cluster 1 in the admixture analysis with number of clusters K=2																		
AnK2_K2			individual assignment probabilities to cluster 2 in the admixture analysis with number of clusters K=2																		
AnK3_K3			individual assignment probabilities to cluster 2 in the admixture analysis with number of clusters K=3																		
ApCen_Ratio		ratio appendage center width to length																		
Brac_appress		bract appression																		
Brac_Col		bract color																		
Brac_Len		bract length																		
Brac_Ratio		ratio bract width to length																		
Brac_RowNr		number of bract rows																		
Brac_Wid		bract width																		
BracApCen_Len		bract appendage center length																		
BracApCen_Wid		bract appendage center width																		
Cap			capitulum identifier																		
Cap_Len			capitulum length																		
Cap_Ratio		ratio capitulum width to length																		
Cap_Wid			capitulum width																		
GenSize			genome size																		
HybClass_LDjacea	taxon to which individual was assigned with highest probability in NewHybrids analysis based on C. jacea marker panel																		
HybClass_LDnigra	taxon to which individual was assigned with highest probability in NewHybrids analysis based on C. nigra marker panel																		
ind			individual identifier including population identifier																		
IndID			individual identifier (nr.) without population identifier																		
PC1			individual score at genomic PCA axis 1																		
PC1_morpho		individual score at morphological PCA axis 1																		
PC2			individual score at genomic PCA axis 2																		
PC2_morpho		individual score at morphological PCA axis 2																		
Perc_SerBrac		percentage pectinate bract rows																		
pop /  Pop		sampled population / sampling location																		
pop_ind_cap		identifier: population, indivdual, capitulum																		
relLen_BracApCen	relative length of bract appendage center																		
relWid_BracApCen	relative width of bract appendage center																		
SerBrac_RowNr		number of pectinate bract rows																		
SerBrac_yes_no		Pectinate bracts yes / no																		
taxon / Taxon		Centaurea cf. jacea, C. cf. nigra, hybrid from New York, hybrid from Vermont																		
Taxon admixture		taxon to which individual was assigned with highest probability in admixture analysis																		


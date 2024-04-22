# A high-resolution model of gene expression during _Gossypium hirsutum_ (cotton) fiber development

The code and files contained within this repo pertain to an analysis of the _Gossypium hirsutum_ cv TM1 transcriptome during fiber development (6 to 24 days post anthesis; DPA). In this project, RNA-seq was generated in triplicate for developing cotton fibers, sampling daily at midday. Differential gene expression (DGE) and network analyses were conducted, focusing on genes pertaining to cellulose synthesis and turgor pressure, both of which are important in developing cotton fibers. 

__Organization of files__:  
- __1_makeReference__  
  - 0_makeReference.sh: commands for making the pseudoAD1 reference from the _G. raimondii_ (Paterson et al. 2012) reference genome, including genome download location  
  - pseudogenome_by_snp.py: python code to replace nucleotides in the _G. raimondii_ genome with species specific SNP  
- __2_performMapping__  
  - mapping: folder of Kallisto pseudomapped counts. Each subfolder corresponds to a sample.  
  - 2_indexKallisto: slurm submission file to index the transcriptome  
  - 3_runKallisto: slurm array to map via Kallisto  
  - 4_prepareDEfiles: bash commands to merge the individual count files  
  - 5_getMappingRate: bash commands to get mapping rate for each file  
- __3_conductAnalyses__: R files used to analyze the data. Each is numbered in order that the analyses were run (e.g., 00.ReadWriteData.R was the first one run). Output files typically begin with 'out-' and are numbered to correspond with the script that generated them. Other files are used within the scripts.  
- __4_FiguresTables__: high-resolution of figures and tables, including supplementary.    

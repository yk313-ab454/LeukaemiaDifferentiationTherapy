# LeukaemiaDifferentiationTherapy
Please follow the following steps before running the code:
- Download gsea software from the following link:
http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/resources/software/gsea2-2.2.4.jar

- Place the downloaded jar file in the gsea folder.

- java 8 installation is a requirement.

- java should be added to the environmental variables in windows or added to the Path in linux.

- Run leukaemiaDifferentiationTherapy.r

- After running the R file, the rank ordered list of compounds would be generated in the results folder (GSE48558.txt). Some detailed gsea output would be generated in the output directory which can be ignored.

- Refer to following link, under “Preparing Data Files for GSEA”, for description of file types in the data folder that are input to gsea software (gct, cls, chip).
https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Other_Ways_to


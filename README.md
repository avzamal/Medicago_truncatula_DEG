
# Medicago_truncatula_DEG

This project is about revealing differentially expressed genes in resistant and suseptible plants in different phases of infection.

DEG identification was provided using R script using edgeR. Packages dplyr, tidyr, BiocManager, edgeR, readr are required.
Pantherdb.org was used for the GO annotation. agriGO (http://bioinfo.cau.edu.cn/agriGO/analysis.php) was used for SEA.

Extracting information from outputs was conducted using python script. Python 3.7 is required.

## Identification of differentially expressed genes

We identified genes using R script, based on edgeR package functions. At the end we received a list of genes that can be saved as a file. Example output:

"Medtr4g012000.2"
"Medtr1g033390.1"
"Medtr3g048910.1"
"Medtr8g102550.1"
"Medtr5g075955.3"

Using online tool (http://bioinformatics.psb.ugent.be/webtools/Venn/) we built Venn diagram of received gene groups. In late infection phase we identified 928 DEG regulated by genotype (resistant or not), 563 regulated by inoculation and 67 regulated by resistance:inoculation (genes that are regulated by inoculation in resistance line only). For the further analysis we selected genes that regulated by only genotype, only inoculation and resistance:inoculation.

![Plant_project_pict 001](https://user-images.githubusercontent.com/70773075/122025845-610cff00-cdd2-11eb-999b-900b4d2a1755.jpeg)

## GO annotation

Received list of genes can be converted to Pantherdb acessible format using first script from jupyter notebook. Obtained list can be downloaded to Pantherdb. Medicago truncatula should be selected as organism. In the result only columns "Gene ID", "Gene Name, Gene Symbol, Ortholog", "GO database MF Complete", "GO database BP Complete", "GO database CC Complete" should be selected. This table should be sent as file. The second script extracts GeneID and GO terms from received file. Obtained list can be used for SEA analysis in AgriGO.

We provided SEA (using AgriGO) for all DEG groups. However, in group of downregulated genes by inoculation was only 3 genes that is not enough for the SEA. For a group of downregulated genes by genotype:inoculation biological process terms were not found.
In the table only the significantly enriched GO terms of deeper level are presented.

![Снимок экрана 2021-05-21 в 10 30 38](https://user-images.githubusercontent.com/70773075/122026336-c5c85980-cdd2-11eb-8c95-3767a95ebb6e.png)








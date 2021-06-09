
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

## GO annotation

Received list of genes can be converted to Pantherdb acessible format using first script from jupyter notebook. Obtained list can be downloaded to Pantherdb. Medicago truncatula should be selected as organism. In the result only columns "Gene ID", "Gene Name, Gene Symbol, Ortholog", "GO database MF Complete", "GO database BP Complete", "GO database CC Complete" should be selected. This table should be sent as file. The second script extracts GeneID and GO terms from received file. Obtained list can be used for SEA analysis in AgriGO.






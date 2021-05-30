library(dplyr)
library(tidyr)
library(BiocManager)
library(edgeR)
library(readr)


MACEVerti_0EIL_MtrV4 <- read_delim("MACEVerti_0EIL_MtrV4.tsv", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)

df <- as_tibble(MACEVerti_0EIL_MtrV4)[,1:30] #selected columns with information

# subsets for different time points
df_timeE <- df[,c(1,3,6,10,13,17,20,24,27)]
df_timeI <- df[,c(1,4,7,11,14,18,21,25,28)]
df_timeL <- df[,c(1,5,8,12,15,19,22,26,29)]

#Creation of DGELists
groupE <- c(2,5,9,12,2,5,9,12)
groupI <- c(3,6,10,13,3,6,10,13)
groupL <- c(4,7,11,14,4,7,11,14)
DGEList_E <- DGEList(counts=df_timeE[,2:9], group = groupE)
DGEList_I <- DGEList(counts=df_timeI[,2:9], group = groupI)
DGEList_L <- DGEList(counts=df_timeL[,2:9], group = groupL)

DGEList_E$genes <- df_timeE[,1]
DGEList_I$genes <- df_timeI[,1]
DGEList_L$genes <- df_timeL[,1]

#Filtering the DGEList

resistance <- factor(c(1,1,0,0,1,1,0,0)) # 1 = R, 0 = S
inoculation <- factor(c(0,1,0,1,0,1,0,1))
replicate <- factor(c(1,1,1,1,2,2,2,2))

# model matrix
designE <- model.matrix(~ replicate + inoculation * resistance, data = DGEList_E)
designI <- model.matrix(~ replicate + inoculation * resistance, data = DGEList_I)
designL <- model.matrix(~ replicate + inoculation * resistance, data = DGEList_L)


keepE <- filterByExpr(DGEList_E, designE, min.count = 5)
DGEList_filtered_E <- DGEList_E[keepE, , keep.lib.sizes=FALSE]
keepI <- filterByExpr(DGEList_I, designI, min.count = 5)
DGEList_filtered_I <- DGEList_I[keepI, , keep.lib.sizes=FALSE]
keepL <- filterByExpr(DGEList_L, designL, min.count = 5)
DGEList_filtered_L <- DGEList_L[keepL, , keep.lib.sizes=FALSE]

#Normalization
DGEList_norm_E <- calcNormFactors(DGEList_filtered_E)
DGEList_norm_I <- calcNormFactors(DGEList_filtered_I)
DGEList_norm_L <- calcNormFactors(DGEList_filtered_L)
# dispersion estimation (multifactor)
DGEList_norm_E <- estimateDisp(DGEList_norm_E, designE)
DGEList_norm_I <- estimateDisp(DGEList_norm_I, designI)
DGEList_norm_L <- estimateDisp(DGEList_norm_L, designL)

# model fitting
fit_E <- glmQLFit(DGEList_norm_E, designE)
fit_I <- glmQLFit(DGEList_norm_I, designI)
fit_L <- glmQLFit(DGEList_norm_L, designL)

# for late dataset I have got DE genes for resistance, inoculation, resistance:inoculation
qlf_L_ino_tr <- glmTreat(fit_L, coef=3, lfc=1)
summary(decideTests(qlf_L_ino_tr))
L_ino <- topTags(qlf_L_ino_tr, n = 563)
qlf_L_resist_tr <- glmTreat(fit_L, coef=4, lfc=1)
summary(decideTests(qlf_L_resist_tr))
L_resist <- topTags(qlf_L_resist_tr, n = 928)
qlf_L_ino_resist_tr <- glmTreat(fit_L, coef=5, lfc=1)
summary(decideTests(qlf_L_ino_resist_tr))
L_ino_resist <- topTags(qlf_L_ino_resist_tr, n = 67)

# I have filtered inoculation only DE genes in L_ino_only and resistance only into L_resist_only
L_ino_only <- anti_join(L_ino$table, L_resist$table, by = 'Geneid')
L_ino_only <- anti_join(L_ino_only, L_ino_resist$table, by = 'Geneid')

L_resist_only <- anti_join(L_resist$table, L_ino$table, by = 'Geneid')
L_resist_only <- anti_join(L_resist_only, L_ino_resist$table, by = 'Geneid')
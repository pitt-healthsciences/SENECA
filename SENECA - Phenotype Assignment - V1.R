###########################################################
# Jason Kennedy                                           #
# jnk28@pitt.edu                                          #
# Uploaded June 27, 2022                                  #
#                                                         #
# Assign Clinical Sepsis Phenotypes by Euclidean Distance #
# (Alpha, Beta, Gamma, Delta)                             #
# Using SENECA Derivation Centroids                       #
#                                                         #
# From "Derivation, Validation, and Potential Treatment   #
# Implications of Novel Clinical Phenotypes for Sepsis"   #
# Seymour CW, Kennedy JN, Wang S, et al.                  #
# JAMA. 2019 May 28;321(20):2003-2017                     #
# DOI: 10.1001/jama.2019.5791                             #
# PMID: 31104070; PMCID: PMC6537818                       #
###########################################################

###NOTES###
#1.) Code assumes units are as follows:
# Age: years
# Albumin: g/dL
# ALT: U/L
# AST: U/L
# BANDS: %
# Bicarbonate: mEq/L
# Bilirubin: mg/dL
# BUN: mg/dL
# Chloride: mEq/L
# Creatinine: mg/dL
# C-Reactive Protein: mg/L
# Elixhauser: 0-31 point scale
# ESR: mm/hour
# GCS: 3-15 point scale
# Glucose: mg/dL
# Hemoglobin: g/dL
# Heart Rate: Beats/min
# INR: ratio
# Lactate (Serum)e: mmol/L
# PaO2: mmHg
# Platelets: x10^9/L
# Respiratory Rate: breaths/minute
# SaO2 - Oxygen Saturation: %
# Sex: 1=male; 0=female
# Sodium: mEq/L
# Systolic Blood Pressure: mmHg
# Temperature: Celcius
# Troponin: ng/mL
# White Blood Cell count: x10^9/L
# Convert units to above prior to phenotype assignment as needed
#
#2.) Expected input has subject identifier in 1st column and model features in remaining columns

#Load Libraries
library(missRanger)

#Load Data
setwd("C:\\SENECA Phenotype\\Analysis") #Change to appropriate directory
data <- read.csv("SENECA - Mock Dataset.csv")

#Assign column names -- NOTE: Columns may need to be reordered or list modified to match your dataset
colnames(data) <- c('id'    ,'age'   ,'alb' ,'alt' ,'ast',
                    'bands' ,'bicarb','bili','bun' ,'cl',
                    'creat' ,'crp'   ,'elix','esr' ,'gcs',
                    'gluc'  ,'hgb'   ,'hr'  ,'inr' ,'lactate',
                    'pao2'  ,'plt'   ,'rr'  ,'sao2','sex',
                    'sodium','sbp'   ,'temp','trop','wbc')

###IMPUTATION BY RANDOM FOREST WITH PMM###
set.seed = 12081023
n.features = 29 # Change to number of features in your data
data.imputed    <- missRanger(data[,c(2:(n.features+1))],pmm.k=5,num.trees=1000,seed=set.seed)
data.imputed$id <- data$id
data.imputed    <- data.imputed[,c((n.features+1),(1:n.features))]

###ASSIGN PHENOTYPES###
#Log Transform / Inverse Transform as appropriate
data.logtrans         <- data.imputed
data.logtrans$alt     <- log(data.imputed$alt)
data.logtrans$ast     <- log(data.imputed$ast)
data.logtrans$bands   <- log(data.imputed$bands)
data.logtrans$bili    <- log(data.imputed$bili)
data.logtrans$bun     <- log(data.imputed$bun)
data.logtrans$creat   <- log(data.imputed$creat)
data.logtrans$crp     <- log(data.imputed$crp)
data.logtrans$esr     <- log(data.imputed$esr)
data.logtrans$gluc    <- log(data.imputed$gluc)
data.logtrans$inr     <- log(data.imputed$inr)
data.logtrans$lactate <- log(data.imputed$lactate)
data.logtrans$plt     <- log(data.imputed$plt)
data.logtrans$sbp     <- log(data.imputed$sbp)
data.logtrans$trop    <- log(data.imputed$trop)
data.logtrans$wbc     <- log(data.imputed$wbc)
data.logtrans$sao2    <- log(101-data.imputed$sao2)

#Z-Transform Data Using Original SENECA Derivation Means and Standard Deviations
#Note: means of ln-transformed and inv-ln transformed variables used where appropriate (e.g. mean and SD of ln(ALT))
data.ztrans <- data.logtrans
data.ztrans$age     <- (data.logtrans$age     - 64.41131) / 17.11103
data.ztrans$alb     <- (data.logtrans$alb     - 2.933706) / 0.723259
data.ztrans$alt     <- (data.logtrans$alt     - 3.538419) / 0.901215
data.ztrans$ast     <- (data.logtrans$ast     - 3.598596) / 0.969935
data.ztrans$bands   <- (data.logtrans$bands   - 1.821413) / 1.118881
data.ztrans$bicarb  <- (data.logtrans$bicarb  - 25.03647)	/ 5.131436
data.ztrans$bili    <- (data.logtrans$bili    + 0.143662) / 0.840652
data.ztrans$bun     <- (data.logtrans$bun     - 3.171860) / 0.712280
data.ztrans$cl      <- (data.logtrans$cl      - 102.7818) / 6.715770
data.ztrans$creat   <- (data.logtrans$creat   - 0.425669) / 0.667512
data.ztrans$crp     <- (data.logtrans$crp     - 1.504626) / 1.858398
data.ztrans$elix    <- (data.logtrans$elix    - 1.817871) / 1.170719
data.ztrans$esr     <- (data.logtrans$esr     - 3.738153) / 0.910863
data.ztrans$gcs     <- (data.logtrans$gcs     - 12.83853) / 3.127180
data.ztrans$gluc    <- (data.logtrans$gluc    - 4.964519) / 0.448293
data.ztrans$hgb     <- (data.logtrans$hgb     - 11.50954) / 2.329503
data.ztrans$hr      <- (data.logtrans$hr      - 97.17861) / 21.92295
data.ztrans$inr     <- (data.logtrans$inr     - 0.367357) / 0.403954
data.ztrans$lactate <- (data.logtrans$lactate - 0.496778) / 0.676092
data.ztrans$pao2    <- (data.logtrans$pao2    - 109.4560) / 76.86320
data.ztrans$plt     <- (data.logtrans$plt     - 5.139547) / 0.654274
data.ztrans$rr      <- (data.logtrans$rr      - 22.16539) / 6.146117
data.ztrans$sao2    <- (data.logtrans$sao2    - 1.818297) / 0.767377
data.ztrans$sex     <- (data.logtrans$sex     - 0.496408) / 0.499999
data.ztrans$sodium  <- (data.logtrans$sodium  - 137.1170) / 5.522682
data.ztrans$sbp     <- (data.logtrans$sbp     - 4.678142) / 0.270757
data.ztrans$temp    <- (data.logtrans$temp    - 36.98350) / 1.006374
data.ztrans$trop    <- (data.logtrans$trop    + 2.288375) / 1.230468
data.ztrans$wbc     <- (data.logtrans$wbc     - 2.237017) / 0.716883

#Generate distances for phenotype mapping
#Alpha
dist.alpha              <- data.frame(data[,1])
dist.alpha$d1sq_age     <- ifelse(is.na(data.ztrans$age),    NA,(data.ztrans$age     + 0.282231680)^2)
dist.alpha$d1sq_alb     <- ifelse(is.na(data.ztrans$alb),    NA,(data.ztrans$alb     - 0.716941788)^2)
dist.alpha$d1sq_alt     <- ifelse(is.na(data.ztrans$alt),    NA,(data.ztrans$alt     - 0.003226264)^2)
dist.alpha$d1sq_ast     <- ifelse(is.na(data.ztrans$ast),    NA,(data.ztrans$ast     + 0.163922218)^2)
dist.alpha$d1sq_bands   <- ifelse(is.na(data.ztrans$bands),  NA,(data.ztrans$bands   + 0.253857263)^2)
dist.alpha$d1sq_bicarb  <- ifelse(is.na(data.ztrans$bicarb), NA,(data.ztrans$bicarb  - 0.311171647)^2)
dist.alpha$d1sq_bili    <- ifelse(is.na(data.ztrans$bili),   NA,(data.ztrans$bili    + 0.002667732)^2)
dist.alpha$d1sq_bun     <- ifelse(is.na(data.ztrans$bun),    NA,(data.ztrans$bun     + 0.659563022)^2)
dist.alpha$d1sq_cl      <- ifelse(is.na(data.ztrans$cl),     NA,(data.ztrans$cl      + 0.028323656)^2)
dist.alpha$d1sq_creat   <- ifelse(is.na(data.ztrans$creat),  NA,(data.ztrans$creat   + 0.556978928)^2)
dist.alpha$d1sq_crp     <- ifelse(is.na(data.ztrans$crp),    NA,(data.ztrans$crp     + 0.603062721)^2)
dist.alpha$d1sq_elix    <- ifelse(is.na(data.ztrans$elix),   NA,(data.ztrans$elix    + 0.255208039)^2)
dist.alpha$d1sq_esr     <- ifelse(is.na(data.ztrans$esr),    NA,(data.ztrans$esr     + 0.580672802)^2)
dist.alpha$d1sq_gcs     <- ifelse(is.na(data.ztrans$gcs),    NA,(data.ztrans$gcs     + 0.022069890)^2)
dist.alpha$d1sq_gluc    <- ifelse(is.na(data.ztrans$gluc),   NA,(data.ztrans$gluc    + 0.201601890)^2)
dist.alpha$d1sq_hgb     <- ifelse(is.na(data.ztrans$hgb),    NA,(data.ztrans$hgb     - 0.631583437)^2)
dist.alpha$d1sq_hr      <- ifelse(is.na(data.ztrans$hr),     NA,(data.ztrans$hr      + 0.147992075)^2)
dist.alpha$d1sq_inr     <- ifelse(is.na(data.ztrans$inr),    NA,(data.ztrans$inr     + 0.338345030)^2)
dist.alpha$d1sq_lactate <- ifelse(is.na(data.ztrans$lactate),NA,(data.ztrans$lactate + 0.254416263)^2)
dist.alpha$d1sq_pao2    <- ifelse(is.na(data.ztrans$pao2),   NA,(data.ztrans$pao2    + 0.121989635)^2)
dist.alpha$d1sq_plt     <- ifelse(is.na(data.ztrans$plt),    NA,(data.ztrans$plt     + 0.057589639)^2)
dist.alpha$d1sq_rr      <- ifelse(is.na(data.ztrans$rr),     NA,(data.ztrans$rr      + 0.316314201)^2)
dist.alpha$d1sq_sao2    <- ifelse(is.na(data.ztrans$sao2),   NA,(data.ztrans$sao2    - 0.005663652)^2)
dist.alpha$d1sq_sex     <- ifelse(is.na(data.ztrans$sex),    NA,(data.ztrans$sex     - 0.025144433)^2)
dist.alpha$d1sq_sodium  <- ifelse(is.na(data.ztrans$sodium), NA,(data.ztrans$sodium  - 0.065611647)^2)
dist.alpha$d1sq_sbp     <- ifelse(is.na(data.ztrans$sbp),    NA,(data.ztrans$sbp     - 0.303289160)^2)
dist.alpha$d1sq_temp    <- ifelse(is.na(data.ztrans$temp),   NA,(data.ztrans$temp    - 0.126497810)^2)
dist.alpha$d1sq_trop    <- ifelse(is.na(data.ztrans$trop),   NA,(data.ztrans$trop    + 0.226681667)^2)
dist.alpha$d1sq_wbc     <- ifelse(is.na(data.ztrans$wbc),    NA,(data.ztrans$wbc     + 0.211849439)^2)

#Beta
dist.beta               <- data.frame(data[,1])
dist.beta$d2sq_age      <- ifelse(is.na(data.ztrans$age),    NA,(data.ztrans$age     - 0.366160979)^2)
dist.beta$d2sq_alb      <- ifelse(is.na(data.ztrans$alb),    NA,(data.ztrans$alb     - 0.058423097)^2)
dist.beta$d2sq_alt      <- ifelse(is.na(data.ztrans$alt),    NA,(data.ztrans$alt     + 0.336533408)^2)
dist.beta$d2sq_ast      <- ifelse(is.na(data.ztrans$ast),    NA,(data.ztrans$ast     + 0.386172035)^2)
dist.beta$d2sq_bands    <- ifelse(is.na(data.ztrans$bands),  NA,(data.ztrans$bands   + 0.267407242)^2)
dist.beta$d2sq_bicarb   <- ifelse(is.na(data.ztrans$bicarb), NA,(data.ztrans$bicarb  - 0.017263985)^2)
dist.beta$d2sq_bili     <- ifelse(is.na(data.ztrans$bili),   NA,(data.ztrans$bili    + 0.380987513)^2)
dist.beta$d2sq_bun      <- ifelse(is.na(data.ztrans$bun),    NA,(data.ztrans$bun     - 0.700632855)^2)
dist.beta$d2sq_cl       <- ifelse(is.na(data.ztrans$cl),     NA,(data.ztrans$cl      - 0.009290949)^2)
dist.beta$d2sq_creat    <- ifelse(is.na(data.ztrans$creat),  NA,(data.ztrans$creat   - 0.788307539)^2)
dist.beta$d2sq_crp      <- ifelse(is.na(data.ztrans$crp),    NA,(data.ztrans$crp     + 0.137781835)^2)
dist.beta$d2sq_elix     <- ifelse(is.na(data.ztrans$elix),   NA,(data.ztrans$elix    - 0.467208085)^2)
dist.beta$d2sq_esr      <- ifelse(is.na(data.ztrans$esr),    NA,(data.ztrans$esr     - 0.281220424)^2)
dist.beta$d2sq_gcs      <- ifelse(is.na(data.ztrans$gcs),    NA,(data.ztrans$gcs     - 0.234033835)^2)
dist.beta$d2sq_gluc     <- ifelse(is.na(data.ztrans$gluc),   NA,(data.ztrans$gluc    - 0.020278755)^2)
dist.beta$d2sq_hgb      <- ifelse(is.na(data.ztrans$hgb),    NA,(data.ztrans$hgb     + 0.286660522)^2)
dist.beta$d2sq_hr       <- ifelse(is.na(data.ztrans$hr),     NA,(data.ztrans$hr      + 0.590184626)^2)
dist.beta$d2sq_inr      <- ifelse(is.na(data.ztrans$inr),    NA,(data.ztrans$inr     + 0.055712128)^2)
dist.beta$d2sq_lactate  <- ifelse(is.na(data.ztrans$lactate),NA,(data.ztrans$lactate + 0.404554650)^2)
dist.beta$d2sq_pao2     <- ifelse(is.na(data.ztrans$pao2),   NA,(data.ztrans$pao2    - 0.019043137)^2)
dist.beta$d2sq_plt      <- ifelse(is.na(data.ztrans$plt),    NA,(data.ztrans$plt     - 0.147947739)^2)
dist.beta$d2sq_rr       <- ifelse(is.na(data.ztrans$rr),     NA,(data.ztrans$rr      + 0.337499772)^2)
dist.beta$d2sq_sao2     <- ifelse(is.na(data.ztrans$sao2),   NA,(data.ztrans$sao2    + 0.272829662)^2)
dist.beta$d2sq_sex      <- ifelse(is.na(data.ztrans$sex),    NA,(data.ztrans$sex     + 0.040713400)^2)
dist.beta$d2sq_sodium   <- ifelse(is.na(data.ztrans$sodium), NA,(data.ztrans$sodium  - 0.073153049)^2)
dist.beta$d2sq_sbp      <- ifelse(is.na(data.ztrans$sbp),    NA,(data.ztrans$sbp     - 0.333667198)^2)
dist.beta$d2sq_temp     <- ifelse(is.na(data.ztrans$temp),   NA,(data.ztrans$temp    + 0.319145547)^2)
dist.beta$d2sq_trop     <- ifelse(is.na(data.ztrans$trop),   NA,(data.ztrans$trop    + 0.187150084)^2)
dist.beta$d2sq_wbc      <- ifelse(is.na(data.ztrans$wbc),    NA,(data.ztrans$wbc     + 0.056699022)^2)

#Gamma
dist.gamma              <- data.frame(data[,1])
dist.gamma$d3sq_age     <- ifelse(is.na(data.ztrans$age),    NA,(data.ztrans$age     - 0.015976043)^2)
dist.gamma$d3sq_alb     <- ifelse(is.na(data.ztrans$alb),    NA,(data.ztrans$alb     + 0.694629256)^2)
dist.gamma$d3sq_alt     <- ifelse(is.na(data.ztrans$alt),    NA,(data.ztrans$alt     + 0.226549811)^2)
dist.gamma$d3sq_ast     <- ifelse(is.na(data.ztrans$ast),    NA,(data.ztrans$ast     + 0.139338316)^2)
dist.gamma$d3sq_bands   <- ifelse(is.na(data.ztrans$bands),  NA,(data.ztrans$bands   - 0.298891480)^2)
dist.gamma$d3sq_bicarb  <- ifelse(is.na(data.ztrans$bicarb), NA,(data.ztrans$bicarb  - 0.037482070)^2)
dist.gamma$d3sq_bili    <- ifelse(is.na(data.ztrans$bili),   NA,(data.ztrans$bili    - 0.000477357)^2)
dist.gamma$d3sq_bun     <- ifelse(is.na(data.ztrans$bun),    NA,(data.ztrans$bun     + 0.104459930)^2)
dist.gamma$d3sq_cl      <- ifelse(is.na(data.ztrans$cl),     NA,(data.ztrans$cl      + 0.203177313)^2)
dist.gamma$d3sq_creat   <- ifelse(is.na(data.ztrans$creat),  NA,(data.ztrans$creat   + 0.260659215)^2)
dist.gamma$d3sq_crp     <- ifelse(is.na(data.ztrans$crp),    NA,(data.ztrans$crp     - 0.702551188)^2)
dist.gamma$d3sq_elix    <- ifelse(is.na(data.ztrans$elix),   NA,(data.ztrans$elix    + 0.102982737)^2)
dist.gamma$d3sq_esr     <- ifelse(is.na(data.ztrans$esr),    NA,(data.ztrans$esr     - 0.685360052)^2)
dist.gamma$d3sq_gcs     <- ifelse(is.na(data.ztrans$gcs),    NA,(data.ztrans$gcs     - 0.164700720)^2)
dist.gamma$d3sq_gluc    <- ifelse(is.na(data.ztrans$gluc),   NA,(data.ztrans$gluc    - 0.053908265)^2)
dist.gamma$d3sq_hgb     <- ifelse(is.na(data.ztrans$hgb),    NA,(data.ztrans$hgb     + 0.456098749)^2)
dist.gamma$d3sq_hr      <- ifelse(is.na(data.ztrans$hr),     NA,(data.ztrans$hr      - 0.543281654)^2)
dist.gamma$d3sq_inr     <- ifelse(is.na(data.ztrans$inr),    NA,(data.ztrans$inr     - 0.084362244)^2)
dist.gamma$d3sq_lactate <- ifelse(is.na(data.ztrans$lactate),NA,(data.ztrans$lactate - 0.185959917)^2)
dist.gamma$d3sq_pao2    <- ifelse(is.na(data.ztrans$pao2),   NA,(data.ztrans$pao2    + 0.146087373)^2)
dist.gamma$d3sq_plt     <- ifelse(is.na(data.ztrans$plt),    NA,(data.ztrans$plt     - 0.037279715)^2)
dist.gamma$d3sq_rr      <- ifelse(is.na(data.ztrans$rr),     NA,(data.ztrans$rr      - 0.511269045)^2)
dist.gamma$d3sq_sao2    <- ifelse(is.na(data.ztrans$sao2),   NA,(data.ztrans$sao2    - 0.266455863)^2)
dist.gamma$d3sq_sex     <- ifelse(is.na(data.ztrans$sex),    NA,(data.ztrans$sex     + 0.042400073)^2)
dist.gamma$d3sq_sodium  <- ifelse(is.na(data.ztrans$sodium), NA,(data.ztrans$sodium  + 0.270658297)^2)
dist.gamma$d3sq_sbp     <- ifelse(is.na(data.ztrans$sbp),    NA,(data.ztrans$sbp     + 0.399313346)^2)
dist.gamma$d3sq_temp    <- ifelse(is.na(data.ztrans$temp),   NA,(data.ztrans$temp    - 0.331493219)^2)
dist.gamma$d3sq_trop    <- ifelse(is.na(data.ztrans$trop),   NA,(data.ztrans$trop    + 0.080529879)^2)
dist.gamma$d3sq_wbc     <- ifelse(is.na(data.ztrans$wbc),    NA,(data.ztrans$wbc     - 0.158379065)^2)

#Delta
dist.delta              <- data.frame(data[,1])
dist.delta$d4sq_age     <- ifelse(is.na(data.ztrans$age),    NA,(data.ztrans$age     + 0.087936040)^2)
dist.delta$d4sq_alb     <- ifelse(is.na(data.ztrans$alb),    NA,(data.ztrans$alb     + 0.499133452)^2)
dist.delta$d4sq_alt     <- ifelse(is.na(data.ztrans$alt),    NA,(data.ztrans$alt     - 1.144945208)^2)
dist.delta$d4sq_ast     <- ifelse(is.na(data.ztrans$ast),    NA,(data.ztrans$ast     - 1.486652343)^2)
dist.delta$d4sq_bands   <- ifelse(is.na(data.ztrans$bands),  NA,(data.ztrans$bands   - 0.579760956)^2)
dist.delta$d4sq_bicarb  <- ifelse(is.na(data.ztrans$bicarb), NA,(data.ztrans$bicarb  + 0.884331532)^2)
dist.delta$d4sq_bili    <- ifelse(is.na(data.ztrans$bili),   NA,(data.ztrans$bili    - 0.793065739)^2)
dist.delta$d4sq_bun     <- ifelse(is.na(data.ztrans$bun),    NA,(data.ztrans$bun     - 0.401287380)^2)
dist.delta$d4sq_cl      <- ifelse(is.na(data.ztrans$cl),     NA,(data.ztrans$cl      - 0.461395728)^2)
dist.delta$d4sq_creat   <- ifelse(is.na(data.ztrans$creat),  NA,(data.ztrans$creat   - 0.280646464)^2)
dist.delta$d4sq_crp     <- ifelse(is.na(data.ztrans$crp),    NA,(data.ztrans$crp     - 0.364269160)^2)
dist.delta$d4sq_elix    <- ifelse(is.na(data.ztrans$elix),   NA,(data.ztrans$elix    + 0.123710575)^2)
dist.delta$d4sq_esr     <- ifelse(is.na(data.ztrans$esr),    NA,(data.ztrans$esr     + 0.522607275)^2)
dist.delta$d4sq_gcs     <- ifelse(is.na(data.ztrans$gcs),    NA,(data.ztrans$gcs     + 0.761415458)^2)
dist.delta$d4sq_gluc    <- ifelse(is.na(data.ztrans$gluc),   NA,(data.ztrans$gluc    - 0.350033751)^2)
dist.delta$d4sq_hgb     <- ifelse(is.na(data.ztrans$hgb),    NA,(data.ztrans$hgb     + 0.055521451)^2)
dist.delta$d4sq_hr      <- ifelse(is.na(data.ztrans$hr),     NA,(data.ztrans$hr      - 0.490428737)^2)
dist.delta$d4sq_inr     <- ifelse(is.na(data.ztrans$inr),    NA,(data.ztrans$inr     - 0.785275737)^2)
dist.delta$d4sq_lactate <- ifelse(is.na(data.ztrans$lactate),NA,(data.ztrans$lactate - 1.092620482)^2)
dist.delta$d4sq_pao2    <- ifelse(is.na(data.ztrans$pao2),   NA,(data.ztrans$pao2    - 0.558641193)^2)
dist.delta$d4sq_plt     <- ifelse(is.na(data.ztrans$plt),    NA,(data.ztrans$plt     + 0.237985698)^2)
dist.delta$d4sq_rr      <- ifelse(is.na(data.ztrans$rr),     NA,(data.ztrans$rr      - 0.450954748)^2)
dist.delta$d4sq_sao2    <- ifelse(is.na(data.ztrans$sao2),   NA,(data.ztrans$sao2    - 0.011792492)^2)
dist.delta$d4sq_sex     <- ifelse(is.na(data.ztrans$sex),    NA,(data.ztrans$sex     - 0.107294741)^2)
dist.delta$d4sq_sodium  <- ifelse(is.na(data.ztrans$sodium), NA,(data.ztrans$sodium  - 0.232320267)^2)
dist.delta$d4sq_sbp     <- ifelse(is.na(data.ztrans$sbp),    NA,(data.ztrans$sbp     + 0.636731112)^2)
dist.delta$d4sq_temp    <- ifelse(is.na(data.ztrans$temp),   NA,(data.ztrans$temp    + 0.323962774)^2)
dist.delta$d4sq_trop    <- ifelse(is.na(data.ztrans$trop),   NA,(data.ztrans$trop    - 1.112482455)^2)
dist.delta$d4sq_wbc     <- ifelse(is.na(data.ztrans$wbc),    NA,(data.ztrans$wbc     - 0.323643150)^2)

#Calculate total distance to phenotype centers for each observation
dist.alpha$total <- sqrt(rowSums(dist.alpha[,2:n.features+1], na.rm = TRUE, dims = 1))
dist.beta$total  <- sqrt(rowSums(dist.beta[,2:n.features+1] , na.rm = TRUE, dims = 1))
dist.gamma$total <- sqrt(rowSums(dist.gamma[,2:n.features+1], na.rm = TRUE, dims = 1))
dist.delta$total <- sqrt(rowSums(dist.delta[,2:n.features+1], na.rm = TRUE, dims = 1))

#Add Distances to original data and select nearest center as phenotype
data.imputed$dist.alpha <- dist.alpha$total
data.imputed$dist.beta  <- dist.beta$total
data.imputed$dist.gamma <- dist.gamma$total
data.imputed$dist.delta <- dist.delta$total

data.imputed$phenotype <- ifelse(data.imputed$dist.alpha < data.imputed$dist.beta  & data.imputed$dist.alpha < data.imputed$dist.gamma & data.imputed$dist.alpha < data.imputed$dist.delta,"Alpha",
                          ifelse(data.imputed$dist.beta  < data.imputed$dist.alpha & data.imputed$dist.beta  < data.imputed$dist.gamma & data.imputed$dist.beta  < data.imputed$dist.delta,"Beta",
                          ifelse(data.imputed$dist.gamma < data.imputed$dist.alpha & data.imputed$dist.gamma < data.imputed$dist.beta  & data.imputed$dist.gamma < data.imputed$dist.delta,"Gamma",
                          ifelse(data.imputed$dist.delta < data.imputed$dist.alpha & data.imputed$dist.delta < data.imputed$dist.beta  & data.imputed$dist.delta < data.imputed$dist.gamma,"Delta",""))))

table(data.imputed$phenotype)
write.csv(data.imputed[,c(1,(n.features+2):(n.features+6))],"SENECA - Pheno Assignments.csv")


SENECA - Assigning Clinical Sepsis Phenotypes
#Jason Kennedy
#jnk28@pitt.edu
#Uploaded June 27, 2022

#Assign Clinical Sepsis Phenotypes (Alpha, Beta, Gamma, Delta) by Euclidean Distance
#Using SENECA Derivation Centroids
#Code presented in R

#From "Derivation, Validation, and Potential Treatment Implications of Novel Clinical Phenotypes for Sepsis"
#Seymour CW, Kennedy JN, Wang S, et al.
#JAMA. 2019 May 28;321(20):2003-2017
#DOI: 10.1001/jama.2019.5791
#PMID: 31104070; PMCID: PMC6537818

# NOTES #
#1.) Code assumes units are as follows:
#Age: years
#Albumin: g/dL
#ALT: U/L
#AST: U/L
#BANDS: %
#Bicarbonate: mEq/L
#Bilirubin: mg/dL
#BUN: mg/dL
#Chloride: mEq/L
#Creatinine: mg/dL
#C-Reactive Protein: mg/L
#Elixhauser: 0-31 point scale
#ESR: mm/hour
#GCS: 3-15 point scale
#Glucose: mg/dL
#Hemoglobin: g/dL
#Heart Rate: Beats/min
#INR: ratio
#Lactate (Serum)e: mmol/L
#PaO2: mmHg
#Platelets: x10^9/L
#Respiratory Rate: breaths/minute
#SaO2 - Oxygen Saturation: %
#Sex: 1=male; 0=female
#Sodium: mEq/L
#Systolic Blood Pressure: mmHg
#Temperature: Celcius
#Troponin: ng/mL
#White Blood Cell count: x10^9/L

#Convert units to above prior to phenotype assignment as needed

#2.) Expected input has subject identifier in 1st column and model features in remaining columns

#Required Libraries:
#missRanger

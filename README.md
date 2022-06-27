# SENECA - Assigning Clinical Sepsis Phenotypes

Jason Kennedy <br />
jnk28@pitt.edu <br />
Uploaded June 27, 2022

Assign Clinical Sepsis Phenotypes (Alpha, Beta, Gamma, Delta) by Euclidean Distance <br />
Using SENECA Derivation Centroids <br />
Code presented in R <br />

From "Derivation, Validation, and Potential Treatment Implications of Novel Clinical Phenotypes for Sepsis" <br />
Seymour CW, Kennedy JN, Wang S, et al. <br />
JAMA. 2019 May 28;321(20):2003-2017 <br />
DOI: 10.1001/jama.2019.5791 <br />
PMID: 31104070; PMCID: PMC6537818 <br />

# NOTES #
1.) Code assumes units are as follows: <br />
Age: years <br />
Albumin: g/dL <br />
ALT: U/L <br />
AST: U/L <br />
BANDS: % <br />
Bicarbonate: mEq/L <br />
Bilirubin: mg/dL <br />
BUN: mg/dL <br />
Chloride: mEq/L <br />
Creatinine: mg/dL <br />
C-Reactive Protein: mg/L <br />
Elixhauser: 0-31 point scale <br />
ESR: mm/hour <br />
GCS: 3-15 point scale <br />
Glucose: mg/dL <br />
Hemoglobin: g/dL <br />
Heart Rate: Beats/min <br />
INR: ratio <br />
Lactate (Serum)e: mmol/L <br />
PaO2: mmHg <br />
Platelets: x10^9/L <br />
Respiratory Rate: breaths/minute <br />
SaO2 - Oxygen Saturation: % <br />
Sex: 1=male; 0=female <br />
Sodium: mEq/L <br />
Systolic Blood Pressure: mmHg <br />
Temperature: Celcius <br />
Troponin: ng/mL <br />
White Blood Cell count: x10^9/L <br />

Convert units to above prior to phenotype assignment as needed

2.) Expected input has subject identifier in 1st column and model features in remaining columns <br />

Required Libraries: <br />
missRanger <br />

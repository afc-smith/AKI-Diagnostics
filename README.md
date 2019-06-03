# AKI-Diagnostics

## AKI-Diagnostics economic model

The files in this repository relate to the AKI Diagnostics economic model, which was developed as part of the National Institute of Health Research (NIHR) funded Health Technology Assessment (HTA) 'AKI Diagnostics' project. Full details of this project and the NIHR HTA monograph report (published in 2018) can be found at the project page: http://www.nets.nihr.ac.uk/projects/hta/1311613 (chapter 5 outlines the economic model structure and parameters). 

### Please cite this work as

> Hall PS, Mitchell ED, Smith AF, Cairns DA, Messenger M, Hutchinson M, Wright J, Vinall-Collier K, Corps C, Hamilton P, Meads D. The future for diagnostic tests of acute kidney injury in critical care: evidence synthesis, care pathway analysis and research prioritisation. Health Technology Assessment. 2018 4;22(32).

The AKI Diagnostics project was funded by the NIHR HTA programme, and supported by the NIHR Diagnostic Evidence Co-operative Leeds. The views expressed are those of the author(s) and not necessarily those of the NHS, the NIHR or the Department of Health.

### Licence for use

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.

The AKI Diagnostics economic model provided in this code was developed by Alison F Smith from the University of Leeds with supervision from Dr David Meads (Leeds) and Dr Peter Hall (University of Edinburgh). This code is provided under the Creative Commons Attribution Share Alike 4.0 International licence: please therefore appropriately cite the HTA monograph report (above), and this code repository if you wish to use this work. 

### Model Code Outline

The model is essentially run via the ‘Master’ R file. The Master file in turn reads in model parameters provided in the three ‘Parameters_’ files: 

1.	Parameters_HOSP.R  [parameters required to run the hospital period of the model]
2.	Parameters_FUP.R  [parameters required to run the follow-up period of the model, post hospital discharge]
3.	Parameters_FUP_tests.R  [parameters specific to the test intervention arms required for the follow-up period of the model, post hospital discharge]

The ‘Parameters_HOSP.R’ file reads in the following transition probability matrices for the model AKI cohort (patients who experience AKI in the ICU): 

1.	tp_norm.RData  
2.	tp_S1.RData
3.	tp_S2.RData
4.	tp_S3.RData
5.	tp_RRT.RData
6.	tp_ward.RData
7.	tp_ward_RRT.RData
8.	tp_ICUmort.RData
9.	tp_ICUdisch_alive.RData

These large data files are too large for github - you need to download them from here: https://www.dropbox.com/sh/lijf3ju0dslqkpp/AADeGN3kj_9BdHzyMoM-gXp4a?dl=0

Details of how these probabilistic transitions were derived using AKI registry data is provided in the HTA report. 
The Master file runs various arms of the model using the code provided in the seven ‘Model_’ files:

1.	Model_HOSP_baseline.R  [standard care baseline arm hospital period model; no intervention tests applied]
2.	Model_HOSP_notest.R  [Intervention arm hospital period model: no intervention test conducted due to patients having pre-existing non-early AKI] 
3.	Model_HOSP_test_noAKI.R  [Intervention arm hospital period model: intervention test applied, no AKI cohort – including true negative and false positive cases]
4.	Model_HOSP_test_FN.R  [Intervention arm hospital period model: intervention test applied, false negative cases]
5.	Model_HOSP_test_TP_treat.R  [Intervention arm hospital period model: intervention test applied, true positive early AKI cases that can benefit from early AKI treatment]
6.	Model_HOSP_test_TP_notreat.R  [Intervention arm hospital period model: intervention test applied, true positive non-early AKI cases that cannot benefit from early AKI treatment]
7.	Model_FUP.R  [Follow up period of the model, post hospital discharge]



# HighFlowSorption
Data analysis and figure generation for King et al. "River phosphorus cycling during high flow constrains Lake Erie cyanobacteria blooms"

Prepared by JMH (March 2022)

------------- 00_wqBackup -------------
Contains original data downloaded from USGS or Heidelberg, but otherwise not manipulated
"blancharddata.xlsx"
"lostcreekdata.xlsx"    
"maumeedata.xlsx"       
"SouthTurkeyFoot.xlsx" 
"westdata.xlsx"         
"zz_turkeydataOLD.xlsx" DELETE

------------- 01_RawData -------------
"00_WK_IsoAmb_JMH.csv"
"00wq_blancharddata.csv": converted  to csv from excel file with same name, no other changes            
"00wq_blancharddata.xlsx": downloaded data with columns adapted for R 
"00wq_lostcreekdata.csv": converted  to csv from excel file with same name, no other changes                        
"00wq_lostcreekdata.xlsx": downloaded data with columns adapted for R
"00wq_maumeedata.csv": converted  to csv from excel file with same name, no other changes                            
"00wq_maumeedata.xlsx": downloaded data with columns adapted for R
"00wq_turkeydata.csv": converted  to csv from excel file with same name, no other changes                            
"00wq_westdata.csv": converted  to csv from excel file with same name, no other changes            
"00wq_westdata.xlsx"    : downloaded data with columns adapted for R              
"MaumeeWatervilleDailyQ1975to2019.csv"
"RawSorptionData.csv": Results from sorption rate and isotherm measurements                
"ScalingValTable.csv": USGS gage #, lat, long, watershed km2
"SedExp_R_WK.csv": Raw sorption and isotherm data                      
"SlopesClean20201015jmh.csv": Slopes and segment lengths for reaches between tributary gage and Waterville, OH	
"SouthTurkeyFoot.xlsx": downloaded data with columns adapted for R                
"zz_00wq_turkeydata_OLD.csv" DELETE
"zz_00wq_turkeydataOLD.xlsx"  DELETE


------------- 02a_ThesisAnalyses - EMPTY delete -------------

------------- 02b_SorptionScaling -------------
"00_EPCredo.R": Fits Langmuir model to isotherm data and merges with other isotherm data.
"01_SorpDataMungeFigS4_S5.R": Cleans up sorption data for scaling. Makes Fig. S4 and S5    

The following three scripts: 
(1) Downloads Q data for target stream and for downstream monitoring stations if available;
(2)Combines and processes Q and solute chem data; 
(3) estimates 75% flow target stream; 
(4) Builds models to predict missing SRP and SS data, then predicts missing values; 
(5) Generates final stream Q, SRP, SS csv for STF            
"02a_WQandQ_STF.R"
"02b_WQandQ_UTLC.R"                    
 "02c_WQandQ_WC.R"
 
"02d_WQandQ_MaumeeWatervilleTake.R": GEts maumee river discharge and water quality data. Fills in gaps for SRP, SS, and TP load and concentration. Combines all data then calculates daily loads. Exports daily load data
    
"03_CleanSegSlopes.R": loads and cleans slope data which was measured by WK in QGIS

"04_ScalingUp_TribsMaumee.R": Assuming colloid-P:DRP = 0.5, combines water quality data for three tributaries and Maumee River and estimates: (1) quantity of P sorption between tributary gaging station and Waterville, OH (Fig 2) and (2) the P sorption upstream of Waterville, OH (Fig. 3). Conducts bootstrapping for each tributary and Maumee River and exports the bootstrap results to a csv.

 "04b_ScalingUp_TribsMaumee_50per.R": Assuming colloid-P:DRP = 0.5, combines water quality data for three tributaries and Maumee River and estimates: (1) quantity of P sorption between tributary gaging station and Waterville, OH (Fig 2) and (2) the P sorption upstream of Waterville, OH (Fig. 3). Conducts bootstrapping for each tributary and Maumee River and exports the bootstrap results to a csv.
 
"05_TribFig2.R": Generates colloid-P:DRP = 0 data for Fig. 2. Also calculates some values used in Results section and will generate a version of Fig. 2 just for colloid-P:DRP = 0.                        
"05b_TribFig2_50perDRP.R": Generates colloid-P:DRP = 0.5 data for Fig. 2. Also calculates some values used in Results section and will generate a version of Fig. 2 just for colloid-P:DRP = 0.

"05b_TribFig2_combined.R": Generates Fig. 2              

"06_MaumeeBootMungePcyanos.R": Estimates cyanoHABs with and without P sorption. Assumes colloid-P:DRP = 0

"06b_MaumeeBootMungePcyanos_50perDRP.R": Estimates cyanoHABs with and without P sorption. Assumes colloid-P:DRP = 0.5

"07_FigS6.R": Generates Figure S6, a version of Fig. 3 assuming colloid-P:DRP = 0

"07b_FigS7_50perDRP.R": Generates Figure S6, a version of Fig. 3 assuming colloid-P:DRP = 0.5                 
"07c_Fig3_combined.R": Generates Fig 3

"08_MaumeeWatervillePloadsFigS1.R": Generates Fig S1

"09_SStimeseriesFigS8.R": Generates Fig S8 (changes in discharge and suspended sediment load during 1975-2019). Also makes some calculations for the Results.

"11a_Fig4_comp.R": Generates Fig 4

"11b_FigS9_Table1_0perDRP.R": Generates Fig. S9 and does analyses for Fig. 4 and Table 1

"11c_FigS10_TableS3_50perDRP.R":  Generates Fig. S10 and does analyses for Fig. 4 and Table S3        
"12_FigS11and12_kd.R": Partitioning analysis (kd) and generates Fig. S11 and S12

------------- 03_Rdata -------------
Rdata associated with script/analyses 

"00_EPCredo_Rdat"
"01_SorpDataMungeFig_rdat"
"02_Turk_Rdat"                                 
"02_West_Rdat"
"02a_CombUTLC_Rdat"
"02d_MaumeeWatervilleWQandQ_Rdat"              
"02d_WQandQ_MaumeeWatervilleTake3_Rdat"
"05_TribFigs_50perDRP_rdat"
"05_TribFigs_rdat"                             
"06_MaumeeBootMungePcyanos_Rdat"
"06b_MaumeeBootMungePcyanos_50perDRP_Rdat"
"06c_ScalingSorption2TribMaumee_50perDRP_Rdata"
"06c_ScalingSorption2TribMaumee_Rdata"
"07_FigS6_Rdata"
"07b_Fig4_50perDRP_Rdata"                      
"07c_Fig3_comp_Rdat"
"08_MaumeeWatervillePloadsFigS1_Rdata"
"09_SStimeseriesFigS6_Rdata"                   
"10_DRP_TPtimeseries"
"11a_Fig4_comp_Rdat"
"11b_FigS9_Table1_0perDRP_Rdat"                
"11c_FigS10_TableS3_50perDRP_Rdat"
"12_FigS10and11_Kd_Rdat"  

------------- 04a_generatedDataOnGit -------------
CSV files generated by script that are small enough to upload to GitHub. 
Prefix number indicates the Rscript that generated the CSV

"01_CleanedSorpIsoDat.csv"
"02d_MaumeeWatervilleDischarge.csv"
"02d_MaumeeWatervilleWaterQual.csv"       
"02d_UTLC.csv"
"02f_Turk.csv""02f_West.csv"                            
"03d_slopes.csv"
"04_WQandTT_STF.csv"
"04_WQandTT_UTLC.csv"                     
"04_WQandTT_WC.csv"
"04c_WQandTT_STF.csv"
"04c_WQandTT_UTLC.csv"                    
"04c_WQandTT_WC.csv"
"04d_SorpDat.csv"
"05_TribDat4Fig2a_0col.csv"               
"05_TribDat4Fig2bc_0col.csv"
"05b_TribDat4Fig2a_50col.csv"
"05b_TribDat4Fig2bc_50col.csv"            
"06_Scale2WVall_MjHf_cyanos.csv"
"06_Scale2WVall_MjHf.csv"
"06b_Scale2WVall_MjHf_50perDRP.csv"       
"06b_Scale2WVall_MjHf_cyanos_50perDRP.csv"
"07_Fig3Dat_maum_0col.csv"
"07_Fig3Dat_maum_50col.csv"               
"07_Fig3Dat_maumCyano_0col.csv"
"07_Fig3Dat_maumCyano_50col.csv"
"09_MaumeeLoadsMarJuneHighFlow.csv"       
"09_SS_QQandYglm"
"11a_Fig4dat_50perDRP_G03.csv"
"11a_Fig4dat_50perDRP_L03.csv"            
"11a_Fig4dat_G03.csv"
"11a_Fig4dat_L03.csv"

------------- 04a_generatedDataTooBigForGit -------------
CSV files that are too large to upload to GitHub.
These are all the bootstrap results
These are available from JMH upon request.
Prefix number indicates the Rscript that generated the CSV

"04_RawBootstrapResults_Maumee.csv"
"04_TribRawBootstrapResults_STF.csv"
"04_TribRawBootstrapResults_UTLC.csv"          
"04_TribRawBootstrapResults_WC.csv"
"04b_RawBootstrapResults_Maumee_50perDRP.csv"
"04b_TribRawBootstrapResults_STF_50perDRP.csv" 
"04b_TribRawBootstrapResults_UTLC_50perDRP.csv"
"04b_TribRawBootstrapResults_WC_50perDRP.csv"
"04c_TribRawBootstrapResults_STF.csv"          
"04c_TribRawBootstrapResults_UTLC.csv"
"04c_TribRawBootstrapResults_WC.csv"    

------------- 05_Figures -------------
Figures
Prefix number indicates the Rscript that generated the CSV

"01_FigS4.png"
"01_FigS4AgTalk.png"
"01_FigS5.png"
"05_Fig2_50perDRP.png"            
"05_Fig2_combined.png"
"05_Fig2.png"
"07_Fig3.png"
"07_FigS6.png"                    
"07_GraphicalAbstractCyanoFig.png"
"07b_FigS7_50perDRP.png"
"08_DRPtoTPratio.png"
"08_FigS1.png"                    
"09_FigS8.png"
"11a_Fig4.png"
"11b_FigS9.png"
"11c_Fig4.png"                    
"11c_FigS10.png"
"12_FigS10.png"
"12_FigS11.png"    


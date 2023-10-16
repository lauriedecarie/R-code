Clear_all = T

if (Clear_all) 
{
  closeAllConnections()
  rm(list=ls())
  graphics.off()
  cat("\f")
} 

###############PROJET Laurie_Nick############

library(rms)        ; library(rms)
library(bigmemory)  ; library(bigmemory)
library(nlme)       ; library(nlme)
library(varhandle)  ; library(varhandle)
library(writexl)    ; library(writexl)
library(readxl)     ; library(readxl)
library(DescTools)  ; library(DescTools)
library(xtable)     ; library(xtable)
library(stringr)    ; library(stringr)
library(foreign)    ; library(foreign)
library(Hmisc)      ; library(Hmisc)
library(purrr)      ; library(purrr)
library(rlist)      ; library(rlist)
library(snpar)      ; library(snpar)
library(gridExtra)  ;library(gridExtra)
library(lme4)       ;library(lme4)
library(ClusterR)    ; library(ClusterR)
library(ClustImpute)  ; library(ClustImpute)
cat("\f")

################
# C:\POST_DOC\CIUSSS\CodeR3_Samira2
# Path_E   = "D:/CIMAQ/RegressionQuatile/CodeR3_SamiraV03" # A MODIFIER AVANT DE ROULER
# Path_F   = "D:/CIMAQ/RegressionQuatile/CodeR3_SamiraV03/FUNCTIONS2"  # A MODIFIER PATH
Path_E   = "C:/POST_DOC/CIUSSS2/Labo_Sylvie_Belleville/Projet_Laurie_Nick/Analyse" # A MODIFIER AVANT DE ROULER
Path_F   = "C:/POST_DOC/CIUSSS2/Labo_Sylvie_Belleville/Projet_Laurie_Nick/Analyse/FUNCTIONS2"  # A MODIFIER PATH

setwd(Path_F)
#---------------------------------------#
source("UnFactorF.R") ;source("AsFactorF.R")
source("Extract_VARF.R") ;source("Summary_F.R")
source("Read_ExcelF.R") ; source("Recode_VectorF.R")
source("MKDIRSF.R")  ; source("SELECT_Var_p_F.R")
source("Keep_NUMERIC2_F.R") ; source("Forma_TABLEF.R")
source("Remove_SpaceF.R") ; source("CapitalizeF.R")
source("Reorder_LevelsF.R") ; source("DichomizeF.R")
source("UnFactorF.R") ;source("AsFactorF.R")
source("MKDIRSF.R")  ; source("SELECT_Var_p_F.R")
source("Age_GROUPS_F.R") ; source("Stat_DESC_VQUAL_F.R")
source("As.vectorF.R") ; source("FILE_PATHF.R")
source("As.data.frameF.R") ; source("DichomizeF.R")
#---------------------------------------#
setwd(Path_E)


######################################
######### PATH ###################
######################################

PATH_R1   = "RESULTATS/R1/" # FOR UNIVARIATE AND BIVARIATE ANALYSIS
MKDIRS(PATH_R1) 
PATH_R3   = "RESULTATS/R3/" # FOR MULTIVARIATE ANALYSIS
MKDIRS(PATH_R3) 


######################################
######### DATABASE ###################
######################################
# cat("\f") 

Simul_DATABASE <- F

if (Simul_DATABASE){
  source("01_Code01_A1.R") # FOR RECODAGE UNIVARIATE AND BIVARIATE ANALYSIS
}

#########
# cat("\f") ; View(DATABASE)

PATH_STATISTIC.DESCRIPTIVE  = paste0(PATH_R1,"01.STATISTIC.DESCRIPTIVE/")
load(file=paste0(PATH_STATISTIC.DESCRIPTIVE,paste0("VARIABLES",".RData")))
load(file=paste0(PATH_STATISTIC.DESCRIPTIVE,paste0("DATABASE_BDF",".RData")))
load(file=paste0(PATH_STATISTIC.DESCRIPTIVE,paste0("Stat_charF",".RData")))
load(file=paste0(PATH_STATISTIC.DESCRIPTIVE,paste0("Stat_numF",".RData")))
Var_Names    = Var_Names ; Var_list=Var_list ; Var_lab_list=Var_lab_list
DATABASE     = DATABASE_BDF ; Stat_charF = Stat_char0F ; Stat_numF=Stat_num0F
Var_CharF    = Stat_charF$Names ; Var_numF    = Stat_numF$Names

###################################
######### COVARIABLES #############
###################################

Covar_list    = list()
Covariables   = Var_lab_list$Covariables
Diagnostics   = Covariables[which(grepl("Diagnostic",Covariables, fixed = TRUE))]
Covar_others  = Covariables[-which(grepl("Diagnostic",Covariables, fixed = TRUE))]
NCovar_others = length(Covar_others) ; NDiagnostics = length(Diagnostics)
NCovar_others_seq = seq(NCovar_others)

Simul_Covar_list = F
if (Simul_Covar_list) 
{
  k = 0
  for (i in 1:NDiagnostics){k = k+1 ; Covar_list[[k]] = Diagnostics[i]}
  for (i in 1:NDiagnostics){
    Diagnostic0   = Diagnostics[i]
    for (j in 1:NCovar_others){
      jIndex_comb0  = combinations(NCovar_others,j)
      NjIndex_comb0 = nrow(jIndex_comb0)
      for (h in 1:NjIndex_comb0){
        jIndex0     = jIndex_comb0[h,]
        k = k+1 ; Covar_list[[k]] = c(Diagnostic0,sort(Covar_others[jIndex0]))
      }
    }
  }
  save(Covar_list,file="Covar_list.RData")
}

if (!Simul_Covar_list){load(file="Covar_list.RData")}


###################################
######### VARIABLES ###############
###################################

Var_YY       = Var_lab_list$Activation
Var_XX       = Var_lab_list$Associative.memory
Covar_list   = Covar_list
NmaxCovar    = max(unlist(lapply(Covar_list,function(w) length(w))))
NCovar_list  = length(Covar_list)

################################
######### MODELS ###############
################################

Model_M1   = list(Model = "M1",
                  Model_names = "Activation.Hippocampus.Left(Atrophie_Left_hip)",
                  Names = "Hippocampe.activation/Hippocampal.volume",
                  Ryx = "Linear negative",
                  Y=Var_YY[1],X=Var_XX[1])
Model_M2   = list(Model = "M2",
                  Model_names = "Activation.Temporal.Inferior.Right(Atrophie_Left_hip)",
                  Names = "Temporal.activation/Hippocampal.volume",
                  Ryx = "Linear negative",
                  Y=Var_YY[2],X=Var_XX[1])
Model_M31  = list(Model = "M31",
                  Model_names = "Activation.Parietal.Superior.Left(cortical_thikness)",
                  Names = "Parietal.activation/Thickness",
                  Ryx   = "Quadratic concave",
                  Y=Var_YY[3],X=Var_XX[2])
Model_M32  = list(Model = "M32",
                  Model_names = "Activation.Parietal.Superior.Left(memory_performance)",
                  Names = "Parietal.activation/Associative.memory",
                  Ryx   = "Quadratic concave",
                  Y=Var_YY[3],X=Var_XX[3])
Model_M33  = list(Model = "M33", 
                  Model_names = "Activation.Parietal.Superior.Left(Atrophie_Left_hip)",
                  Names = "Parietal_activation/Hippocampal_Volume",
                  Ryx = "Quadratic concave",
                  Y=Var_YY[3],X=Var_XX[1])

MODELS     = list(Model_M1=Model_M1,Model_M2=Model_M2,
                  Model_M31=Model_M31,Model_M32=Model_M32,Model_M33=Model_M33)
NMODELS    = length(MODELS)

#########
# Pch.v         = c(20,2,3,4,5,8)
# Pch.v2        = c("●","???","+","?","???","*")

Diagnostics.v = c("Prodromal_AD","HC+","HC","SCD","SCD+","MCI")
Colors.v      = c("blue", "gold", "blueviolet", "chocolate4","gray49","mediumturquoise")
Pch.v         = rep(20,6)
Pch.v2        = rep("●",6)
Param.diagnostics = data.frame(Diagnostics= Diagnostics.v,Colors= Colors.v,
                              Pch=Pch.v,Pch2=Pch.v2)

#####################################
######### COVARIABLES SELECTION #####
#####################################

Cov_retain      = c("Age","Sex")
Covar_list_used = list()
k0              = 0

for (k in 1:length(Covar_list)) 
{
  # cat("\f") ; k = 10
  Covar_list0 = Covar_list[[k]]
  test_1      = (length(Covar_list0)%in%c(1,3))
  test_2      = (all(Cov_retain%in%Covar_list0)==TRUE)
  ########
  if (length(Covar_list0)==1){
    if (test_1){k0=k0+1;Covar_list_used[[k0]]=Covar_list0}
  }
  ########
  if (length(Covar_list0)>1){
    if (test_1&test_2){k0=k0+1;Covar_list_used[[k0]]=Covar_list0}
  }
}


#####################################
######### REGRESSION.QR.PLOT ########
#####################################

PATH_REGRESSION.QR.PLOT  = paste0(PATH_R3,"03.REGRESSION.QR.PLOT/")
MKDIRS(PATH_REGRESSION.QR.PLOT)
NCov  = length(Covar_list_used)

###############################################################################
#######Seulement les participants qui ont un SCD+ (triangle) et un MCI#########
# cat("\f"); View(Data.Diagnostics)

Diagnostics.vars = paste0("Diagnostic.",seq(6))
Data.Diagnostics = DATABASE[,Diagnostics.vars]
Diagnostics.v    = sort(unique(na.omit(as.vector(unlist(as.vector.data.frame(DATABASE[,Diagnostics.vars]))))))

index.Diagnostics = c()

for (iCov in 1:NCov)
{
  # cat("\f") ; iCov = 10
  Covar_X0 = Covar_list_used[[iCov]] ; Covar_X0_init = Covar_X0
  Covar_X0 = Covar_X0[which(grepl("Diagnostic", Covar_X0, fixed = TRUE))]
  Data_Covar_X0 = DATABASE[,Covar_X0]
  stopifnot(is.factor(Data_Covar_X0)==T)
  #######
  With_other = F
  if(("Age"%in%Covar_X0_init)&("Sex"%in%Covar_X0_init)){With_other=T}
  #######
  
  # if((With_other==T)&("HC"%in%levels(Data_Covar_X0))&("MCI"%in%levels(Data_Covar_X0))&("SCD+"%in%levels(Data_Covar_X0))){
  #  index.Diagnostics=c(index.Diagnostics,iCov)
  #  break
  # }
  
  #######
  
  if((With_other==T)&(length(levels(Data_Covar_X0))==1)&("Prodromal_AD"%in%levels(Data_Covar_X0))){
    index.Diagnostics=c(index.Diagnostics,iCov)
   break
  }
  
}




##################################  
#########REGRESSION###############


for (iCov in index.Diagnostics) # 1:NCov
{
  # cat("\f") ; iCov = 7
  Covar_X0 = Covar_list_used[[iCov]]
  
  for (iM in 1:NMODELS) # NMODELS
  {
    # cat("\f") ; iM = 1
    #########
    reV.x.Axis0 = T
    MODELS0  = MODELS[[iM]]
    if (MODELS0$Model=="M32"){reV.x.Axis0 = F}
    ########
    Negative.x.Axis0 = F
    if (MODELS0$Model=="M32"){Negative.x.Axis0 = T}
    
    ########
    Model0   = MODELS0$Model ; Model_names0 = MODELS0$Model_names
    Names0   = MODELS0$Names ; Rxy = MODELS0$Ryx
    Var_Y    = MODELS0$Y     ; Var_X   = MODELS0$X ; Var_X2 = paste0(Var_X,"2")
    Covar_X  = Covar_X0
    if ("Diagnostic.1"%in%Covar_X0){
      Covar_X = Covar_X0[!(Covar_X0%in%"Diagnostic.1")]
    }
    
    Var_Model0 = c("ID",Var_Y,Var_X,Var_X2,Covar_X)
    #########
    DATA_REG0 = Extract_VAR(DATABASE,c("ID",Var_Y,Var_X,Covar_X))
    DATA_REG0[,Var_X2] = (DATABASE[,Var_X])^2
    DATA_REG0 = DATA_REG0[,Var_Model0]
    #########
    DATA_REG0$Diagnostic.5 = DATABASE$Diagnostic.5
    
    #########
    Diagnostic0  = Covar_X0[which(grepl("Diagnostic", Covar_X0, fixed = TRUE))]
    DATA_REG0    = DATA_REG0[!is.na(DATABASE[,Diagnostic0]),]
    
    #########
    CoVar_X_num0 = NULL ; CoVar_X_cat0 = NULL
    if (length(intersect(Var_numF,Covar_X))>0){CoVar_X_num0=intersect(Var_numF,Covar_X)} 
    if (length(intersect(Var_CharF,Covar_X))>0){CoVar_X_cat0=intersect(Var_CharF,Covar_X)}
    Var_XF0     = c(Var_X,Var_X2)
    if (length(which(grepl("Linear", Rxy, fixed = TRUE)))>0){Var_XF0 = Var_X}
    #########
    Covar_X0_save      = paste(Covar_X0,collapse = "-")
    Path.Model0 = paste0(PATH_REGRESSION.QR.PLOT,"GRAPHS","/",Model0,"-",Model_names0,"/")
    MKDIRS(Path.Model0)
    Path.graph0   = paste0(Path.Model0,paste0(Model0,"-",Covar_X0_save,".tif"))
    # Path.graph0   = paste0(Path.Model0,paste0(Model0,"-",Covar_X0_save,".png"))
    Diagnostic0  = Covar_X0[which(grepl("Diagnostic", Covar_X0, fixed = TRUE))]
    Diagnostic0F = gsub("[.]","-",Diagnostic0)
    Diagnostic0_txt = paste(levels(DATABASE[,Diagnostic0]),collapse = "-")
    Title0      = paste0(Model0,"-",Names0," (",Diagnostic0," : ",Diagnostic0_txt,")")
    
    #########
    
    if (!(Diagnostic0%in%names(DATA_REG0))){
      if ("Diagnostic.1"%in%Covar_X0){
        DATA_REG0[,Diagnostic0] = DATABASE[!is.na(DATABASE$Diagnostic.1),Diagnostic0]
      }
      if (!("Diagnostic.1"%in%Covar_X0)){
        DATA_REG0[,Diagnostic0] = DATABASE[,Diagnostic0]
      }
    }
    
    #########
    # Factors.legend0     = levels(DATABASE[,Diagnostic0])
    
    Factors.legend0     = levels(DATABASE[,"Diagnostic.5"])
    Param.diagnostics0  = Param.diagnostics[Param.diagnostics$Diagnostics%in%Factors.legend0,]
    Data.diagnostics0   = data.frame(Diagnostics=Factors.legend0)
    Param.diagnostics0F = merge(x=Data.diagnostics0, y=Param.diagnostics0,
                                by = "Diagnostics",sort = FALSE)
    #########
    
    DATA_REGG          = DATA_REG0
    activation_HC      = DATABASE[(DATABASE$Diagnostic.3=="HC"),Var_Y]
    Mean_activation_HC = mean(activation_HC,na.rm=T) 
    Sd_activation_HC   = sd(activation_HC,na.rm=T)
    n_activation_HC    = length(activation_HC[!is.na(activation_HC)])
    z_alpha            = 1.96
    Mean_activation_HC_IC_low = Mean_activation_HC-z_alpha*Sd_activation_HC/((n_activation_HC)^0.5)
    Mean_activation_HC_IC_up  = Mean_activation_HC+z_alpha*Sd_activation_HC/((n_activation_HC)^0.5)
    
    ##########
    Data_activation_HC = data.frame(Mean_activation_HC=Mean_activation_HC,
                                    IC_low=Mean_activation_HC_IC_low,
                                    IC_up=Mean_activation_HC_IC_up)
    ###########
    # cat("\f"); View(DATA_REGG_F)
    
    Data_Diagnostic0 = DATA_REGG[,Diagnostic0]
    DATA_REGG_F      = DATA_REGG[which(Data_Diagnostic0%in%c("MCI","SCD+")),]
    
    #############
   
    source("REGRESSION_Q_F.R")
    REGQ0      = REGRESSION_Q(Data           = DATA_REGG ,
                              Var_ID         = "ID",
                              Var_y          = Var_Y,
                              Var_X          = Var_XF0,
                              CoVar_X_num    = CoVar_X_num0 ,
                              CoVar_X_cat    = CoVar_X_cat0 ,
                              Cov_effects    ="Linear",
                              Var_X_NULL     = F ,
                              auto.weight    = T , 
                              relative.error = T ,
                              ImputeNA       = T ,
                              Tau_c          = c(0.25,0.50,0.75),
                              add.plot       = T ,
                              add.Quantreg   = T ,
                              Path.graph     = Path.graph0,
                              rev.x.Axis     = reV.x.Axis0, # = F si pas renverser
                              Negative.x.Axis = Negative.x.Axis0, # = T si multiplier par -1
                              param.plot     = list(Title=Title0,Var_factors="Diagnostic.5",
                                                    x_lab =NULL,y_lab=NULL,add.intercept=T,
                                                    Factors.names="Diagnostic patients",bdw=0.9,
                                                    Factors.legend.as = Param.diagnostics0F$Diagnostics,
                                                    Factors.colors.as = Param.diagnostics0F$Colors,
                                                    Factors.pch.as    = Param.diagnostics0F$Pch,
                                                    Factors.pch2.as   = Param.diagnostics0F$Pch2,
                                                    plot.ClusterQR    = F,
                                                    Add_OLS = F))# Var_factors=Diagnostic0
    
    REGQ0$Data_plot$Mean_activation_HC = Mean_activation_HC
    
    ##################
    REG_QR0   = REGQ0$OUTPUT_QR$REG_QR_TAB
    names(REG_QR0)[1] = paste0(Model0,"-",Var_Y," (",gsub("[.]"," ",Diagnostic0),")")
    REG_QR0_f = Forma_TABLE(REG_QR0,type.format=2,
                            Title = paste0("Quantile regression : ",
                                           Model0,"-",Names0,
                                           " (",Diagnostic0," : ",Diagnostic0_txt,")"))
    ##################
    
    IQQ_TEST0   = REGQ0$OUTPUT_QR$INTER_QQ_TEST
    names(IQQ_TEST0)[1] = paste0(Model0,"-",Var_Y," (",gsub("[.]"," ",Diagnostic0),")")
    names(IQQ_TEST0) = gsub("Q50","Qmean",names(IQQ_TEST0))
    IQQ_TEST0_f = Forma_TABLE(IQQ_TEST0,type.format=2,
                              Title = paste0("Inter-Quantile regression test : ",
                                             Model0,"-",Names0,
                                             " (",Diagnostic0," : ",Diagnostic0_txt,")"))
    
    
    ##################
    Path.REG.Model0    = paste0(PATH_REGRESSION.QR.PLOT,"EQUATIONS","/",Model0,"-",Model_names0,"/")
    MKDIRS(Path.REG.Model0)
    REG_QR0_f %>% cat(., file = paste0(Path.REG.Model0,paste0(Model0,"-",Covar_X0_save,".html")))
    write_xlsx(REG_QR0, paste0(Path.REG.Model0,paste0(Model0,"-",Covar_X0_save,".xlsx")))
    # save(REG_QR0,file = paste0(Path.REG.Model0,paste0(Model0,"-",Covar_X0_save,".RData")))
    ##################
    Path.DATA.PLOT.Model0    = paste0(PATH_REGRESSION.QR.PLOT,"DATA.PLOT","/",Model0,"-",Model_names0,"/")
    MKDIRS(Path.DATA.PLOT.Model0)
    write_xlsx(REGQ0$Data_plot, paste0(Path.DATA.PLOT.Model0,paste0(Model0,"-",Covar_X0_save,".xlsx")))
    
    ##################
    Path.INTER_QUANTILES_TEST.Model0    = paste0(PATH_REGRESSION.QR.PLOT,"INTER_QUANTILES_TEST","/",Model0,"-",Model_names0,"/")
    MKDIRS(Path.INTER_QUANTILES_TEST.Model0)
    IQQ_TEST0_f %>% cat(., file = paste0(Path.INTER_QUANTILES_TEST.Model0,paste0(Model0,"-",Covar_X0_save,".html")))
    write_xlsx(IQQ_TEST0,paste0(Path.INTER_QUANTILES_TEST.Model0,paste0(Model0,"-",Covar_X0_save,".xlsx")))
    
    ##################
    Path.CLUSTER_QR.Model0    = paste0(PATH_REGRESSION.QR.PLOT,"CLUSTER_QR","/",Model0,"-",Model_names0,"/")
    MKDIRS(Path.CLUSTER_QR.Model0)
    write_xlsx(REGQ0$OUTPUT_QR$CLUSTER_QR, paste0(Path.CLUSTER_QR.Model0,paste0(Model0,"-",Covar_X0_save,".xlsx")))
    
    ##################
    Path.DATA.MEAN.ACT.Model0    = paste0(PATH_REGRESSION.QR.PLOT,"DATA.MEAN.ACT.HC","/",Model0,"-",Model_names0,"/")
    MKDIRS(Path.DATA.MEAN.ACT.Model0)
    write_xlsx(Data_activation_HC, paste0(Path.DATA.MEAN.ACT.Model0,paste0(Model0,"-",Covar_X0_save,".xlsx")))
  }
}



cat("\f")









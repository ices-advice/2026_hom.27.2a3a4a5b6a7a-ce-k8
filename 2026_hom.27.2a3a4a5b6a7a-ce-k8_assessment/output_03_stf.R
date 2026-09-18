
library(icesTAF)
library(r4ss)
library(dplyr)
library(tidyverse)
#*****************************
#Reference Points
#************************
MSYBtrigger = 787443
Fmsy = 0.08
Blim=566678
Bpa= MSYBtrigger
Flim=0.22
Fp05=0.08
Fpa=Fp05 

################## Directories

stf.dir <- paste0(getwd(),"/model/stf")
scenarios.dir <- paste0(getwd(),"/model/stf/Scenarios")
SS3_exe_path <- paste0(getwd(),"/boot/software/ss3.exe")
mkdir("output/stf")
stf_out<-paste0(getwd(),"/output/stf")


################# Get the outputs from the stf
#retrieve and summarise outputs

################# Get the outputs
#retrieve and summarise outputs

FMult_names <- c("Fmsy","Fmsy","F0","Fsq","Btrigger3","Blim10")


forecastModels <- r4ss::SSgetoutput(dirvec = file.path(stf.dir,"Scenarios", c(FMult_names)), getcovar = FALSE)
save(forecastModels,file=file.path(scenarios.dir, "forecast.RData"))
forecastSummary <- SSsummarize(forecastModels)

#SSB, F, Recr
SSB <- as.data.frame(forecastSummary$SpawnBio)
SSB.SD <- as.data.frame(forecastSummary$SpawnBioSD)   #SD
SSB.Lower <- as.data.frame(forecastSummary$SpawnBioLower)   #SSB + 1.96SD
SSB.Upper <- as.data.frame(forecastSummary$SpawnBioUpper)   #SSB - 1.96SD
Fvalue <- as.data.frame(forecastSummary$Fvalue)
Recr <- as.data.frame(forecastSummary$recruits)

#summarise results
all.scen <-FMult_names
num.scen <- length(all.scen)
dfSummary <- 
  data.frame(Scenario.Type = c("Fmsy","Fpa","F0","Fsq","Btrigger","Blim"),
             Val = all.scen,
             Catch_2027 = NA, Catch_2028 = NA, 
             F_2027 = as.numeric(Fvalue[Fvalue$Yr==2027,paste0('replist',seq(1,num.scen))]),
             F_2028 = as.numeric(Fvalue[Fvalue$Yr==2028,paste0('replist',seq(1,num.scen))]),
             SSB_2027 = as.numeric(SSB[SSB$Yr==2027,paste0('replist',seq(1,num.scen))]),
             SSB_2027_SD = as.numeric(SSB.SD[SSB.SD$Yr==2027,paste0('replist',seq(1,num.scen))]),
             SSB_2028 = as.numeric(SSB[SSB$Yr==2028,paste0('replist',seq(1,num.scen))]),
             SSB_2028_SD = as.numeric(SSB.SD[SSB.SD$Yr==2028,paste0('replist',seq(1,num.scen))]),
             Rec_2027 = as.numeric(Recr[Recr$Yr==2027,paste0('replist',seq(1,num.scen))]),
             Rec_2028 = as.numeric(Recr[Recr$Yr==2028,paste0('replist',seq(1,num.scen))]),
             Fmsy=Fmsy, Blim=Blim, Change_advice=NA,Change_ssb=NA)

#probability of being below Blim
dfSummary$pBlim_2028 <- pnorm(Blim, dfSummary$SSB_2028, dfSummary$SSB_2028_SD)

#catches
for (i in 1:num.scen){
  output = forecastModels[[i]]
  catch <- output$timeseries %>% filter(Era == "FORE" ) %>% 
    dplyr::select("Yr", starts_with("dead(B)")) %>% 
    mutate(FMult = FMult_names[i])
  names(catch) <- c("Year","Catch","FMult")
  catch <- catch %>% pivot_wider(names_from = Year, values_from = Catch, names_prefix = "Catch_")
  dfSummary$Catch_2027[i] <- catch$Catch_2027
  dfSummary$Catch_2028[i] <- catch$Catch_2028
}

## Change in advice and ssb
last_adv=74214

for (i in 1:num.scen){
  
  dfSummary$Change_ssb[i]<-(dfSummary$SSB_2028[i]-dfSummary$SSB_2027[i])/dfSummary$SSB_2027[i]*100
  dfSummary$Change_advice[i]<-(dfSummary$Catch_2027[i]-last_adv)/last_adv*100
}

write.csv(dfSummary, file = paste0(stf_out,"/catch_option_table.csv"))


##########
#Stock development included the forecast period for the FMSY scenario
##########

stf_dir <- paste0(stf.dir,"/Scenarios", "/Fmsy")
model_run <- SS_output(stf_dir)

iniYr<-1982
FinYr<-2028
year.range<-iniYr:FinYr 
########
# F values
#***********
fbar_raw<-model_run$derived_quants[substr(model_run$derived_quants$Label,1,2) == "F_",]
fbar_lo <- fbar_raw$Value - qnorm(1-(1-0.95)/2)*fbar_raw$StdDev 
fbar_hi <- fbar_raw$Value + qnorm(1-(1-0.95)/2)*fbar_raw$StdDev 
fbar <- cbind.data.frame(year=substr(fbar_raw$Label,3,6),Value=fbar_raw$Value, fbar_lo, fbar_hi)
fbar$datatype<-ifelse(fbar$year%in%tail(year.range,3),"forecast","assessment")

ggplot(fbar[fbar$year>2009,],aes(year,Value,group=1,ymin=fbar_lo,ymax=fbar_hi,col=datatype))+geom_point()+geom_line()+labs(y="Fbar (1-10)",x="")+geom_linerange(aes(x=year))+geom_vline(xintercept="2026",lty=3,col="blue")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "")
dev.print(jpeg,paste0(stf_out,"/stf_fbar.jpg"), width = 7, height = 5, res = 300, units = "in")


# SSB values
#************
ssb_raw <- rbind(model_run$derived_quants[substr(model_run$derived_quants$Label,1,5)	=="SSB_1",],
                 model_run$derived_quants[substr(model_run$derived_quants$Label,1,5)	=="SSB_2",])

#ssb <- cbind.data.frame(Value=ssb_raw$Value, lo=ssb_raw$Value-1.96*ssb_raw$StdDev, hi=ssb_raw$Value+1.96*ssb_raw$StdDev)
#plot(ssb$Value, type="l")
#lines(ssb$lo)
#lines(ssb$hi)

ssb_lo <- ssb_raw$Value - qnorm(1-(1-0.95)/2)*ssb_raw$StdDev   
ssb_hi <- ssb_raw$Value + qnorm(1-(1-0.95)/2)*ssb_raw$StdDev
ssb <- cbind.data.frame(year=substr(ssb_raw$Label,5,8),Value=ssb_raw$Value, ssb_lo, ssb_hi)
 
ssb$datatype<-ifelse(ssb$year%in%tail(year.range,3),"forecast","assessment")

ggplot(ssb[ssb$year>2009,],aes(year,Value,group=1,ymin=ssb_lo,ymax=ssb_hi,col=datatype))+geom_point()+geom_line()+labs(y="SSB (t)",x="")+geom_linerange(aes(x=year))+geom_vline(xintercept="2026",lty=3,col="blue")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "")
dev.print(jpeg,paste0(stf_out,"/stf_ssb.jpg"), width = 7, height = 5, res = 300, units = "in")





rec_raw <- rbind(model_run$derived_quants[substr(model_run$derived_quants$Label,1,6)	=="Recr_1",],
                 model_run$derived_quants[substr(model_run$derived_quants$Label,1,6)	=="Recr_2",])

## recruitment is assumed to have a lognormal distribution:
rec_hi <- exp(log(rec_raw$Value) + qnorm(1-(1-0.95)/2) * sqrt(log(1 + (rec_raw$StdDev/rec_raw$Value) * (rec_raw$StdDev/rec_raw$Value))))  
rec_lo <- exp(log(rec_raw$Value) - qnorm(1-(1-0.95)/2) * sqrt(log(1 + (rec_raw$StdDev/rec_raw$Value) * (rec_raw$StdDev/rec_raw$Value)))) 

rec <- cbind.data.frame(year=substr(rec_raw$Label,6,9),Value=rec_raw$Value, rec_lo, rec_hi)
rec$datatype<-ifelse(rec$year%in%tail(year.range,3),"forecast","assessment")
#remove CI for the forecast period (it is based in geometric mean):
rec$rec_lo<-ifelse(rec$datatype=="assessment",rec$rec_lo,NA)
rec$rec_hi<-ifelse(rec$datatype=="assessment",rec$rec_hi,NA)

ggplot(rec[rec$year>2009,],aes(year,Value,group=1,ymin=as.numeric(rec_lo),ymax=as.numeric(rec_hi),col=datatype))+geom_point()+geom_line()+labs(y="Recruitment (millions)",x="")+geom_linerange(aes(x=year))+geom_vline(xintercept="2026",lty=3,col="blue")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(legend.position = "")
dev.print(jpeg,paste0(stf_out,"/stf_rec.jpg"), width = 7, height = 5, res = 300, units = "in")


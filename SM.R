SOLANUM <- function(weather,
  crop_parameters, 
  dates, 
  plant_density = 3.7,
  env_parameters = list(TN=4, TO=17, TX=28, Pc=12, w=0.5),
  potential_growth_module = TRUE,
  lateblight_limited_module = TRUE,
  water_limited_module = FALSE,
  frost_limited_module = FALSE,
  severity_data,
  umbral=0.5) {

  if(!(potential_growth_module)){stop("Activate potential growth module")}
  
  # Checking weather-data names
  if(!("date" %in% colnames(weather))){stop("There is no 'date' data")}
  if(!("tmin" %in% colnames(weather))){stop("There is no 'tmin' data")}
  if(!("tmax" %in% colnames(weather))){stop("There is no 'tmax' data")}
  if(!("srad" %in% colnames(weather))){stop("There is no 'srad' data")}

#  if(water_limited_module){
#    
#    stop("Activate potential growth module")}
  
  
##################        INICIALIZACION DE VARIABLES           ############################

sowing = as.Date(dates$sowing)
harvest = as.Date(dates$harvest)
EDay = as.Date(dates$emergence) - sowing + 1  # Emergence day
EDay = as.numeric(EDay)

Tb_0 = env_parameters$TN  # Minimum Temperature for tuber initiation
To_0 = env_parameters$TO  # Optimum temperature for tuber initiation
Tu_1 = env_parameters$TX  # Maximum Temperature for tuber initiation

Pc = env_parameters$Pc    # Critical Photperiod
 w = env_parameters$w     # Photperiod Sensivity
 
 A_0 = crop_parameters$A
Tu_0 = crop_parameters$tu
   b = crop_parameters$b
   
wmax = crop_parameters$wmax
tm = crop_parameters$tm 
te = crop_parameters$te

plantDensity = plant_density
RUE = crop_parameters$RUE
DMcont = crop_parameters$DMc

severity_data$date = severity_data$Date 
weather = full_join(weather, severity_data,"date")
 
v=0.0  # variablity
TDM = 0.0
TDM_LB = 0.0
cTII = 0.0
Tac = 0.0
tt<-0.0
PDay<--1
ttfixed<--1
a = (log(2))/log((Tu_1-Tb_0)/(To_0-Tb_0))
#############################################################################################
#                                                                                           #
#                                 Inicio de la simulacion                                   #
#                                                                                           #
#############################################################################################
rnd<-runif(1)
rnd<-2*rnd-1
rdm = (log((1+rnd)/(1-rnd)))/1.82 ## log: es logaritmo natural
# levantamos toda la base de datos de clima en un data frame
climate <- weather   # abrir un archivo ascii
climate$date<-as.Date(climate$date)

# seleccionamos datos de clima correspondiente al periodo de simulacion
climate <- climate[climate$date >= sowing & climate$date <= harvest, ]
# retiro las columnas Date, Tmin y Tmax del dataframe y los guardo en vectores por separado
DATE<-as.Date(climate$date)
TMIN<-as.numeric(climate$tmin)
TMAX<-as.numeric(climate$tmax)
TT<-thermalTime(DATE,TMIN,TMAX,sowing,harvest,EDay,c(0, 12, 24, 35)) # se calcula el TT
climate<-data.frame(climate,TT=TT$tt) # agrego la columna TT calculada al df climate
timeDuration<-as.numeric(as.Date(harvest)-as.Date(sowing) +1)

# OUTPUTS
tdm <- vector(mode="numeric", length=timeDuration)
dty <- vector(mode="numeric", length=timeDuration)
fty <- vector(mode="numeric", length=timeDuration)
cc <- vector(mode="numeric", length=timeDuration)


tdm_LB <- vector(mode="numeric", length=timeDuration)
dty_LB <- vector(mode="numeric", length=timeDuration)
fty_LB <- vector(mode="numeric", length=timeDuration)
cc_LB <- vector(mode="numeric", length=timeDuration)


# OUTPUTS
Pindex<-vector()
Tindex<-vector()
TII<-vector()
Part<-vector()
ctii<-vector()


for(i in 1:timeDuration)
{
  N = climate$sunsh[i]
  Tmax = climate$tmax[i]
  Tmin = climate$tmin[i]
  Tav = (Tmax+Tmin)/2
  tt = climate$TT[i];
  Tindex[i] = ifelse(Tav<Tb_0,0,ifelse(Tav>Tu_1,0,(2*((Tav-Tb_0)^a)*((To_0-Tb_0)^a)-((Tav-Tb_0)^(2*a)))/((To_0-Tb_0)^(2*a))))
  Pindex[i] = ifelse(N>Pc,exp(-1*(w*(N-Pc))),1)
  TII[i] = Tindex[i]*Pindex[i]
  cTII = cTII+TII[i];
  ctii[i]<-cTII
  if(cTII>20 && ttfixed==-1) ttfixed<-tt
  Tac<-tt
  Part[i] = A_0*exp(-1*(exp((-1*(Tac-Tu_0))/b)))
  #print(cTII)
}
if(cTII>=20){
  dato<-ttfixed
}else{
  dato<-climate$TT[timeDuration];
}
Tu_cTII=(dato+b)/Tu_0

# OUTPUTS
Part_cTII<-vector()
HI_cTII<-vector()
dW<-vector()
dty1<-vector()
dty2<-vector()

for(i in 1:timeDuration)
{
  tt = climate$TT[i];
  Tac<-tt
  Part_cTII[i] = A_0*exp(-1*(exp((-1*(Tac-Tu_0*Tu_cTII))/b)))
  HI_cTII[i] = ifelse(Tu_cTII<=1,Part[i],Part_cTII[i])
  
  canopy = wmax*exp(-1*(tm/(Tac*plantDensity)))*(1+(te-Tac)/(te-tm))*((Tac/te)^(te/(te-tm)))
  canopy1 = rdm*v*canopy+canopy
  PDay=PDay+1
  DAE = ifelse(PDay>=EDay,PDay-EDay,0)
  CC = ifelse(DAE<=0,0,ifelse(canopy1>0,canopy1,10^(-6)))  
  PAR = climate$srad[i]*0.5
  dW[i] = (RUE*CC*PAR)/100
  TDM = TDM+dW[i];
  
  DTY1 = TDM*Part[i]
  DTY2 = TDM*HI_cTII[i]
  DTY = ifelse(Tu_cTII<=1,DTY1,DTY2)
  FTY = DTY/DMcont
  
  tdm[i]<-TDM
  dty[i]=DTY;
  fty[i]=FTY;
  cc[i]=CC;  
  dty1[i]=DTY1;
  dty2[i]=DTY2;
  
###################################################################
  
  CC_LB = ifelse(is.na(weather$SimSeverity[i]), CC, CC - umbral*(weather$SimSeverity[i]/100))
  dW_LB = (RUE*CC_LB*PAR)/100
  TDM_LB = TDM_LB+dW_LB;
  
  DTY1_LB = TDM_LB*Part[i]
  DTY2_LB = TDM_LB*HI_cTII[i]
  DTY_LB = ifelse(Tu_cTII<=1,DTY1_LB,DTY2_LB)
  FTY_LB = DTY_LB/DMcont
  
  tdm_LB[i]<-TDM_LB
  dty_LB[i]=DTY_LB
  fty_LB[i]=FTY_LB
  cc_LB[i]=CC_LB;  


}
return(out=data.frame(Fecha=climate$date,TT=climate$TT,cc,Tmin=climate$tmin,Tmax=climate$tmax,N=climate$sunsh,sr=climate$srad,dW,Tindex,Pindex,TII,ctii,Part,Part_cTII,HI_cTII,
                      dty1,dty2,dty,fty,tdm,
                      dty_LB,fty_LB,tdm_LB,cc_LB))
#write.csv(df,"out_cameroon.csv")
#write.csv(df,"out_congo.csv")


}




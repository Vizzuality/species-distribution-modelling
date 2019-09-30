require(vspt)




library(raster)

#climatic data
bio4=raster("/Users/gretacvega/Documents/Copernicus/forest/biovar_Sweden/biovar_historical/biovar04_era5_Sweden_1986-2005.tif")
bio6=raster("/Users/gretacvega/Documents/Copernicus/forest/biovar_Sweden/biovar_historical/biovar06_era5_Sweden_1986-2005.tif")
bio16=raster("/Users/gretacvega/Documents/Copernicus/forest/biovar_Sweden/biovar_historical/biovar16_era5_Sweden_1986-2005.tif")


bios=stack(bio4,bio6,bio16)
names(bios)=c("bio4","bio6","bio16")

library(rgdal)
swePoly=readOGR("/Users/gretacvega/Documents/Copernicus/forest/swedenBasemap/", "SWE_adm0")

bios_sweStack=stack(mask(bios,swePoly))
names(bios_sweStack)



library(rgbif)

####queing downloads

keys=numeric()
for(sp in 1:length(sweden_species)){
  key=name_suggest(q=sweden_species[sp], rank='species')[,1]
  keys=rbind(keys, key)
}


queries <- list()
for (i in 2:length(keys$key)) {
  queries[[i]] <- occ_download_prep(
    paste0("taxonKey = ", keys$key[i]),'country = SE', 'hasCoordinate = TRUE', user = "greta89", pwd = "h82YGcUEaqakSBG", email = "gretacvega@gmail.com"

  )
}
out <- occ_download_queue(.list = queries)
out

occ_download_list("greta89", pwd = "h82YGcUEaqakSBG" )$results[,c('key','created', 'status','size','totalRecords')]

dat <- occ_download_get("0004350-190813142620410")
df=occ_download_import(dat)
summary(df)
checkDwld=occ_download_list("username", pwd = "****")


dat <- occ_download_get("0004335-190813142620410")
df=occ_download_import(dat)
summary(df)
plot(df$decimalLatitude~df$decimalLongitude)

spData=subset(df,select=c("genus", "specificEpithet", "decimalLatitude", "decimalLongitude") )
head(spData)
summary(spData)
spData$sciname=paste(spData$genus, spData$specificEpithet)

presence.absence.raster <- function (mask.raster,species.data,raster.label="") {
  require(raster)

  # set the background cells in the raster to 0
  mask.raster[!is.na(mask.raster)] <- 0

  #set the cells that contain points to 1
  speciesRaster <- rasterize(species.data,mask.raster,field=1)
  speciesRaster <- merge(speciesRaster,mask.raster)

  #label the raster
  names(speciesRaster) <- raster.label
  return(speciesRaster)
}

pa.rasterWorld<-presence.absence.raster(mask.raster=raster(bios_sweStack,1), species.data=spData[,c("decimalLongitude","decimalLatitude")], raster.label=spData$sciname[1])

plot(pa.rasterWorld)

species<-as.data.frame(rasterToPoints(pa.rasterWorld))
head(species)
summary(species)

presences <- which(species[,3]==1)
resp.var <- as.numeric(species[presences,3]) # species presences vector (only 'ones')
resp.xy <- species[presences,1:2] # coordinates of presence cells

library(biomod2)
setwd("~/Documents/Copernicus/forest/biomod_output")

# prepare presence absence data
biomodDataPOs <- BIOMOD_FormatingData(
  resp.var=resp.var,
  expl.var=bios_sweStack,
  resp.xy=resp.xy,
  resp.name=spData$sciname[1],
  PA.nb.rep=2,
  PA.nb.absences=300,
  PA.strategy="sre",
  PA.sre.quant=0.25)

# set modelling options
ModOptions <- BIOMOD_ModelingOptions() #saving object with default options


#calibrate the models

ModelOut <- BIOMOD_Modeling(biomodDataPOs,
                            models = c('CTA', 'RF','GLM'),
                            models.options = ModOptions,
                            DataSplit = 80,
                            NbRunEval = 2,
                            Yweights = NULL,
                            VarImport = 3,
                            models.eval.meth = c('TSS','ROC'),
                            SaveObj = T,
                            rescal.all.models = T,
                            do.full.models = T,
                            modeling.id = paste(spData$sciname[1], biomodDataPOs@PA.strategy, as.character(format(Sys.time(), '%y%m%d')), sep="_"))

#Total number of model runs : 18 = 3 algo * 2 PA * (2 rep + 1 fullrep)


#project the models to current climate

new.env <- bios_sweStack
proj.name='baseline_current'

ProjBas <- BIOMOD_Projection(
  modeling.output = ModelOut,
  selected.models = "all",
  new.env = new.env ,
  xy.new.env = NULL,
  proj.name = proj.name,
  binary.meth = 'ROC',
  filtered.meth = NULL,
  build.clamping.mask = T,
  compress = T)

#ensemble modelling

EnsMod <- BIOMOD_EnsembleModeling(
  modeling.output = ModelOut,
  chosen.models = "all",
  em.by = "all", # 5 available options
  VarImport = 1,
  eval.metric = "ROC",
  eval.metric.quality.threshold = NULL, # here we can set up a threshold to choose the models to use for the ensemble
  models.eval.meth = c("TSS","ROC"),
  prob.mean = T,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = F,
  prob.mean.weight.decay = 'proportional',
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05)


#ensemble forecasting

EnsBas <- BIOMOD_EnsembleForecasting(
  EM.output =  EnsMod,
  projection.output = ProjBas,
  binary.meth = "ROC",
  filtered.meth = NULL)

str(EnsBas)

filelist <- list.files(paste0(gsub(" ",".",spData$sciname[1]),"/proj_",proj.name) )

enslist <- filelist[grep("ensemble", filelist)]  # selecting the files for the consensus projections (those that have the word "ensemble" in their name; the remaining ones in the same folder are the single.model projections or the Clamping mask )

probEns=enslist[grep(paste0("proj_", proj.name, "_", gsub(" ", ".", spData$sciname[1]), "_ensemble.grd"), enslist)]

ensp <- stack(paste(gsub(" ", ".", spData$sciname[1]),paste0("proj_",proj.name), probEns, sep="/"))

#ensp has the probabilistic projections, from here we want the ca and the cv
names(ensp)
EM.methods=c("EMca", "EMcv","EMmean", "EMmedian")

ca_raster=raster(ensp, grep(EM.methods[1], names(ensp)))
cv_raster=raster(ensp, grep(EM.methods[2], names(ensp)))

par(mfrow=c(2,3))
plot(ensp)



binEns=enslist[grep(paste0("proj_", proj.name, "_", gsub(" ", ".", spData$sciname[1]), "_ensemble_ROCbin.grd"), enslist)]

ensb <- stack(paste(gsub(" ", ".", spData$sciname[1]),paste0("proj_",proj.name), binEns, sep="/"))

#ensb has the binary projections, from here we want the mean and the median
names(ensb)
mean_bin_raster=raster(ensb, grep(EM.methods[3], names(ensb)))
median_bin_raster=raster(ensb, grep(EM.methods[4], names(ensb)))

#par(mfrow=c(2,3))
#plot(ensb)
#


outputCurrent=stack(ca_raster, cv_raster, mean_bin_raster, median_bin_raster )

plot(outputCurrent)
names(outputCurrent)=EM.methods


dfEns=as.data.frame(rasterToPoints(outputCurrent))

currentProp=sum(dfEns$EMca>=500)/length(dfEns$EMca)

head(dfEns)


library(ggplot2)
sweFort=fortify(swePoly)
head(sweFort)

pa_df=as.data.frame(rasterToPoints(pa.rasterWorld))
names(pa_df)[3]="presences"
presencesDF=pa_df[pa_df$presences==1,]

ggplot()+
  #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="#276419", colour="transparent")+
  geom_raster(data=dfEns, aes(x=x, y=y,alpha =as.factor(EMmean), fill=EMcv), hjust=0.5, vjust=0.5,
              interpolate=FALSE)+
  scale_alpha_manual(values=c(0.75,1))+
  #scale_fill_gradientn(colours=c("black","white","red"),breaks=c(0,500,1000), limits=c(0,1000))+
  #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="transparent", colour="black")+
  geom_point(data=presencesDF, aes(x=x,y=y), colour="red")+
  coord_fixed(ratio=1.3)+
  theme_bw()

ggplot()+
  #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="#276419", colour="transparent")+
  geom_raster(data=dfEns, aes(x=x, y=y, fill=EMca), hjust=0.5, vjust=0.5,
              interpolate=FALSE)+
  #scale_alpha_manual(values=c(0.75,1))+
  scale_fill_gradientn(colours=c("green","white","red"),breaks=c(0,500,1000), limits=c(0,1000))+
  #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="transparent", colour="black")+
  #geom_point(data=presencesDF, aes(x=x,y=y))+
  coord_fixed(ratio=1.3)+
  theme_bw()+
  ggtitle(paste(spData$sciname[1],"Current modeled distribution" ))

###Projecting to the future
#4 6 16

decades=seq(2020,2090,10)

#for each decade
d=2

futPath="/Users/gretacvega/Documents/Copernicus/forest/biovar_Sweden/biovar_future"

propSuitable_rcp45=numeric()
propSuitable_rcp85=numeric()

dfs_rcp45=list()
dfs_rcp85=list()

for (d in 3:length(decades)){
  scenario45_bios=list.files(path=paste(futPath, decades[d], sep="/"), pattern ="rcp45")

  scenario85_bios=list.files(path=paste(futPath, decades[d], sep="/"), pattern ="rcp85")

  #here it would be good to check that both scenarios have the same gcm models

  bio4_rcp45=scenario45_bios[grep("biovar04", scenario45_bios)]
  bio6_rcp45=scenario45_bios[grep("biovar06", scenario45_bios)]
  bio16_rcp45=scenario45_bios[grep("biovar16", scenario45_bios)]

  bio4_rcp45_stack=stack(paste(futPath, decades[d], bio4_rcp45,sep="/"))
  bio4_rcp45_stack_swe=stack(mask(bio4_rcp45_stack,swePoly))
  bio4_rcp45_stack_swe_mean=mean(bio4_rcp45_stack_swe)

  bio6_rcp45_stack=stack(paste(futPath, decades[d], bio6_rcp45,sep="/"))
  bio6_rcp45_stack_swe=stack(mask(bio6_rcp45_stack,swePoly))
  bio6_rcp45_stack_swe_mean=mean(bio6_rcp45_stack_swe)

  bio16_rcp45_stack=stack(paste(futPath, decades[d], bio16_rcp45,sep="/"))
  bio16_rcp45_stack_swe=stack(mask(bio16_rcp45_stack,swePoly))
  bio16_rcp45_stack_swe_mean=mean(bio16_rcp45_stack_swe)



  bio4_rcp85=scenario85_bios[grep("biovar04", scenario85_bios)]
  bio6_rcp85=scenario85_bios[grep("biovar06", scenario85_bios)]
  bio16_rcp85=scenario85_bios[grep("biovar16", scenario85_bios)]

  bio4_rcp85_stack=stack(paste(futPath, decades[d], bio4_rcp85,sep="/"))
  bio4_rcp85_stack_swe=stack(mask(bio4_rcp85_stack,swePoly))
  bio4_rcp85_stack_swe_mean=mean(bio4_rcp85_stack_swe)

  bio6_rcp85_stack=stack(paste(futPath, decades[d], bio6_rcp85,sep="/"))
  bio6_rcp85_stack_swe=stack(mask(bio6_rcp85_stack,swePoly))
  bio6_rcp85_stack_swe_mean=mean(bio6_rcp85_stack_swe)

  bio16_rcp85_stack=stack(paste(futPath, decades[d], bio16_rcp85,sep="/"))
  bio16_rcp85_stack_swe=stack(mask(bio16_rcp85_stack,swePoly))
  bio16_rcp85_stack_swe_mean=mean(bio16_rcp85_stack_swe)

  bios_rcp45_mean=stack(bio4_rcp45_stack_swe_mean, bio6_rcp45_stack_swe_mean, bio16_rcp45_stack_swe_mean)
  bios_rcp85_mean=stack(bio4_rcp85_stack_swe_mean, bio6_rcp85_stack_swe_mean, bio16_rcp85_stack_swe_mean)


  names(bios_rcp45_mean)=names(bios)
  names(bios_rcp85_mean)=names(bios)

  proj.name.rcp45=paste0("rcp45_", decades[d])
  proj.name.rcp85=paste0("rcp85_", decades[d])

  ProjFut_rcp45 <- BIOMOD_Projection(
    modeling.output = ModelOut,
    selected.models = 'all',
    new.env = bios_rcp45_mean,
    xy.new.env = NULL,
    proj.name = proj.name.rcp45,
    binary.meth = "ROC",
    filtered.meth = NULL,
    build.clamping.mask = T,
    compress = T)




  ProjFut_rcp85 <- BIOMOD_Projection(
    modeling.output = ModelOut,
    selected.models = 'all',
    new.env = bios_rcp85_mean,
    xy.new.env = NULL,
    proj.name = proj.name.rcp85,
    binary.meth = "ROC",
    filtered.meth = NULL,
    build.clamping.mask = T,
    compress = T)


  EnsFut_rcp45 <- BIOMOD_EnsembleForecasting(
    projection.output = ProjFut_rcp45,
    EM.output = EnsMod,
    binary.meth = "ROC")


  EnsFut_rcp85 <- BIOMOD_EnsembleForecasting(
    projection.output = ProjFut_rcp85,
    EM.output = EnsMod,
    binary.meth = "ROC")


  filelist_rcp45 <- list.files(paste0(gsub(" ",".",spData$sciname[1]),"/proj_",proj.name.rcp45) )
  enslist_rcp45 <- filelist_rcp45[grep("ensemble", filelist_rcp45)]  # selecting the files for the consensus projections (those that have the word "ensemble" in their name; the remaining ones in the same folder are the single.model projections or the Clamping mask )
  probEns_rcp45=enslist_rcp45[grep(paste0("proj_", proj.name.rcp45, "_", gsub(" ", ".", spData$sciname[1]), "_ensemble.grd"), enslist_rcp45)]
  ensp_rcp45 <- stack(paste(gsub(" ", ".", spData$sciname[1]),paste0("proj_",proj.name.rcp45), probEns_rcp45, sep="/"))

  ca_rcp45_raster=raster(ensp_rcp45, grep(EM.methods[1], names(ensp_rcp45)))
  cv_rcp45_raster=raster(ensp_rcp45, grep(EM.methods[2], names(ensp_rcp45)))

  binEns_rcp45=enslist_rcp45[grep(paste0("proj_", proj.name.rcp45, "_", gsub(" ", ".", spData$sciname[1]), "_ensemble_ROCbin.grd"), enslist_rcp45)]
  ensb_rcp45 <- stack(paste(gsub(" ", ".", spData$sciname[1]),paste0("proj_",proj.name.rcp45), binEns_rcp45, sep="/"))
  mean_bin_rcp45_raster=raster(ensb_rcp45, grep(EM.methods[3], names(ensb_rcp45)))
  median_bin_rcp45_raster=raster(ensb_rcp45, grep(EM.methods[4], names(ensb_rcp45)))


  output_rcp45=stack(ca_rcp45_raster, cv_rcp45_raster, mean_bin_rcp45_raster, median_bin_rcp45_raster )

  names(output_rcp45)=EM.methods


  dfEns_rcp45=as.data.frame(rasterToPoints(output_rcp45))

  dfs_rcp45[[d]]=dfEns_rcp45

  propSuitable_rcp45[d]=sum(dfEns_rcp45$EMca>=500)/length(dfEns_rcp45$EMca)

  # ggplot()+
  #   #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="#276419", colour="transparent")+
  #   geom_raster(data=dfEns_rcp45, aes(x=x, y=y,alpha =as.factor(EMmean), fill=EMcv), hjust=0.5, vjust=0.5,
  #               interpolate=FALSE)+
  #   scale_alpha_manual(values=c(0.75,1))+
  #   #scale_fill_gradientn(colours=c("black","white","red"),breaks=c(0,500,1000), limits=c(0,1000))+
  #   #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="transparent", colour="black")+
  #   #geom_point(data=presencesDF, aes(x=x,y=y), colour="red")+
  #   coord_fixed(ratio=1.3)+
  #   theme_bw()+
  #   ggtitle(paste("rcp45", decades[d]))

  gg45=ggplot()+
    #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="#276419", colour="transparent")+
    geom_raster(data=dfEns_rcp45, aes(x=x, y=y, fill=EMca), hjust=0.5, vjust=0.5,
                interpolate=FALSE)+
    #scale_alpha_manual(values=c(0.75,1))+
    scale_fill_gradientn(colours=c("green","white","red"),breaks=c(0,500,1000), limits=c(0,1000))+
    #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="transparent", colour="black")+
    #geom_point(data=presencesDF, aes(x=x,y=y))+
    coord_fixed(ratio=1.3)+
    theme_bw()+
    ggtitle(paste(spData$sciname[1], "scenario: rcp45, decade:", decades[d]))

  ggsave(plot=gg45, filename=paste0(gsub(" ",".",spData$sciname[1]),"/plots/",gsub(" ","_",spData$sciname[1]), "_rcp45_", decades[d], ".png"), width = 15, height = 15, units = "cm", dpi = 300 )


  filelist_rcp85 <- list.files(paste0(gsub(" ",".",spData$sciname[1]),"/proj_",proj.name.rcp85) )
  enslist_rcp85 <- filelist_rcp85[grep("ensemble", filelist_rcp85)]  # selecting the files for the consensus projections (those that have the word "ensemble" in their name; the remaining ones in the same folder are the single.model projections or the Clamping mask )
  probEns_rcp85=enslist_rcp85[grep(paste0("proj_", proj.name.rcp85, "_", gsub(" ", ".", spData$sciname[1]), "_ensemble.grd"), enslist_rcp85)]
  ensp_rcp85 <- stack(paste(gsub(" ", ".", spData$sciname[1]),paste0("proj_",proj.name.rcp85), probEns_rcp85, sep="/"))

  ca_rcp85_raster=raster(ensp_rcp85, grep(EM.methods[1], names(ensp_rcp85)))
  cv_rcp85_raster=raster(ensp_rcp85, grep(EM.methods[2], names(ensp_rcp85)))

  binEns_rcp85=enslist_rcp85[grep(paste0("proj_", proj.name.rcp85, "_", gsub(" ", ".", spData$sciname[1]), "_ensemble_ROCbin.grd"), enslist_rcp85)]
  ensb_rcp85 <- stack(paste(gsub(" ", ".", spData$sciname[1]),paste0("proj_",proj.name.rcp85), binEns_rcp85, sep="/"))
  mean_bin_rcp85_raster=raster(ensb_rcp85, grep(EM.methods[3], names(ensb_rcp85)))
  median_bin_rcp85_raster=raster(ensb_rcp85, grep(EM.methods[4], names(ensb_rcp85)))


  output_rcp85=stack(ca_rcp85_raster, cv_rcp85_raster, mean_bin_rcp85_raster, median_bin_rcp85_raster )

  names(output_rcp85)=EM.methods


  dfEns_rcp85=as.data.frame(rasterToPoints(output_rcp85))
  dfs_rcp85[[d]]=dfEns_rcp85
  propSuitable_rcp85[d]=sum(dfEns_rcp85$EMca>=500)/length(dfEns_rcp85$EMca)


  gg85=ggplot()+
    #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="#276419", colour="transparent")+
    geom_raster(data=dfEns_rcp85, aes(x=x, y=y, fill=EMca), hjust=0.5, vjust=0.5,
                interpolate=FALSE)+
    #scale_alpha_manual(values=c(0.75,1))+
    scale_fill_gradientn(colours=c("green","white","red"),breaks=c(0,500,1000), limits=c(0,1000))+
    #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="transparent", colour="black")+
    #geom_point(data=presencesDF, aes(x=x,y=y))+
    coord_fixed(ratio=1.3)+
    theme_bw()+
    ggtitle(paste(spData$sciname[1], "scenario: rcp85, decade:", decades[d]))

  ggsave(plot=gg85, filename=paste0(gsub(" ",".",spData$sciname[1]),"/plots/",gsub(" ","_",spData$sciname[1]), "_rcp85_", decades[d], ".png"), width = 15, height = 15, units = "cm", dpi = 300 )

}


df_props=data.frame(decades, propSuitable_rcp45, propSuitable_rcp85)
df_props=rbind(df_props, c(2010, currentProp, currentProp))
library(reshape2)
df_props_long=melt(df_props, id.vars = "decades")

ggplot(data=df_props_long, aes(x=decades, y=value, colour=variable))+
  geom_point()+
  geom_line()+
  scale_y_continuous(name="proportion of suitable area", limits = c(0,1))+
  scale_x_continuous(name = "time", breaks = c(2010, decades), labels = c("current", decades))+
  theme_bw()


##

# library(hexbin)
# # ggplot(raster_df, aes(x, y, fill=blabla)) + geom_hex(stat="identity")
# #
# hx=hexbin(x=dfEns_rcp85$x,y=dfEns_rcp85$y,IDs = TRUE)
# plot(hx)
# xx=dfEns_rcp85
# coordinates(xx) <- xx[,c('x','y')]
#
# hex_points <- spsample(xx, type = "hexagonal", cellsize = 0.25)
# hex_grid <- HexPoints2SpatialPolygons(hex_points, dx = 0.25)
# plot(hex_grid)
#
# library(spatialEco)
# pts.poly <- point.in.poly(xx, hex_grid)


# ggplot()+
#   #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="#276419", colour="transparent")+
#   geom_hex(data=dfEns_rcp85, aes(x=x, y=y, fill=EMca),stat="identity")+
#   #scale_alpha_manual(values=c(0.75,1))+
#   scale_fill_gradientn(colours=c("green","white","red"),breaks=c(0,500,1000), limits=c(0,1000))+
#   #geom_polygon(data=sweFort,aes(group = group, x=long, y=lat),fill="transparent", colour="black")+
#   #geom_point(data=presencesDF, aes(x=x,y=y))+
#   coord_fixed(ratio=1.3)+
#   theme_bw()+
#   ggtitle(paste(spData$sciname[1], "scenario: rcp85, decade:", decades[d]))

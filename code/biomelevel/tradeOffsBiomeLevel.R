# get pareto front for the biome level

require(emoa)
require(foreign)
require(geometry)
require(lattice)
require(sp)
require(maptools)
require(hexbin)
require(raster)
require(sf)

library(tidyverse)
getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

ws <- getCurrentFileLocation()
setwd(ws)

outFolder <- "../../results/biomes/"
dataFolder <- "../../data/"

if(!dir.exists(outFolder))
  dir.create(outFolder, recursive = TRUE)

########################
# read input data
#########################

###
# plant/crop functional types
###
# contains pasture now
pfts <- read.table(paste0(dataFolder, "pft_cft.txt"), header=TRUE, sep=";")

# new constraint
constraints <- read.table(paste0(dataFolder, "constraints.csv"), header=TRUE)

###
# biome map
###
biome <- read.table(paste0(dataFolder, "biome.txt"), header=TRUE)
biome.legend <- read.table(paste0(dataFolder,"legend_biomes_port_AT.txt"), header=TRUE)

# get area for biomes

getBiomeAreasFromDataFrame <- function(theDf)
{
  theDf.sp <- theDf
  # make spatial
  coordinates(theDf.sp) <- ~ Lon + Lat

  gridded(theDf.sp) <- TRUE

  # for llgridlines
  proj4string(theDf.sp) <- "+proj=longlat +ellps=WGS84"
  newproj <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  #solutions.sp.prj <- spTransform(solutions.sp, newproj)

  theDf.spGrid <- as(theDf.sp, "SpatialGridDataFrame")

  #head(food.spGrid@data)
  theDf.raster <- brick(theDf.spGrid)
  crs(theDf.raster) <- "+proj=longlat +ellps=WGS84"
  theDf.poly <- rasterToPolygons(theDf.raster, dissolve = TRUE)
  #spplot(theDf.poly, zcol="ATbiome")
  #project polygons
  theDf.polyMollweide <- spTransform(theDf.poly, newproj)
  spplot(theDf.polyMollweide, zcol="ATbiome")
  # convert to sf for simpler area calculation
  SfpolyMollweide <- st_as_sf(theDf.polyMollweide)
  # calculate area
  SfpolyMollweide$area_qm <- st_area(SfpolyMollweide)
  SfpolyMollweide$area_qkm <- SfpolyMollweide$area_qm
  units(SfpolyMollweide$area_qkm) <- with(units::ud_units, km^2)

  return(st_drop_geometry(SfpolyMollweide[c("ATbiome", "area_qkm")]))
}

biomeArea <- getBiomeAreasFromDataFrame(biome)

biomeArea <- merge(biome.legend, biomeArea, by.x="biomID", by.y="ATbiome")
biomeArea$area_qkm <- as.numeric(biomeArea$area_qkm)
biome <- merge(biome, biomeArea, by.x="ATbiome", by.y="biomID")

biome$LONLAT <- paste(sep="/", as.character(biome$Lon), as.character(biome$Lat))

#***************************************************
# handling of protected areas --------------------
#***************************************************
protectedAreas <- read.table(paste0(dataFolder,"gridlist_hurtt_RNDM_1deg_protarea.txt"), header=FALSE)
names(protectedAreas) <- c("Lon", "Lat", "Protected")

# correct for the shift of .5°

protectedAreas$Lon <- protectedAreas$Lon +.5
protectedAreas$Lat <- protectedAreas$Lat +.5

protectedAreas$LONLAT <- paste(sep="/", as.character(protectedAreas$Lon), as.character(protectedAreas$Lat))


###
# objective values
###

# store them as a list inside a list to facilitate a loop over the two rcps and the two time periods
timePeriods <- c("2033-2042", "2090-2099")
rcps <- c("26", "60")

for (aRcp in rcps)
{
  for (aTimePeriod in timePeriods)
  {
    #get corresponding constraints
    idx_constraints <- which(constraints$rcp == aRcp & constraints$timePeriod == aTimePeriod)
    theConstraints <- constraints[idx_constraints,]

    baseFn <- paste0(dataFolder, "/ESquant/ESquant_rcp",aRcp, "_", aTimePeriod, "_INT_9412_v8b_4gcms_")
    food <- read.table(paste0(baseFn,  "FoPro.txt"), header=TRUE)
    #summary(food)
    cstor <- read.table(paste0(baseFn,"CStor.txt"), header=TRUE)
    #summary(cstor)
    water <- read.table(paste0(baseFn,"WaSu.txt"), header=TRUE)
    #summary(water)
    forage <- read.table(paste0(baseFn,"Forage.txt"), header=TRUE)
    #summary(forage)

    cstor$LONLAT <- paste(sep="/", as.character(cstor$Lon), as.character(cstor$Lat))
    water$LONLAT <- paste(sep="/", as.character(water$Lon), as.character(water$Lat))
    food$LONLAT <- paste(sep="/", as.character(food$Lon), as.character(food$Lat))
    forage$LONLAT <- paste(sep="/", as.character(forage$Lon), as.character(forage$Lat))

        # reorder biomes
    idx <- match(cstor$LONLAT, biome$LONLAT)
    head(cbind(cstor$LONLAT, forage$LONLAT, water$LONLAT, food$LONLAT, biome$LONLAT[idx]))
    biome <- biome[idx,]
    
    # run simple tests
    test1 <- which(!(biome$LONLAT %in% cstor$LONLAT))
    test2 <- which(!(water$LONLAT %in% cstor$LONLAT))
    test3 <- which(!(food$LONLAT %in% cstor$LONLAT))
    test4 <- which(!(forage$LONLAT %in% cstor$LONLAT))


    if(length(c(test1, test2, test3, test4)) > 0)
    {
      print(paste("Problem matching objective data frames for rcp", aRcp, " and time period ", aTimePeriod))
      print("Skipping pareto search...")
      print("The following rows did not match:")

      if(length(test1)> 0)
      {
        print("biome and cstor:")
        print(cstor[test1,])
      }
      if(length(test2)> 0)
      {
        print("water and cstor:")
        print(cstor[test2,])
      }
      if(length(test3)> 0)
      {
        print("food and cstor:")
        print(cstor[test3,])
      }
      if(length(test4)> 0)
      {
        print("forage and cstor:")
        print(cstor[test4,])
      }
      next
    }

    # sum(biome$LONLAT %in% cstor$LONLAT) - nrow(biome)
    # sum(water$LONLAT %in% cstor$LONLAT) - nrow(water)
    # sum(food$LONLAT %in% cstor$LONLAT) - nrow(food)

    idx_nonmatch <- which(!(protectedAreas$LONLAT %in% food$LONLAT))

    if (length(idx_nonmatch) > 0)
    {
      print("Protected area cells without a match in food")
      print(paste("for rcp", aRcp, "and time period", aTimePeriod))
      print(dim(protectedAreas[idx_nonmatch,]))
    }
    # now (dec 2019) 14794 - 9412 =  5382


    # being conservative, join tables first
    idx_matched <- match(food$LONLAT, protectedAreas$LONLAT )
    protectedAreas <- protectedAreas[idx_matched,]
    dim(protectedAreas)
    dim(food)

    ###
    # now calculate the ES adjusted for the protected areas
    ###

    idx_prot <- which(protectedAreas$Protected > 0)
    foodCor <- food
    cstorCor <- cstor
    waterCor <- water
    forageCor <- forage

    # for each protected cell, assume that the protected part is covered with PNV
    for(aCropName in names(foodCor)[4:12])
    {
      foodCor[idx_prot, aCropName] <- foodCor[idx_prot, aCropName] * (1- protectedAreas$Protected[idx_prot])

      cstorCor[idx_prot, aCropName] <- cstorCor[idx_prot, aCropName] * (1- protectedAreas$Protected[idx_prot]) + cstorCor$PNV[idx_prot]* protectedAreas$Protected[idx_prot]
      waterCor[idx_prot, aCropName] <- waterCor[idx_prot, aCropName] * (1- protectedAreas$Protected[idx_prot]) + waterCor$PNV[idx_prot]* protectedAreas$Protected[idx_prot]
      forageCor[idx_prot, aCropName] <- forageCor[idx_prot, aCropName] * (1- protectedAreas$Protected[idx_prot]) + forageCor$PNV[idx_prot]* protectedAreas$Protected[idx_prot]
    }

    summary(foodCor)
    summary(food)

    # write to disk for use in python
    # add ID for merge in python
    cstorCor$ID <- 1:nrow(cstorCor)
    waterCor$ID <- 1:nrow(waterCor)
    foodCor$ID <- 1:nrow(foodCor)
    forageCor$ID <- 1:nrow(forageCor)
    # drop LatLon column befor export
    (idxLatLon <- which(names(cstorCor) == "LONLAT"))
    write.table(cstorCor[,-idxLatLon], paste0( baseFn, "CStor_cor4ProtArea.txt"), row.names=FALSE, sep="\t", quote = FALSE)
    write.table(waterCor[,-idxLatLon], paste0( baseFn,"WaSu_cor4ProtArea.txt"), row.names=FALSE, sep="\t", quote = FALSE)
    write.table(foodCor[,-idxLatLon], paste0( baseFn,"FoPro_cor4ProtArea.txt"), row.names=FALSE, sep="\t", quote = FALSE)
    write.table(forageCor[,-idxLatLon], paste0( baseFn,"Forage_cor4ProtArea.txt"), row.names=FALSE, sep="\t", quote = FALSE)
    # now the real optimization

    # first create the data frame with all possible combinations
    colOffset <- 2

    thePfts <- sort(unique(pfts$ID))
    theBiomes <- sort(unique(biome$ATbiome))

    nPft <- length(thePfts)
    nBiome <- length(theBiomes)

    theCombinations <- expand.grid(biome1=thePfts, biome2=thePfts, biome3=thePfts, biome4=thePfts, biome5=thePfts, biome6=thePfts, biome7=thePfts, biome8=thePfts)

    theCombinations$Cstorage <- 0
    theCombinations$Water_supply <- 0
    theCombinations$Food_provisioning <- 0
    theCombinations$Forage <- 0
    # the first element contains the indices for all cells that belong to biome 1
    # the second the indices for biome 2 and so on
    biomesIdx <- lapply(1:8, FUN= function(x) which(biome$ATbiome == x))
    # number of cells that belong to biomes
    nCellsByBiome <- sapply(biomesIdx, FUN=function(x) length(x))
    data.frame(biome.legend$biomName, nCellsByBiome)
    sum(nCellsByBiome)

    # calculate the ES provisioning for the different combinations
    startTime <- Sys.time()
    for(aBiome in theBiomes)
    {
      for(aPft in thePfts)
      {
        # select all entries in theCombinations which have this Pft in that biom

        idx_Combinations <-which(theCombinations[,aBiome] == aPft)
        theIndicesOfTheCurrentBiome <- biomesIdx[[aBiome]]

        theCombinations$Cstorage[idx_Combinations] <- theCombinations$Cstorage[idx_Combinations] + sum(cstorCor[theIndicesOfTheCurrentBiome, colOffset + aPft+1])
        theCombinations$Water_supply[idx_Combinations] <- theCombinations$Water_supply[idx_Combinations] + sum(waterCor[theIndicesOfTheCurrentBiome, colOffset + aPft+1])
        theCombinations$Food_provisioning[idx_Combinations] <- theCombinations$Food_provisioning[idx_Combinations] + sum(foodCor[theIndicesOfTheCurrentBiome, colOffset + aPft+1])
        theCombinations$Forage[idx_Combinations] <- theCombinations$Forage[idx_Combinations] + sum(forageCor[theIndicesOfTheCurrentBiome, colOffset + aPft+1])

      }
    }
    endTime <- Sys.time()


    print(paste("Objectives for all biome pft combination ready. Time required: ", round(endTime - startTime,2), "sec."))
    print(paste("for rcp", aRcp, "and time period", aTimePeriod))

    dim(theCombinations)

    summaryCombinations <- summary(theCombinations)
    summaryCombinations[,-(1:8)]
    # select only those Combinations that are not to far from the reference situation
    selectFactor <- 0.9
    idx_suitable <- which( (theCombinations$Cstorage > selectFactor * constraints$carbon * 10^3) & (theCombinations$Food_provisioning > selectFactor* constraints$food*10^9) & (theCombinations$Water_supply > selectFactor* constraints$water*10^3) &  (theCombinations$Forage > selectFactor* constraints$forage*10^3))
    length(idx_suitable)
    length(idx_suitable) / nrow(theCombinations)

    theCombinationsSuited <- theCombinations[idx_suitable,]
    save(theCombinationsSuited, file = paste(sep="", outFolder, "theCombinations_rcp", aRcp, "_", aTimePeriod, "_correctedProtected.Rdata"))

    print(paste0("Number of solutions that are better than ", selectFactor * 100, "% of all constraints: ", length(idx_suitable)))
    # now identify pareto front
    # step I

    startTime <- Sys.time()
    # new and important: multiply by -1 since nds_rank is minimzing...
    obj <- -1* as.matrix(theCombinationsSuited[, c("Cstorage", "Water_supply", "Food_provisioning", "Forage")])

    summary(obj / 10^6)

    # memory problems, split problem into parts
    getFront0Idx <- function(mat)
    {
      ranks <-nds_rank(t(mat))
      idx_ranks <- which(ranks ==1)
      return(idx_ranks)
    }
    #nds_rank(t(maxgen.mat))

    # 10^6 seems to be a maximum, better 10^5
    idx_optimal <- NULL
    splitPoints <- seq(1, nrow(theCombinationsSuited), by=10^4)
    length(splitPoints)

    print("First round of optimization")
    for(i in 2:length(splitPoints))
    {
      if (i %%25 ==0)
        print(paste(round(i / length(splitPoints) * 100 ,1), "% of parts done"))
      startSplitIdx <- splitPoints[i-1]
      endSplitIdx <- splitPoints[i] -1

      obj.sub <- obj[startSplitIdx:endSplitIdx,]
      idx_tmp <- getFront0Idx(obj.sub)
      idx_optimal <- c(idx_optimal, idx_tmp + startSplitIdx -1) # bug fixed! here -1 was missing, which resulted in an offset of 1 for the index..., fixed 13.08.2018

      #save(idx_optimal,  file = paste(sep="", outFolder, "idx_optimal.Rdata"))
    }

    length(idx_optimal)
    # make a copy that contains only the pareto-optimal solutions from the split parts of the whole combinations
    obj.optimParts <- obj[idx_optimal,]

    save(obj.optimParts, file = paste(sep="", outFolder, "obj_sub_", aRcp, "_", aTimePeriod, ".Rdata"))


    #########################
    # a second round

    print("Starting second phase of Pareto search")
    theCombination.pareto.sub <- theCombinationsSuited[idx_optimal,]

    print(paste("Number of solutions to be considered:", dim(theCombination.pareto.sub) ) )

    summary(theCombination.pareto.sub[, c("Cstorage", "Water_supply", "Food_provisioning", "Forage")] / 10^6)

    obj <- -1* as.matrix(theCombination.pareto.sub[, c("Cstorage", "Water_supply", "Food_provisioning", "Forage")])

    idx_optimal <- NULL
    splitPoints <- seq(1, nrow(theCombination.pareto.sub), by=5*10^4)
    nSplitPoints <- length(splitPoints)

    if(nSplitPoints > 2)
    {
      for(i in 2:length(splitPoints))
      {
        if (i %%25 ==0)
          print(paste(round(i / nSplitPoints * 100, 1), "% of parts done"))
        startSplitIdx <- splitPoints[i-1]
        endSplitIdx <- splitPoints[i] -1

        obj.sub <- obj[startSplitIdx:endSplitIdx,]
        idx_tmp <- getFront0Idx(obj.sub)
        idx_optimal <- c(idx_optimal, idx_tmp + startSplitIdx -1) #m again, -1 was missing here, fixed 13.08.2018

        #save(idx_optimal,  file = paste(sep="", outFolder, "idx_optimal_2.Rdata"))
      }

      print(paste("Number of solutions at second stage identified: ", length(idx_optimal) ) )

      obj.sub <- obj[idx_optimal,]
      summary(obj.sub)
      dim(obj.sub)

      theCombination.paretoStage2 <- theCombination.pareto.sub[idx_optimal,]
      obj <- -1* as.matrix(theCombination.paretoStage2[, c("Cstorage", "Water_supply", "Food_provisioning", "Forage")])

      save(theCombination.paretoStage2, file = paste(sep="", outFolder, "theCombination_paretoStage2_", aRcp, "_", aTimePeriod, ".Rdata"))
    } else {
      theCombination.paretoStage2 <- theCombination.pareto.sub
      obj <- -1* as.matrix(theCombination.paretoStage2[, c("Cstorage", "Water_supply", "Food_provisioning", "Forage")])
    }




    ##################################
    # final round

    print("Final stage of pareto front search...")

    idx_pareto <- getFront0Idx(obj)
    print(paste("Number of solutions at final stage: ", length(idx_pareto)))

    print(paste("Sum of ranks (should be equal to number of solutions):", sum(nds_rank(t(obj[idx_pareto,]))) ))

    theCombination.paretoFinal <- theCombination.paretoStage2[idx_pareto,]
    # beeing conservative, testing again if really pareto optimal
    #obj.pareto.test <- -1* as.matrix(theCombination.pareto[, c("Cstorage", "Water_supply", "Food_provisioning")])
    #sum(nds_rank(t(obj.pareto.test)))
    #dim(obj.pareto.test)
    # yes they are
    endTime <- Sys.time()
    print(paste("Time needed for all stages of pareto front search:", round(endTime-startTime,2), "sec."))

    save(theCombination.paretoFinal, file = paste(sep="", outFolder, "theCombination_paretoBiomeFinal_", aRcp, "_", aTimePeriod, ".Rdata"))

    write.table(theCombination.paretoFinal, file = paste(sep="", outFolder, "theCombination_paretoBiomeFinal_", aRcp, "_", aTimePeriod, ".csv"),  row.names = FALSE, sep=";")
  }
}

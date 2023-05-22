# Data

## Ecosystem sercive values - ESquant folder

Contains the ecosystem services calculatged by LPJ-GUESS. Each file contains the ecosystem service value for each grid cell, for each scenario for each land use vegetation class.

File naming convention:

ESquant_rcp*rcpID*\_*time-horizon*\_INT\_9412\_v8b\_4gcms\_*ecosystemservice*.txt

For example
ESquant_rcp26_2033-2042_INT_9412_v8b_4gcms_CStor.txt

Options 

  - RCP: *26* | *60*
  - time-horizon: *2033-2042* | *2090-2099*
  -ecosystem service: *CStor* | *FoPro* | *WaSu* | *Forage* (carbon storang, food provisioning, water supply, forage provisioning)
  
 The preprocessing scripts correct the values for those cells that contain protected areas and svae results in a new file with added "*cor4ProtArea*" to the file name.
 
 File format: tabulator (or white space) separated values, first line contains header:
      Lon      Lat             PNV      C3cer_TeWW     C3cer_TeWWi     C3crps_TeSW    C3crps_TeSWi     C4crps_TeCo    C4crps_TeCoi       Rice_TrRi      Rice_TrRii             Pas
 Longitude, Lattitude of the cell center (WGS-1984), followed by ecosystem service value per cell for the following crop functional types:
 
   - *PNV* - potential natural vegetation
   - *C3cer\_TeWW*, C3 cereals, rainfeed
   - *C3cer\_TeWWi*, C3 cereals, irrigated
   - *C3crps\_TeSW*, non-C3 cereals, rainfeed
   - *C3crps\_TeSWi*, non-C3 cereals, irrigated
   - *C4crps\_TeCo*, corn, rainfeed
   - *C4crps\_TeCoi*, corn, irrigated
   - *Rice\_TrRi*, rice, rainfeed
   - *Rice\_TrRii* - rice, irrigated
   - *Pas* - pasture

## Biome

- *biome.txt* - provides for each grid cell the information to which biome it belongs
- *legend_biomes_port_AT.txt* biomeID, desctiption and RGB bvalues for legend

## Protected areas

- gridlist_hurtt_RNDM_1deg_protarea.txt - lattitude, longitude, share of area protected
## Optimization constraints

- contraints.txt - for each time-horizon and RCP the ecosystem serive provisioning for the current land use configuration

## Link between food producing units and grid cells

- fpuCentroids.txt and fpuCentroids.Rdata - the food producing unit ID (fpu) for each grid cell as a lattitude, longitude, fpuID table, for *fpu_final.txt* in addition the combination of lattitue and longitude with the slash as an separator
- fpu2biomeGrid.csv - for each grid cell the information on biome (ATbiome), fpu (FPU_CODE, original and newFPU after intersection with biomes), grid cell ID (ID) and the lattitude and longitude information
- 
- fpu2biomeTable.csv

## Crop functional types

- abbrevation, color code and short description of the crop functional type as used by the optimization

## Biodiversity

- *Pouzols_global2000_int1deg.txt* the biodiversity information from Pouzols et al.



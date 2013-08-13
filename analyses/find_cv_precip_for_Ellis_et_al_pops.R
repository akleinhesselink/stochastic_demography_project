#### Ellis et al. 2012, get climate data

rm(list = ls() )
#install.packages('raster')
#install.packages('SDMTools')
library(raster)
library(SDMTools)

setwd("~/Documents/Courses/Matrix Population Models/final_project/analyses/")

species_info = read.table('Ellis_et_al_Species_Information.txt', sep = '\t', header = T)
species_info

pop_info = read.table('Ellis_et_al_Population_data.txt', sep = '\t', header = T)

######## ERROR in LAT_LONGS of sites in the DATA!!!!!!!!!!! 
######## FOR TRGR and for LAVE I need to switch lat longs 
TRGR_lat_longs = pop_info[pop_info$SPP == 'TRGR', c(4,3)]   #### switch lat and long
TRGR_lat_longs[, 2] = -TRGR_lat_longs[, 2]                  #### Set the longitude negative 
pop_info[pop_info$SPP == 'TRGR', c(3,4)] = TRGR_lat_longs   #### reset the values in the DF
##########

######## ERROR on HYCU_91, population as well, longitude sould be negative
pop_info[pop_info$POP == 'HYCU_91', 4] = - pop_info[pop_info$POP == 'HYCU_91', 4]
########

pop_info = pop_info[order(pop_info$SPP, pop_info$POP), ]

pop_locales = cbind(pop_info$Long, pop_info$Lat)
pop_locales

prec_cv_data = raster("prec_var_values.tif")
locales = SpatialPoints(pop_locales, proj4string = CRS(projection(prec_cv_data)) )

plot(prec_cv_data)
points(locales)
cv_vals = extract(prec_cv_data, locales)
pop_info[cv_vals < -100, ]
cv_vals[cv_vals <  -20] = (17.15+17.65)/2 #### Set the missing Baltic SEA to average of the two values
prec_cv_at_all_sites = cbind(pop_info[, c(1,2,3,4)], cv_vals)
prec_cv_at_all_sites
write.csv(prec_cv_at_all_sites, "prec_cv_data.csv", row.names = F)



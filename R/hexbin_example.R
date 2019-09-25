library(raster)
library(sf)
library(exactextractr)
library(dggridR)
library(geojsonsf)

# Pull municipal boundaries for Brazil
sweden_st <- getData('GADM', country='SWE', level=0)
sweden <- st_as_sf(sweden_st)

# Write to shape file
st_write(sweden, "~/Downloads/sweden.shp", delete_dsn=TRUE)

# Pull gridded spp data
s <- stack("/home/edward/Downloads/proj_rcp45_2020_Pinus.mugo_ensemble.grd")
r <- raster(s, layer=6) 
crs(r) <- CRS("+init=epsg:4326") 

crs(r)
res(r)
#newproj <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"
#r2 <- projectRaster(r, crs=newproj, res=10000)
#crs(r2)
#res(r2)
#r2
stk <- stack(list(ensemble_mean=r, area=area(r)))


# Create hexgrid
dggs      <- dggridR::dgconstruct(res=11, metric=TRUE, resround='nearest', show_info = T)
hex_sweden <- dgshptogrid(dggs, "~/Downloads/sweden.shp", frame=F)
hs <- st_as_sf(hex_sweden)


#prec <- getData('worldclim', var='prec', res=10)[[12]]

# Find the mean precipitation amount for each municipality
hs$weighted_ensemble_mean <- exact_extract(stk, hs, function(values, coverage_frac)
  weighted.mean(values$ensemble_mean, values$area*coverage_frac, na.rm=TRUE))

st_write(hs, "~/Downloads/sweden_pinus_wem_hex.geojson", delete_dsn=TRUE) 

# Find min and max precipitation amount in a single pass
#hs[, c('min_prec', 'max_prec')] <- exact_extract(prec, hs, c('min', 'max'))

# Plot
require(tmap)
tmap_mode("view")
## tmap mode set to interactive viewing
#tm_shape(sweden) + tm_borders("white", lwd = .5) + 
tm_shape(hs) + 
tm_fill("weighted_ensemble_mean", palette = sf.colors(3)) +
tm_legend(show = FALSE) +  
tm_shape(sweden) +
tm_borders("black", lwd = .5) +
tm_legend(show = FALSE)


#sf::plot(hs["mean_prec"])

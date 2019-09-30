require(vspt)
require(raster)

bv <- create_biovar_stack(
  region = 'biovar_Sweden',
  biovar_list = c("biovar04", "biovar06", "biovar16"),
  time_period = "biovar_historical",
  path = "inst/ext_data/bioclimatic"
)
#plot(s$biovar_stack[[1]])
#plot(s$land_poly)

r <- bv$biovar_stack
raster_object <- raster::stack(r, raster::area(r))
names(raster_object) <- c(names(r), 'area')

polygon_object <- create_grid_from_poly(bv$land_poly, grid_res=11, set_seed=999)
names(polygon_object)

# Define function to apply
fun <- function(values){
  nms <- names(values)
  nms <- names(values)[1:(length(nms) -2)]
  return(sapply(nms, function(nm){
      if(!all(is.na(values[,nm]))){
        out <- weighted.mean(values[,nm], values$area*values$coverage_frac, na.rm=T)}
      else{
        out <- NA
        }
  }))
  }

zsp <- zonal_stats_poly(fun,
                             polygon_object,
                             raster_object,
                             remove_na = T,
                             out_path=F)


names(zsp)

if(require(tmap)){
  tmap_mode("view")
  pv <- c("biovar04","biovar06")
  #adjustcolor( '#f7f7f7', alpha.f = 0.1)
  cs <- rev(colorRampPalette(c('#ef8a62','#F7F7F71A','#67a9cf'))(5))
  alpha <- 0.7
  tm_shape(zsp) +
    tm_borders(alpha=0.1) +
    tm_fill("biovar04"
            , style="fixed"
            , breaks=c(200, 300, 400, 500, 600, 700, 800, 900, 1000)
            , alpha=alpha, palette=cs, popup.vars = pv ) +
    #tm_shape(zsp) +
    #tm_borders(alpha=0.5) +
    #tm_fill("biovar06"
    #        , style="fixed"
    #        , breaks = c(-15.0, -10.0, -5.0, 0.0, 5.0)
    #        , alpha=alpha, palette=cs, popup.vars = pv) +
    tm_scale_bar()

   }else{plot(out)}



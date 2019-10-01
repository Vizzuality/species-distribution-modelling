require(vspt)
require(raster)

# Define lists
region='biovar_Sweden'
biovar_list= c("biovar01","biovar02","biovar03","biovar04","biovar05","biovar06",
               "biovar07","biovar08","biovar09","biovar10","biovar11","biovar12",
               "biovar13","biovar14","biovar15","biovar16")
scenario_list=c("rcp45", "rcp85", "rcp26", "rcp60")


# Create biovar stack object (everything bioclimatic we need for the modelling)
bvso <- create_biovar_stack(region, biovar_list, scenario_list,
  path="inst/ext_data/bioclimatic",
  mask_by_land = T,
  pca_transform = T
)
names(bvso)
names(bvso[[1]])
save(bvso, file="inst/ext_data/SWE_bioclimatic.rda")


# Create vector polygon grid object
load("inst/ext_data/SWE_bioclimatic.rda")
vpgo <- create_grid_from_poly(bvso[[1]]$land_poly, grid_res=11, set_seed=999)
names(vpgo)

# Do zonal statistics on bioclimatic rasters
do_zs <- function(so){
  soa <- raster::stack(so, raster::area(r))
  names(soa) <- c(names(r), 'area')
  # Define function to apply
  f1 <- function(values){
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
  return(
    zonal_stats_poly(fi, vpgo, so, remove_na = T, return_df = T, out_path=F)
  )

}
zstat <- list()
zstat[[1]]$historical <- do_zs(bvso$biovar_Sweden$historical$raw)


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



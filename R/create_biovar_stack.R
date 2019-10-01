#' Create masked bioclimatic variable raster stack
#'
#' Select historic and future  bioclimatic variable rasters from directory,
#' create ensembles, create a stacks, optionally convert to principle components
#' and optionally mask by land. Note this function can take a while....
#'
#' @param region String representing the region directory
#' @param biovar_list List of string representing the biovar codes
#' @param path String representing the path to the bioclimatic directory
#'
#' @return
#' A named list with the region as first level containing lists of;
#' data file paths, processed historical and future biovariable rasters,
#' polygon of land, the biovariable and scenario lists used, and optionally
#' a pca model of the historical biovariables.
#' @export
#'
#' @examples
#' bvso <- create_biovar_stack(
#' region='biovar_Sweden',
#' biovar_list= c("biovar01","biovar02","biovar03","biovar04","biovar05","biovar06",
#'             "biovar07","biovar08","biovar09","biovar10","biovar11","biovar12",
#'                         "biovar13","biovar14","biovar15","biovar16"),
#'                         scenario_list=c("rcp45", "rcp85", "rcp26", "rcp60"),
#'                         path="inst/ext_data/bioclimatic",
#'                         mask_by_land = T,
#'                         pca_transform = T
#'                         )
#'                         names(bvso)
#'                         names(bvso[[1]])
create_biovar_stack <- function(
  region='biovar_Sweden',
  biovar_list= c("biovar01","biovar02","biovar03","biovar04","biovar05","biovar06",
                 "biovar07","biovar08","biovar09","biovar10","biovar11","biovar12",
                 "biovar13","biovar14","biovar15","biovar16"),
  scenario_list=c("rcp45", "rcp85", "rcp26", "rcp60"),
  path="inst/ext_data/bioclimatic",
  mask_by_land = T,
  pca_transform = T
){
  # Set raster::rater options to show progress as text
  rasterOptions(progress = 'text',timer=TRUE)

  # Start cluster
  raster::beginCluster()

  # Create the output biovar stacks object for a single region
  bvso <- list(
    list(
      df_file_paths = NA
      , historical = list(biovar=NA, pc=NA, raw_masked=NA, pc_masked=NA)
      , future = list(biovar=NA, pc=NA, raw_masked=NA, pc_masked=NA)
      , pca_model = NA
      , land_poly = NA
      , biovar_list = biovar_list
    , scenario_list = scenario_list
    )
  )
  names(bvso) <- region

  # Parse the bioclimatic directory
  print(paste('Parsing bioclimatic directory:', path))
  df_bc <- vspt::parse_bioclimatic(path=path)
  df <- df_bc[[region]]
  bvso[[1]]$df_file_paths <- df

  # Historical
  dfh <- df[df$time_period == 'biovar_historical' & df$bio_variable %in% biovar_list,]

  # Future
  dff <- df[df$time_period == 'biovar_future' &
              df$bio_variable %in% biovar_list &
              df$scenario %in% scenario_list,]

  # Create raster stacks. For future create enemsbles and then stacks as a list
  process_future <- function(dff){
    fl <- split(dff[, "file_name"], dff[, c("time_interval", "scenario", "bio_variable")])
    # make ensemble mean per year, scenario, biovar
    s_list <- lapply(fl, function(x){
      # create stack
      out <- raster::stack(unlist(x))
      raster::crs(out) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
      # ensemble mean

      out <- raster::clusterR(out, raster::calc, args=list(fun=mean))
      return(out)
    })

    # make stacks per year per scenario
    sdf <- data.frame(matrix(unlist(strsplit(names(s_list), "[.]")), ncol = 3, byrow = T))
    sdf <- sdf[order(sdf$X1, sdf$X2, sdf$X3),]
    ss <- list()
    years <- unique(sdf[,1])
    scenes <- unique(sdf[,2])
    for(y in years){
      for(s in scenes){
        sl <- paste(y, s , biovar_list, sep=".")
        #print(sl)
        #print(na.exclude(match(names(s_list), sl, incomparables = F)))
        ssl <-  s_list[names(s_list) %in% sl]
        st <- raster::stack(ssl)
        ss[[y]][[s]] <- list(st)

        }
    }
    return(ss)
  }
  process_historical <- function(dfh){
    fl <- unlist(dfh$file_name)
    s <- raster::stack(fl)
    raster::crs(s) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    names(s) <- biovar_list
    return(s)
  }

  print('Creating historical biovar ensembles and raster stack(s)')
  bvso[[1]]$historical$raw <- process_historical(dfh)
  print('Creating future biovar raster stack(s)')
  bvso[[1]]$future$raw <- process_future(dff)

  # Optionally convert to principle components
  if(pca_transform){
    print('Converting all stacks to principle components')
    bvso[[1]]$pca_model <- prcomp(raster::values(bvso[[1]]$historical$raw), scale=T, tol = 0.5)
    bvso[[1]]$historical$pca <- raster::clusterR(bvso[[1]]$historical$raw,
                                             raster::predict,
                                             args = list(model=bvso[[1]]$pca_model,
                                                         index=1:ncol(bvso[[1]]$pca_model$x))
    )
    names(bvso[[1]]$historical$pca) <- colnames(bvso[[1]]$pca_model$x)
    bvso[[1]]$future$pca <- lapply(
      bvso[[1]]$future$raw,
      function(l){
        return(
          lapply(l,
                 function(x){
                   #print(class(x[[1]]))
                   names(x[[1]]) <- bvso[[1]]$biovar_list
                   out <- raster::clusterR(x[[1]],
                                    raster::predict,
                                    args = list(model=bvso[[1]]$pca_model,
                                                index=1:ncol(bvso[[1]]$pca_model$x))
                                    )
                   names(out) <- colnames(bvso[[1]]$pca_model$x)
                   return(out)
                   }
                 )
        )
      }
    )
  }

  # Add land polygon
  p <- sf::st_read("inst/ext_data/boundary-vectors/ne_10m_land/ne_10m_land.shp")
  bvso[[1]]$land_poly <- sf::st_crop(p, sf::st_bbox(bvso[[1]]$historical$raw))

  if(mask_by_land){
    print('Masking all stacks by land')
    bvso[[1]]$historical$raw_masked <- raster::mask(bvso[[1]]$historical$raw, bvso[[1]]$land_poly)
    try(
      bvso[[1]]$historical$pc_masked <- raster::mask(bvso[[1]]$historical$pc, bvso[[1]]$land_poly)
    )
    bvso[[1]]$future$raw_masked <- lapply(bvso[[1]]$future$raw, function(l){
      return(
        lapply(l, function(x) raster::mask(x, bvso[[1]]$land_poly))
      )})
    try(bvso[[1]]$future$pc_masked <- lapply(bvso[[1]]$future$pc, function(l){
      return(
        lapply(l, function(x) raster::mask(x, bvso[[1]]$land_poly))
      )}))
      }
  # End cluster
  raster::endCluster()
  return(bvso)
}

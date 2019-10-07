#' Zonal stats for each polygon feature
#'
#' @param fun Function to apply to values in each polygon
#' @param vpgo A sf polygon feature collection
#' @param ro A raster::raster object
#' @param remove_na Boolean, remove features with NA values?
#' @param out_path Write object to this path
#'
#' @return
#' A sf::sf object identical to the input polugon object with zonal statistics added
#' or a dataframe of the zonal statistics.
#' @export
#'
#' @examples
#' # add example
zonal_stats_poly <- function(fun,
                             vpgo,
                             ro,
                             remove_na = T,
                             return_df = F,
                             out_path=F
){
  nmsr <- names(ro)
  bb <- sf::st_bbox(ro)

  # Crop polygon object to raster bbox
  crs_poly <-  sf::st_crs(vpgo)
  hg <- sf::st_crop(vpgo, bb)
  sf::st_crs(hg) <- crs_poly
  nmsp <- names(hg)

  # Get zonal stats for polygons
  zstat <- exactextractr::exact_extract(x=ro, y=hg, fun=NULL)
  if(class(zstat)=='list'){zstat <- as.data.frame(do.call(rbind, lapply(zstat, fun)))}

  # Add results to polygon
  if(!return_df){out <- sf::st_sf(data.frame(hg, zstat))}else{out <- zstat}
  #names(out) <-  c(nmsp, paste0(c('zstat', nmsr), collapse = '_'))

  # Remove NA values
  if(remove_na == T){ out <- na.omit(out) }

  # Optionally write to file
  if(out_path != F){ sf::st_write(out, out_path, delete_dsn=TRUE) }

  return(out)
}

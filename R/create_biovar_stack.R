#' Create masked bioclimatic variable raster stack
#'
#' Select bioclimatic variable rasters from directory, create a stack and mask by land
#'
#' @param region String representing the region directory
#' @param biovar_list List of string representing the biovar codes
#' @param time_period String representing the time-period code
#' @param path String representing the path to the bioclimatic directory
#'
#' @return
#' A list with 'biovar_stack' and 'land_poly' representing the raster::stack
#' of the selected bioclimatic variables, for the region and time-period, masked by land mass
#' , and a polygon of land, respectively.
#' @export
#'
#' @examples
#' s <- create_biovar_stack(region='biovar_Sweden',
#' biovar_list= c("biovar04","biovar06","biovar16"),
#' time_period="biovar_historical",
#' path="inst/ext_data/bioclimatic")
#' raster::plot(s$biovar_stack)
create_biovar_stack <- function(
  region='biovar_Sweden',
  biovar_list= c("biovar04","biovar06","biovar16"),
  time_period="biovar_historical",
  path="inst/ext_data/bioclimatic"
){

  # Parse the bioclimatic directory
  df_bc <- vspt::parse_bioclimatic(path=path)

  # Create stack
  df <- df_bc[[region]]
  df <- df[df$bio_variable %in% biovar_list &df$time_period == time_period,]
  df <- df[df$time_period == time_period,]
  df$scenario <- NULL
  fl <- unlist(file.path(path, apply(df[,c("biovar_region", "time_period", "file_name")], 1, paste0, collapse = "/")))
  s <- raster::stack(fl)
  raster::crs(s) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  bb <- sf::st_bbox(s)
  names(s) <- biovar_list

  # Mask raster by land
  # get global land
  poly <- sf::st_read("inst/ext_data/boundary-vectors/ne_10m_land/ne_10m_land.shp")
  poly_cropped <- sf::st_crop(poly, bb)

  # mask
  s_masked <- raster::mask(s, poly_cropped)
  #raster::plot(s_masked)

  return(list(biovar_stack=s_masked, land_poly=poly_cropped))
}

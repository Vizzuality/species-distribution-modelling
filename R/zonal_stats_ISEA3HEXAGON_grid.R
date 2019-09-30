#' Get zonal stats binned in a hexagon grid
#'
#' Given a geometry (multi) polygon object, create a ISEA-3-HEXAGON grid at the given resolution
#' covering the geometry, apply the zonal statistics function(s) to each feature of the grid
#' feature collection, and return the grid with calculated properties. Optionally clipping the grid
#' by the geometry object, and saving to file.
#'
#' @param fun An optional function or character vector, as described in exact_extract {exactextractr}
#' @param polygon_object A sf::sf object
#' @param raster_object A raster::raster object
#' @param grid_res The zoom resolution of the ISEA-3-HEXAGON grid, see dggridR::dgconstruct
#' @param clip Boolean, clip by the polygon? Note this can be slow
#' @param out_path Path to write the grid using sf::st_write
#'
#' @return
#' A sf object representing the ISEA-3-HEXAGON grid covering the geometry with zonal statistics as
#' defined by fun
#' @export
#'
#' @examples
#' polygon_path <- system.file("ext_data", "boundary-vectors", "sweden.shp", package="vspt")
#' polygon_object <- sf::st_as_sf(sf::st_read(polygon_path))
#' raster_path <- system.file("ext_data/predicted-spp-occurence-rasters/proj_rcp45_2020_Pinus.mugo_ensemble.grd", package="vspt")
#' raster_object <- raster::raster(raster_path, band=6, crs="+proj=longlat +datum=WGS84")
#' raster_object
#' # Make stack with area
#' raster_object <- raster::stack(raster_object, raster::area(raster_object))
#' out_path <- file.path(system.file("ext_data/predicted-spp-occurence-vectors", package="vspt"), "SWE_pinus-mugo_rcp45-2020_ensemble.geojson")
#' # Define function to apply
#' area_weighted_mean <- function(values, coverage_frac){return(weighted.mean(values[,1], values[,2]*coverage_frac, na.rm=TRUE))}
#' fun = area_weighted_mean
#' out <- zonal_stats_ISEA3HEXAGON_grid(fun, polygon_object,raster_object,grid_res = 11, remove_na = T, clip=F, out_path=out_path)
#' names(out)
#' if(require(tmap)){
#'   tmap_mode("view")
#'   tm_shape(out) + tm_fill("zstat_Pinus.mugo_EMcaByROC_mergedAlgo_mergedRun_mergedData", palette = sf::sf.colors(3), alpha=0.7, colorNA=NULL)
#'   }else{plot(out)}
zonal_stats_ISEA3HEXAGON_grid <- function(fun,
                                          polygon_object,
                                          raster_object,
                                          grid_res = 11,
                                          remove_na = T,
                                          clip=F,
                                          out_path=F
                                          ){
  # Convert raster object into raster object
  #raster_object <- raster::raster(raster_object)
  #raster::crs(raster_object) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  #plot(raster_object)
  nm <- names(raster_object)[1]
  bb <- sf::st_bbox(raster_object)

  # Convert polygon object into sf object
  crs_poly <-  sf::st_crs(polygon_object)
  poly <- sf::st_crop(polygon_object, bb)
  sf::st_crs(polygon_object) <- crs_poly

  # Convert to temporary shape file
  # required as dggs only accepts shape files
  f <- file.path(tempdir(), "poly.shp")
  sf::st_write(poly, f, delete_dsn=TRUE)

  # Create ISEA-3-HEXAGON grid
  hg <- vspt::create_grid_from_poly(polygon_object = poly, grid_res = 11, set_seed = 666)

  # Important to ensure bbox of hg and raster match!
  hg <- sf::st_crop(hg, bb)

  # Get zonal stats for grid
  zstat <- exactextractr::exact_extract(x=raster_object, y=hg, fun=fun)

  # Add results to polygon
  out <- sf::st_sf(data.frame(hg, zstat))
  names(out) <-  c('uuid', 'area_m2', 'perimeter_m', paste0(c('zstat', nm), collapse = '_'), 'geometry')

  # Remove NA values
  if(remove_na == T){ out <- na.omit(out) }

  # Optionally clip grid
  # This is silly slow!!!
  if(clip == T){ out <- sf::st_as_sf(sf::st_intersection(out, poly)) }

  # Optionally write to file
  if(out_path != F){ sf::st_write(out, out_path, delete_dsn=TRUE) }

  return(out)
}

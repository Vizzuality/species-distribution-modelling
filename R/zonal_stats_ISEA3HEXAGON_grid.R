#' Get zonal stats binned in a hexagon grid
#'
#' Given a geometry (multi) polygon object, create a ISEA-3-HEXAGON grid at the given resolution
#' covering the geometry, apply the zonal statistics function(s) to each feature of the grid
#' feature collection, and return the grid with calculated properties. Optionally clipping the grid
#' by the geometry object, and saving to file.
#'
#' @param fun An optional function or character vector, as described in exact_extract {exactextractr}
#' @param polygon_path Path to geometry file compatible with sf::st_read
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
#' fun = c("mean", "mode")
#' polygon_path <- system.file("ext_data", "boundary-vectors", "sweden.shp", package="vspt")
#' # FIXME: why does it not read via system.file?
#' raster_path <- "/home/edward/Downloads/proj_rcp45_2020_Pinus.mugo_ensemble.grd"
#' #system.file("ext_data/predicted-spp-occurence-rasters/proj_rcp45_2020_Pinus.mugo_ensemble.grd", package="vspt")
#' raster_object <- raster::raster(raster_path, band=6, crs="+proj=longlat +datum=WGS84")
#' raster_object
#' out_path <- file.path(system.file("ext_data/predicted-spp-occurence-vectors", package="vspt"), "SWE_pinus-mugo_rcp45-2020_ensemble.geojson")
#' out <- zonal_stats_ISEA3HEXAGON_grid(fun,polygon_path,raster_object,grid_res = 11, remove_na = T, clip=F, out_path=out_path)
#' if(require(tmap)){
#'   tmap_mode("view")
#'   tm_shape(out) + tm_fill("mode", palette = sf.colors(3), alpha=0.7, colorNA=NULL)
#'   }else{plot(out)}
zonal_stats_ISEA3HEXAGON_grid <- function(fun,
                                          polygon_path,
                                          raster_object,
                                          grid_res = 11,
                                          remove_na = T,
                                          clip=F,
                                          out_path=F
                                          ){

  # Read polygon path into sf object
  poly <- sf::st_read(polygon_path)

  # Convert to temporary shape file
  # required as dggs only accepts shape files
  f <- file.path(tempdir(), "poly.shp")
  sf::st_write(poly, f, delete_dsn=TRUE)

  # Create ISEA-3-HEXAGON grid definition
  dggs <- dggridR::dgconstruct(res=grid_res, metric=TRUE, resround='down', show_info = T)

  # Apply to polygon and convert to sf object
  hg <- sf::st_as_sf(
    dggridR::dgshptogrid(dggs, f, frame=F)
  )

  # Add random uuids to each feature
  hg$uuid <- apply(hg, 1, FUN=function(x){return(uuid::UUIDgenerate())})

  # Read raster path into raster object
  #raster_object <- raster::stack(raster_path)

  # Create stack with cell area
  # Not needed for equal area heaxagons?
  #s <- raster::stack(list(stat=r, area=area(r)))

  # Get zonal stats for grid
  zstat <- exactextractr::exact_extract(x=raster_object, y=hg, fun=fun)

  # Add results to polygon
  out <- sf::st_sf(data.frame(hg, zstat))

  # Remove NA values
  if(remove_na == T){ out <- na.omit(out) }

  # Optionally clip grid
  # This is silly slow!!!
  if(clip == T){ out <- sf::st_as_sf(sf::st_intersection(out, poly)) }

  # Optionally write to file
  if(out_path != F){ sf::st_write(out, out_path, delete_dsn=TRUE) }

  return(out)
}

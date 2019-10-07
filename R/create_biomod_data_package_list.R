#' Create biomod2 data package list
#'
#' TODO: set aside data for model evaluation!
#'
#' @param spp_occurence_raster_list List of species occurence rasters
#' @param biovar_stack raster::stack of bioclimatic explanationary values
#'
#' @return
#' A list of biomod2 data packages
#' @export
#'
#' @examples
create_biomod_data_package_list <- function(spp_occurence_raster_list, biovar_stack){

  return(
    lapply(spp_occurence_raster_list, function(x){

      # spp name
      spp_name <- names(x)

      # convert spp raster to points
      species <- as.data.frame(raster::rasterToPoints(x))
      presences <-
      resp_var <- as.numeric(species[which(species[,3]==1) , 3]) # species presences vector (only 'ones')
      resp_xy <- species[which(species[,3]==1) , 1:2] # coordinates of presence cells

      # create biomod data object
      return(
        biomod2::BIOMOD_FormatingData(
          resp.var = resp_var,
          expl.var = raster::stack(biovar_stack),
          resp.xy = resp_xy,
          resp.name = spp_name,
          PA.nb.rep = 2,
          PA.nb.absences = 300,
          PA.strategy = "sre",
          PA.sre.quant = 0.25)
      )

    })
  )
}

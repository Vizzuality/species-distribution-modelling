#' Project biomod ensemble model onto biovar stack
#'
#' @param biomod_ensemble_model_list List of calibrated biomod2 ensemble models
#' @param bimod_model_projection_list List of raster::stacks representing the model projections on biovariables
#'
#' @return
#' @export
#'
#' @examples
create_biomod_ensemble_model_projection_list <- function(biomod_ensemble_model_list, bimod_model_projection_list){
  return(
    mapply(function(x, y){
             spp_name <- names(x)
             return(
               biomod2::BIOMOD_EnsembleForecasting(
                 EM.output =  x,
                 projection.output = y,
                 binary.meth = "ROC",
                 filtered.meth = NULL)

             )
           }, biomod_ensemble_model_list, bimod_model_projection_list
    )
  )
}

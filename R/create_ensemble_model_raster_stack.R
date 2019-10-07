#' Create ensemble model projections raster stack
#'
#' @param spp_list List of species names used in modelling
#' @param proj_nm The projection name
#'
#' @return
#' A raster::stack with probalistic (CA, CV), and mean and median binary predictions.
#' @export
#'
#' @examples
create_ensemble_model_raster_stack <- function(spp_list, proj_nm){
  out_list <- lapply(spp_list, function(spp){
      fl <- list.files(paste0(gsub(" ",".", spp),"/proj_", proj_nm))
      # selecting the files for the consensus projections, those that have the word "ensemble"
      # in their name; the remaining ones in the same folder are the single.model projections
      # or the Clamping mask.
      enslist <- fl[grep("ensemble", fl)]
      probEns <- enslist[
        grep(paste0("proj_", proj_nm, "_", gsub(" ", ".", spp), "_ensemble.grd"), enslist)]
      ensp <- raster::stack(paste(gsub(" ", ".", spp),paste0("proj_", proj_nm), probEns, sep="/"))
      # ensp has the probabilistic projections, from here we want the ca and the cv
      EM.methods=c("EMca", "EMcv","EMmean", "EMmedian")
      ca_raster = raster::raster(ensp, grep(EM.methods[1], names(ensp)))
      cv_raster = raster::raster(ensp, grep(EM.methods[2], names(ensp)))
      binEns = enslist[
        grep(paste0("proj_", proj_nm, "_", gsub(" ", ".", spp), "_ensemble_ROCbin.grd"), enslist)]
      ensb <- raster::stack(paste(gsub(" ", ".", spp), paste0("proj_", proj_nm), binEns, sep="/"))
      # ensb has the binary projections, from here we want the mean and the median
      mean_bin_raster = raster::raster(ensb, grep(EM.methods[3], names(ensb)))
      median_bin_raster = raster::raster(ensb, grep(EM.methods[4], names(ensb)))
      out <- raster::stack(ca_raster, cv_raster, mean_bin_raster, median_bin_raster)
      names(out) <- EM.methods
      return(out)
    })
  names(out_list) <- spp_list
  return(out_list)

}

#' Create biomod2 ensemble of models
#'
#' @param biomod_model_list List of biomod2 model objects
#'
#' @return
#' A biomod2 ensemble model object
#' @export
#'
#' @examples
create_biomod_ensemble_model_list <- function(biomod_model_list){
  return(
    lapply(biomod_model_list,
           function(x){
             spp_name <- names(x)
             return(
               biomod2::BIOMOD_EnsembleModeling(
                 modeling.output = x,
                 chosen.models = "all",
                 em.by = "all", # 5 available options
                 VarImport = 1,
                 eval.metric = "ROC",
                 eval.metric.quality.threshold = NULL, # here we can set up a threshold to choose the models to use for the ensemble
                 models.eval.meth = c("TSS","ROC"),
                 prob.mean = T,
                 prob.median = T,
                 committee.averaging = T,
                 prob.mean.weight = F,
                 prob.mean.weight.decay = 'proportional',
                 prob.cv = T,
                 prob.ci = T,
                 prob.ci.alpha = 0.05)

             )
           })
  )
}

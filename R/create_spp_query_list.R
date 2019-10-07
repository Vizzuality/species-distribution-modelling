#' Create list of GBIF queries
#'
#' @param spp_list Species names
#' @param continent The continent
#' geom WTK polygon bounding box to limit search [THIS DOES NOT SEEM TO WORK]
#' @param creds gbif credentials
#' @param start_year Start year  to limit search
#' @param end_year End year to limit search
#'
#' @return
#' List of rgbif::occ_download_prep objects named by spp_list
#' @export
#'
#' @examples
#spp_list <- c('Pinus mugo')
#geom <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(bvso$biovar_Sweden$historical$raw)))
#create_spp_query_list(spp_list, geom, creds, start_year=1986, end_year=2020)
create_spp_query_list <- function(spp_list, continent="EUROPE", creds, start_year=1986, end_year=2020){
  keys <- lapply(spp_list, rgbif::name_backbone, rank='species')

  c_list <- lapply(keys, function(k){
    rgbif::occ_count(
      taxonKey = k$speciesKey
      , georeferenced = T
      , basisOfRecord = "HUMAN_OBSERVATION,OBSERVATION,MACHINE_OBSERVATION"
      , datasetKey = NULL
      , date = NULL
      , typeStatus = NULL
      , country = NULL
      , year = NULL
      , from = start_year
      , to = end_year
      , type = "countries"
      , publishingCountry = "US"
      , protocol = NULL
      , curlopts = list())
  })
  #c_list[[1]]$SWEDEN
  #c_list[[1]]$DENMARK
  #c_list[[1]]$NORWAY

   q_list <- lapply(keys, function(k){
     print(k$species)
     rgbif::occ_download_prep(
    paste("taxonKey =", k$speciesKey),
    "basisOfRecord = HUMAN_OBSERVATION,OBSERVATION,MACHINE_OBSERVATION",
    "hasCoordinate = true",
    paste("continent =", continent),
    #paste("geometry =", geom),
    "hasGeospatialIssue = false",
    paste("year >=", start_year),
    paste("year <=", end_year),
    email =creds$email,
    pwd =creds$pwd,
    user =creds$user
    )
   })
   names(q_list) <- spp_list
   return(q_list)
}

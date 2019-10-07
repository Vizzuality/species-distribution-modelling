
#' Title
#'
#' @param spp_list List of scientific species names
#' @param bbox sf::bbox object, used to crop extent of occurences
#' @param creds GBIF user credentials
#' @param continent The continent to limit search
#' @param start_year Start year for occurence query
#' @param end_year End year for occurence query
#'
#' @return
#' A named list of sf::sf data.frame(s) per species name.
#' @export
#'
#' @examples
create_spp_occurence_sf_list <- function(spp_list, bbox, creds, path="inst/ext_data/species_gbif", continent="EUROPE", start_year=1986, end_year=2020){

  # Create spp occurence GBIF query list
  q_list <- create_spp_query_list(spp_list, continent=continent, creds, start_year=start_year, end_year=end_year)

  # Create spp occurence GBIF data list, converted to sf::sf
  download_list <- rgbif::occ_download_queue(.list = q_list)
  names(download_list) <- spp_list
  dat_list <- lapply(download_list, function(x){
    dat <- rgbif::occ_download_get(x, path=path)
    df <- rgbif::occ_download_import(dat)[, c("genus", "specificEpithet", "decimalLongitude", "decimalLatitude")]
    return(
      sf::st_crop(
      sf::st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant"),
      bbox)
    )
    })
  names(dat_list)
  return(dat_list)
}

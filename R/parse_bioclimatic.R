#' Parse directory structure of bioclimatic rasters
#'
#' @param path Path to directory with sub-directories of regions each containing historical and future raster files
#'
#' @return A data frame with the directory names, file info and base filename
#' @export
#'
#' @examples
#' df_bc <- parse_bioclimatic(path="inst/ext_data/bioclimatic")
#' summary(df_bc[[1]])
#' head(df_bc[[1]])
parse_bioclimatic <- function(path="inst/ext_data/bioclimatic"){

  # Parse future directory
  parse_future_dir <- function(path="inst/ext_data/bioclimatic/biovar_Sweden/biovar_future"){
    options(stringsAsFactors=F)
    fl <- list.files(path, full.names = T, recursive = T)
    df <- as.data.frame(t(as.data.frame(sapply(fl, strsplit, split="/"))))
    row.names(df) <- NULL
    bn <- unlist(strsplit(df[,ncol(df)], split=".tif"))
    row.names(bn) <- NULL
    #head(bn)
    bn_df <- as.data.frame(t(as.data.frame(sapply(bn, strsplit, split="_"))))
    row.names(bn_df) <- NULL
    names(bn_df) <- c("bio_variable", "model", "scenario", "region", "time_interval")
    #head(bn_df)
    df <- df[,3:ncol(df)]
    names(df) <- c("data_type", "biovar_region", "time_period", "time_interval", "file_name")
    out <- cbind(
      df[,c("data_type", "biovar_region", "time_period")],
      bn_df[, c("bio_variable", "model", "scenario", "region", "time_interval")]
    )
    out$file_name <- fl
    row.names(out) <- NULL
    #head(out)
    return(out)
  }

  # Parse historical directory
  parse_historical_dir <- function(path="inst/ext_data/bioclimatic/biovar_Sweden/biovar_historical"){
    options(stringsAsFactors=F)
    fl <- list.files(path, full.names = T, recursive = T)
    df <- as.data.frame(t(as.data.frame(sapply(fl, strsplit, split="/"))))
    row.names(df) <- NULL
    bn <- unlist(strsplit(df[,ncol(df)], split=".tif"))
    row.names(bn) <- NULL
    #head(bn)
    bn_df <- as.data.frame(t(as.data.frame(sapply(bn, strsplit, split="_"))))
    row.names(bn_df) <- NULL
    names(bn_df) <- c("bio_variable", "model", "region", "time_interval")
    #head(bn_df)
    df <- df[,3:ncol(df)]
    names(df) <- c("data_type", "biovar_region", "time_period", "file_name")
    #head(df)
    out <- cbind(
      df[,c("data_type", "biovar_region", "time_period")],
      bn_df[, c("bio_variable", "model", "region", "time_interval")]
    )
    out$file_name <- fl
    row.names(out) <- NULL
    #head(out)
    return(out)
  }

  # Get list of region dirs
  dl <- list.files(path)

  # Rerurn df
  out <- lapply(dl, function(x){
    return(
      merge(
        parse_historical_dir(path=file.path(path, x, "biovar_historical")),
        parse_future_dir(path=file.path(path, x, "biovar_future")),
        all=T
        )
    )
  })
  names(out) <- dl
  return(out)
}

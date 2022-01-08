#' Get distances between survey points
#'
#' @param survey_df A dataframe of sampling locations with Lon and Lat. Should already be ordered by date and hauljoin, so that the sites are in chronological order.
#' @details If the dataframe is from a historical dataset, Id has to be added and lat/lon converted to title case. If it's from an optimized survey design, Id is the unique identifier for that location.
#' @return a dataframe with the distances (in km) between each pair of points. Essentially a longform version of the distance matrix.
#' @export
#'
#' @examples
#'
library(tidyverse)

get_distances <- function(survey_df = yeardat) {
  # Add Id if the survey data does not contain one (i.e., if it's a historical dataset)
  if (!"Id" %in% names(survey_df)) {
    print("This is a historical dataset; Id column has been added and lat/lon columns have been capitalized")
    
    survey_df <- survey_df %>%
      rename(Lon = "lon", Lat = "lat") %>%
      arrange(DATE, HAULJOIN) %>% # chronological order
      mutate(Id = 1:nrow(survey_df))
  }

  # Add columns for distances
  survey_df <- survey_df %>%
    add_column(
      distance_from_prev = NA,
      cumu_distance = 0
    )

  survey_sf <- sf::st_as_sf(
    x = survey_df,
    coords = c("Lon", "Lat"),
    crs = 4326, agr = "constant"
  )

  distance_matrix_km <- matrix(as.numeric(sf::st_distance(survey_sf) / 1000),
    nrow = nrow(survey_sf)
  )
  rownames(distance_matrix_km) <-
    colnames(distance_matrix_km) <-
    survey_sf$Id
  distance_df <- as.data.frame(distance_matrix_km) %>%
    add_column(surveyId = colnames(distance_matrix_km)) %>%
    pivot_longer(cols = colnames(distance_matrix_km))
  ####

  df_out <- survey_df
  df_out$distance_from_prev[1] <- 0

  for (i in 2:nrow(df_out)) {
    df_out$distance_from_prev[i] <- distance_df %>%
      filter(
        surveyId == df_out$Id[i],
        name == df_out$Id[i - 1]
      ) %>%
      dplyr::select(value) %>%
      as.numeric()
    df_out$cumu_distance[i] <- sum(df_out$distance_from_prev[1:i])
  }

  cumu_dist <- max(df_out$cumu_distance)
  max_dist <- max(df_out$distance_from_prev)
  ###

  return(list(
    distance_df = distance_df,
    cumu_distance = cumu_dist,
    max_distance = max_dist
  ))
}

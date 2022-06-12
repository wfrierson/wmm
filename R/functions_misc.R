#' Derive WMM version based on given time
#'
#' @noRd
#'
#' @param t Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#'
#' @return List of reference year and compatible WMM versions inferred from \code{time}.
.DeriveVersionInfo <- function(t) {
  output <- if(t >= 2025) {
    stop('Time value not supported in current version of wmm package.')
  } else if(t >= 2020 & t < 2025) {
    list(
      'year' = 2020,
      'version' = 'WMM2020'
    )
  } else if(t >= 2015 & t < 2020) {
    list(
      'year' = 2015,
      'version' = c('WMM2015v2', 'WMM2015')
    )
  } else if(t >= 2010) {
    list(
      'year' = 2010,
      'version' = 'WMM2010'
    )
  } else if(t >= 2005){
    list(
      'year' = 2005,
      'version' = 'WMM2005'
    )
  } else {
    stop('Time value not supported in current version of wmm package.')
  }

  return(output)
}

#' Check if given time is consistent with available WMM versions
#'
#' @noRd
#'
#' @param t Annualized date time. E.g., 2015-02-01 = (2015 + 32/365) = 2015.088
#' @param wmmVersion String representing WMM version to use. Must be consistent with \code{time} and one of the following: 'derived', 'WMM2000', 'WMM2005', 'WMM2010', 'WMM2015', 'WMM2015v2'.
.CheckVersionWMM <- function(t, wmmVersion) {
  # Get what WMM versions are compatible with the input time
  derivedVersionInfo <- .DeriveVersionInfo(t)

  if(wmmVersion != 'derived') {
    # Check if the input WMM version is consistent with time-compatible versions
    if(!(wmmVersion %in% derivedVersionInfo[['version']])) {
      stop(
        paste(
          paste(
            'WMM version can only be 1 of the following for the time value',
            'provided:'
          ),
          paste0(shQuote(derivedVersionInfo[['version']]), collapse = ', ')
        )
      )
    }
  }

  return(derivedVersionInfo)
}


#' Check if given horizontal intensity triggers a blackout zone
#'
#' @noRd
#'
#' @param h horizontal intensity, \code{numeric}
#'
#' @return \code{warning} if warranted
#'
.CheckBlackoutZone <- function(h) {
  if (h < 2000) {
    warning(
      'Location is in the blackout zone around the magnetic pole as defined by the WMM military specification (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF89500B.pdf). Compass accuracy is highly degraded in this region.'
    )
  } else if (h < 6000) {
    warning(
      'Location is approaching the blackout zone around the magnetic pole as defined by the WMM military specification (https://www.ngdc.noaa.gov/geomag/WMM/data/MIL-PRF89500B.pdf). Compass accuracy may be degraded in this region.'
    )
  }
}

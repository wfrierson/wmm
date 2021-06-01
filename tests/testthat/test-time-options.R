testthat::context("Testing 'time' alternatives")

testthat::test_that("'time' options are the same", {
  dat <- testData[ h > 6000, ][ 1, ]

  orig <- expect_silent(
    dat[
      , GetMagneticFieldWMM(
        lon = lon,
        lat = lat,
        height = height * 1e3,
        time = floor(year),
        wmmVersion = wmmVersion
      )]
  )

  posix <- expect_silent(
    dat[
      , GetMagneticFieldWMM(
        lon = lon,
        lat = lat,
        height = height * 1e3,
        time = as.POSIXct(paste0(floor(dat$year), "-01-01")),
        wmmVersion = wmmVersion
      )]
  )

  date <- expect_silent(
    dat[
      , GetMagneticFieldWMM(
        lon = lon,
        lat = lat,
        height = height * 1e3,
        time = as.Date(paste0(floor(dat$year), "-01-01")),
        wmmVersion = wmmVersion
      )]
  )

  expect_identical(orig, posix)
  expect_identical(orig, date)
})

library(testthat)
library(data.table)
library(devtools)

devtools::load_all('.')

test_check('wmm')

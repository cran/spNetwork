## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_kfunc.rda",
                           package = "spNetwork", mustWork = TRUE))

## ---- fig.show='hold', fig.align = 'center'-----------------------------------
library(spNetwork)

networkgpkg <- system.file("extdata", "networks.gpkg",
                           package = "spNetwork", mustWork = TRUE)
eventsgpkg <- system.file("extdata", "events.gpkg",
                          package = "spNetwork", mustWork = TRUE)

main_network_mtl <- rgdal::readOGR(networkgpkg,layer="main_network_mtl", verbose = FALSE)
mtl_libraries <- rgdal::readOGR(eventsgpkg,layer="mtl_libraries", verbose = FALSE)
mtl_theatres <- rgdal::readOGR(eventsgpkg,layer="mtl_theatres", verbose = FALSE)

par(mar = c(0, 0, 0, 0))
sp::plot(main_network_mtl)
sp::plot(mtl_libraries, col = "blue", add=T, pch = 20)
sp::plot(mtl_theatres, col = "red", add=T, pch = 20)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  kfun_theatre <- kfunctions(main_network_mtl, mtl_theatres,
#                             start = 0, end = 5000, step = 50,
#                             width = 1000, nsim = 50,
#                             verbose = FALSE, conf_int = 0.05)
#  kfun_theatre$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
kfun_theatre$plotk

## ---- fig.show='hold', fig.align = 'center'-----------------------------------
kfun_theatre$plotg

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  kfun_biblio <- kfunctions(main_network_mtl, mtl_libraries,
#                            start = 0, end = 5000, step = 50,
#                            width = 1000, nsim = 50, verbose = FALSE)
#  
#  kfun_biblio$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
kfun_biblio$plotk

## ---- fig.show='hold', fig.align = 'center'-----------------------------------
kfun_biblio$plotg

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  cross_biblio_theatre <- cross_kfunctions(main_network_mtl, mtl_libraries,
#                          mtl_theatres, start = 0, end = 5000, step = 50,
#                          width = 1000, nsim = 50, verbose = FALSE)
#  
#  cross_biblio_theatre$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
cross_biblio_theatre$plotk

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  cross_theatre_biblio <- cross_kfunctions(main_network_mtl, mtl_theatres,
#                        mtl_libraries, start = 0, end = 5000,
#                        step = 50, width = 1000, nsim = 50, verbose = FALSE)
#  
#  cross_theatre_biblio$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
cross_theatre_biblio$plotk

## ----include = FALSE----------------------------------------------------------
# reset all the user parameters
options(backup_option)
setwd(base_wd)
oldpar <- par(mfrow = c(1,2))
par(oldpar)


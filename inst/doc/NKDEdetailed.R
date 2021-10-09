## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()

## ----message=FALSE, warning=FALSE, fig.height = 3, fig.width  = 3-------------
options(rgl.useNULL=TRUE)
library(spNetwork)
library(rgeos)
library(sp)
library(maptools)
library(rgl)
library(FNN)

# we define a set of lines
wkt_lines <- c(
    "LINESTRING (0 0, 1 0)", "LINESTRING (1 0, 2 0)", "LINESTRING (2 0, 3 0)",
    "LINESTRING (0 1, 1 1)", "LINESTRING (1 1, 2 1)", "LINESTRING (2 1, 3 1)",
    "LINESTRING (0 2, 1 2)", "LINESTRING (1 2, 2 2)", "LINESTRING (2 2, 3 2)",
    "LINESTRING (0 3, 1 3)", "LINESTRING (1 3, 2 3)", "LINESTRING (2 3, 3 3)",
    "LINESTRING (0 0, 0 1)", "LINESTRING (0 1, 0 2)", "LINESTRING (0 2, 0 3)",
    "LINESTRING (1 0, 1 1)", "LINESTRING (1 1, 1 2)", "LINESTRING (1 2, 1 3)",
    "LINESTRING (2 0, 2 1)", "LINESTRING (2 1, 2 2)", "LINESTRING (2 2, 2 3)",
    "LINESTRING (3 0, 3 1)", "LINESTRING (3 1, 3 2)", "LINESTRING (3 2, 3 3)"
    )

# and create a spatial lines dataframe
linesdf <- data.frame(wkt = wkt_lines,
                      id = paste("l",1:length(wkt_lines),sep=""))

geoms <- do.call(rbind,lapply(1:nrow(linesdf),function(i){
  txt <- as.character(linesdf[i,]$wkt)
  geom <- rgeos::readWKT(txt,id=i)
  return(geom)
}))

all_lines <- sp::SpatialLinesDataFrame(geoms, linesdf,match.ID = F)

# we define now one event
event <- data.frame(x=c(1),
                    y=c(1.5))
sp::coordinates(event) <- cbind(event$x,event$y)

# and map the situation
par(mar=c(0.1,0.1,0.1,0.1))
sp::plot(all_lines)
sp::plot(event,add=T,col="red",pch = 19)

## ----message=FALSE, warning=FALSE, echo = FALSE-------------------------------
d3_plot_situation <- function(lines, events, pt_samples, densities, scales){
  
  open3d(scale=scales)
  
  densities < round(densities,7)
  
  ## finding for each event its closest samples
  XY1 <- coordinates(pt_samples)
  XY2 <- coordinates(events)
  idx <- knnx.index(XY1, query = XY2, k = 1)
  
  events$dens <- densities[idx]
  events$x <- XY2[,1]
  events$y <- XY2[,2]
  eidx <- do.call(c,lapply(1:nrow(events), function(i){c(i,i)}))
  vert_lines <- events@data[eidx,]
  vert_lines$dens <- ifelse(1:nrow(vert_lines)%%2 == 0, vert_lines$dens,0)
  
  
  ## plotting the situation
  line_coords <- do.call(rbind,unlist(coordinates(lines), recursive = FALSE))
  sample_coords <- coordinates(pt_samples)
  
  segments3d(x=line_coords[,1],
             y=line_coords[,2],
             z=rep(0, nrow(line_coords)))
  
  segments3d(x=vert_lines$x,
             y=vert_lines$y,
             z=vert_lines$dens,
             col = "red")
  
  coords_events <- coordinates(events)
  
  points3d(
    x = coords_events[,1],
    y = coords_events[,2],
    z = rep(0,nrow(event)),
    col = "red",
    size = 5
  )
  
  points3d(
    x = sample_coords[,1],
    y = sample_coords[,2],
    z = densities,
    col = "blue"
  )
  
  axes3d()
  title3d(xlab="X",ylab="Y",zlab="Z")
}

knitr::knit_hooks$set(webgl = hook_webgl)

## ----message=FALSE, warning=FALSE---------------------------------------------
samples_pts <- lines_points_along(all_lines,0.01)

simple_kernel <- nkde(all_lines, event, w = 1,
                  samples = samples_pts, kernel_name = "quartic", 
                  bw = 2, method = "simple", div = "bw", 
                  digits = 3, tol = 0.001, grid_shape = c(1,1),
                  check = FALSE,
                  verbose = FALSE)

## ----test-rgl, message=FALSE, warning=FALSE, webgl = TRUE, out.width= "75%"----
d3_plot_situation(all_lines, event, samples_pts, simple_kernel, scales = c(1,1,3))

## ----message=FALSE, warning=FALSE---------------------------------------------
discontinuous_kernel <- nkde(all_lines, event, w = 1,
                  samples = samples_pts, kernel_name = "quartic",
                  bw = 2, method = "discontinuous", div = "bw",
                  digits = 3, tol = 0.001, grid_shape = c(1,1),
                  check = FALSE,
                  verbose = FALSE)

## ----test-rgl2, message=FALSE, warning=FALSE, webgl = TRUE, out.width= "75%"----
d3_plot_situation(all_lines, event, samples_pts, discontinuous_kernel, scales = c(1,1,3))

## ----message=FALSE, warning=FALSE---------------------------------------------
continuous_kernel <- nkde(all_lines, event, w = 1,
                  samples = samples_pts, kernel_name = "quartic",
                  bw = 2, method = "continuous", div = "bw",
                  digits = 3, tol = 0.001, grid_shape = c(1,1),
                  check = FALSE,
                  verbose = FALSE)

## ----test-rgl3, message=FALSE, warning=FALSE, webgl = TRUE, out.width= "75%"----
d3_plot_situation(all_lines, event, samples_pts, continuous_kernel, scales = c(1,1,3))

## ----message=FALSE, warning=FALSE---------------------------------------------

# we define now two events
events <- data.frame(x=c(1,2),
                    y=c(1.5,1.5))
sp::coordinates(events) <- cbind(events$x,events$y)

continuous_kernel <- nkde(all_lines, events, w = 1,
                  samples = samples_pts, kernel_name = "quartic",
                  bw = 1.3, method = "continuous", div = "bw",
                  digits = 3, tol = 0.001, grid_shape = c(1,1),
                  check = FALSE,
                  verbose = FALSE)

## ----test-rgl4, message=FALSE, warning=FALSE, webgl = TRUE, out.width= "75%"----
d3_plot_situation(all_lines, events, samples_pts, continuous_kernel, scales = c(1,1,3))

## ----test-rgl5, message=FALSE, warning=FALSE, webgl = TRUE, out.width= "75%"----

continuous_kernel <- nkde(all_lines, events, w = 1,
                  samples = samples_pts, kernel_name = "quartic",
                  bw = 1.6, method = "continuous", div = "bw",
                  digits = 3, tol = 0.001, grid_shape = c(1,1),
                  check = FALSE,
                  verbose = FALSE)

d3_plot_situation(all_lines, events, samples_pts, continuous_kernel, scales = c(1,1,3))

## ----include = FALSE----------------------------------------------------------
# reset all the user parameters
options(backup_option)
setwd(base_wd)
oldpar <- par(mfrow = c(1,2))
par(oldpar)


## ----message=FALSE, warning=FALSE---------------------------------------------
# first load data and packages
library(sp)
library(maptools)
library(rgeos)
library(spNetwork)
library(raster)
library(tmap)
library(FNN)

networkgpkg <- system.file("extdata", "networks.gpkg",
                           package = "spNetwork", mustWork = TRUE)
eventsgpkg <- system.file("extdata", "events.gpkg",
                          package = "spNetwork", mustWork = TRUE)

mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network",verbose = FALSE)
bike_accidents <- rgdal::readOGR(eventsgpkg,layer="bike_accidents", verbose = FALSE)

# then plotting the data
plot(mtl_network)
plot(bike_accidents,add=T,col='red',pch = 19)

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_network_build.rda",
                           package = "spNetwork", mustWork = TRUE))

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # calculating the density values
#  densities <- nkde(mtl_network,
#                    events = bike_accidents,
#                    w = rep(1,nrow(bike_accidents)),
#                    samples = bike_accidents,
#                    kernel_name = "quartic",
#                    bw = 300, div= "bw",
#                    method = "discontinuous", digits = 2, tol = 0.5,
#                    grid_shape = c(1,1), max_depth = 8,
#                    agg = 5,
#                    sparse = TRUE,
#                    verbose = FALSE)
#  
#  bike_accidents$density <- densities * 1000

## ----message=FALSE, warning=FALSE---------------------------------------------
# mapping the density values
tm_shape(mtl_network) + 
  tm_lines(col = "black") +
  tm_shape(bike_accidents) + 
  tm_dots(col = "density", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE) 

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  bike_accidents$weight <- 1
#  agg_points <- aggregate_points(bike_accidents, maxdist = 5)
#  
#  agg_points$OID <- 1:nrow(agg_points)
#  mtl_network$LineID <- 1:nrow(mtl_network)
#  
#  snapped_accidents <- snapPointsToLines2(agg_points,
#                                          mtl_network,
#                                          "LineID")

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  new_lines <- split_lines_at_vertex(mtl_network,
#                                     snapped_accidents,
#                                     snapped_accidents$nearest_line_id,
#                                     mindist = 0.1)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  new_lines$OID <- 1:nrow(new_lines)
#  new_lines$length <- gLength(new_lines, byid = TRUE)
#  
#  graph_result <- build_graph(new_lines, 2, "length", attrs = TRUE)

## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  btws <- igraph::betweenness(graph_result$graph, directed = FALSE,
#                              normalized = TRUE)
#  vertices <- graph_result$spvertices
#  vertices$btws <- btws

## ----message=FALSE, warning=FALSE---------------------------------------------
# mapping the betweenness
tm_shape(vertices) + 
  tm_dots(col = "btws", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE) 


## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  # first: nn merging between snapped points and nodes
#  xy1 <- coordinates(snapped_accidents)
#  xy2 <- coordinates(vertices)
#  corr_nodes <- get.knnx(xy2, xy1, k=1)
#  
#  snapped_accidents$btws <- vertices$btws[corr_nodes[[1]]]
#  
#  # second: nn merging between original points and snapped points
#  xy1 <- coordinates(bike_accidents)
#  xy2 <- coordinates(snapped_accidents)
#  
#  corr_nodes <- get.knnx(xy2, xy1, k=1)
#  bike_accidents$btws <- snapped_accidents$btws[corr_nodes[[1]]]

## ----message=FALSE, warning=FALSE---------------------------------------------
# mapping the results
tm_shape(bike_accidents) + 
  tm_dots(col = "btws", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE)

tm_shape(bike_accidents) + 
  tm_dots(col = "density", style = "kmeans",
          n = 6, size = 0.1, palette = "-RdYlBu")+
  tm_layout(legend.outside = TRUE) 

## ----message=FALSE, warning=FALSE---------------------------------------------
cor.test(bike_accidents$density, bike_accidents$btws)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### spatial indexing ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Build a quadtree
#'
#' @description Generate a quadtree object from package SearchTrees, useful to speed up
#' spatial requesting (INTERNAL).
#'
#' @param data a feature collection of linestrings or a feature collection of points
#' @return quadtree object from package SearchTrees
#' @export
#' @examples
#' #This is an internal function, no example provided
build_quadtree <- function(data){
  #step1 : extracting the bbox of the geometrie
  if(class(data)[[1]] != "sf"){
    stop("The data argument must be a feature collection from the package sf")
  }
  geom_type <- unique(st_geometry_type(data))

  if(geom_type == "LINESTRING"){

    coords <- st_coordinates(data)
    ids <- coords[,3]
    coords_list <- split(coords[,1:2], f = ids)

    bbox_coords <- lapply(coords_list,function(i){
      line_coords <- matrix(i, ncol = 2, byrow = FALSE)
      maxs <- apply(line_coords,MARGIN=2,FUN = max)
      mins <- apply(line_coords,MARGIN=2,FUN = min)
      row <- c(maxs,mins)
      return(row)
    })

    bbox_coords <- do.call(rbind,bbox_coords)
    #step2 : generate the spatial index
    spIndex <- SearchTrees::createTree(bbox_coords,dataType = "rect")

  }else if (geom_type == "POINT"){
    coords <- st_coordinates(data)
    spIndex <- SearchTrees::createTree(coords,dataType = "point")
  }else {
    stop("The supported geometry types are POINT and LINESTRING")
  }

  return(spIndex)
}


#' @title Spatial request
#'
#' @description Use a quadtree index to perform spatial request.
#'
#' @param geometry sf like object (feature collection or simple geometry)
#' @param tree a tree object from package SearchTrees
#' @param data the original data used to build the tree object
#' @return a subset of data, intersecting geometry
#' @importFrom sf st_intersects
#' @export
#' @examples
#' #This is an internal function, no example provided
spatial_request <- function(geometry,tree,data){
  ## step1 : find candidates
  #box <- t(raster::bbox(geometry))
  box <- matrix(st_bbox(geometry), nrow = 2, byrow = TRUE)
  idx <- SearchTrees::rectLookup(tree,box[1,],box[2,])
  candidates <- data[idx,]
  if(nrow(candidates) > 0){
    ## step2 : find real intersection
    final_vector <- st_intersects(candidates, geometry, sparse = FALSE)[,1]
    final_data <- subset(candidates,final_vector)
    return(final_data)
  }else{
    return(candidates)
  }
}


#' @title Find closest points
#'
#' @description Solve the nearest neighbour problem for two feature collections of points
#' This is a simple wrap-up of the dbscan::kNN function
#'
#' @param origins a feature collection of points
#' @param targets a feature collection of points
#' @return for each origin point, the index of the nearest target point
#' @export
#' @examples
#' #This is an internal function, no example provided
#' eventsgpkg <- system.file("extdata", "events.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_libraries <- sf::st_read(eventsgpkg,layer="mtl_libraries")
#' mtl_theatres <- sf::st_read(eventsgpkg,layer="mtl_theatres")
#' close_libs <- closest_points(mtl_theatres, mtl_libraries)
closest_points <- function(origins, targets){

  xy_origins <- st_coordinates(origins)
  xy_targets <- st_coordinates(targets)

  # idx <- FNN::knnx.index(data = xy_targets,
  #                       query = xy_origins,
  #                       k = 1)
  if(nrow(xy_targets) > 1){
    idx <- dbscan::kNN(x = xy_targets, query = xy_origins, k = 1)$id
  }else if(nrow(xy_targets) == 1){
    return(rep(1, nrow(xy_origins)))
  }else{
    stop("Error in the function closest_points, less than one target...")
  }

  return(as.vector(idx[,1]))
}




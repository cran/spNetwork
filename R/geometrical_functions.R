
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### helper functions ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#defining some global variables to prevent check error (weird flex but ok)
utils::globalVariables(c("Xs", "Ys"))

#' @title Coordinates to unique character vector
#'
#' @description Generate a character vector based on a coordinates matrix and
#' the maximum number of digits to keep.
#'
#' @param coords A n * 2 matrix representing the coordinates
#' @param digits The number of digits to keep from the coordinates
#' @return A vector character vector of length n
#' @importFrom data.table data.table tstrsplit :=
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
sp_char_index <- function(coords, digits) {
    tempdf <- data.frame(Xs = as.character(coords[, 1]),
                         Ys = as.character(coords[, 2]))
    tempdt <- data.table(tempdf)

    tempdt[,  c("xint", "xdec") := tstrsplit(Xs, ".", fixed = TRUE)]
    tempdt[,  c("yint", "ydec") := tstrsplit(Ys, ".", fixed = TRUE)]

    X <- paste(tempdt$xint, substr(tempdt$xdec, start = 1, stop = digits),
        sep = ".")
    X <- gsub("NA","0",X,fixed = TRUE)
    Y <- paste(tempdt$yint, substr(tempdt$ydec, start = 1, stop = digits),
        sep = ".")
    Y <- gsub("NA","0",Y,fixed = TRUE)
    return(paste(X, Y, sep = "_"))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### lines manipulations ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Get lines extremities
#'
#' @description Generate a SpatialPointsDataFrame with the first and last vertex
#'   of each line in a SpatialLinesDataFrame.
#'
#' @param lines A SpatialLinesDataFrame
#' @return A SpatialPointsDataFrame
#' @importFrom sp coordinates
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
lines_extremities <- function(lines) {
    coords <- coordinates(lines)
    ptcoords <- lapply(coords, function(line) {
        linecoords <- line[[1]]
        ptscoords <- linecoords[c(1, nrow(linecoords)), ]
        return(ptscoords)
    })
    newpts <- do.call(rbind, ptcoords)
    ids <- lapply(1:nrow(lines), function(i) {
        return(c(i, i))
    })
    types <- lapply(1:nrow(lines), function(i) {
        return(c("start", "end"))
    })
    ids <- do.call("c", ids)
    types <- do.call("c", types)
    data <- lines@data[ids, ]
    data$X <- newpts[, 1]
    data$Y <- newpts[, 2]
    data$pttype <- types
    sp::coordinates(data) <- cbind(data$X, data$Y)
    raster::crs(data) <- raster::crs(lines)
    return(data)
}

#' @title Remove loops
#'
#' @description Remove from a SpatialLinesDataFrame the lines that have the
#' same starting and ending point.
#'
#' @param lines A SpatialLinesDataFrame
#' @param digits An integer indicating the number of digits to keep for the
#' spatial coordinates
#' @return A SpatialLinesDataFrame
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
remove_loop_lines <- function(lines, digits){
    lines_ext <- lines_extremities(lines)
    starts <- subset(lines_ext,lines_ext$pttype=="start")
    ends <- subset(lines_ext,lines_ext$pttype=="end")
    starts$spoid <- sp_char_index(coordinates(starts),digits)
    ends$spoid <- sp_char_index(coordinates(ends),digits)
    test <- starts$spoid != ends$spoid
    return(subset(lines,test))
}

#' @title Unify lines direction
#'
#' @description A function to deal with the directions of lines.
#'
#' @param lines A SpatialLinesDataFrame
#' @param field Indicate a field giving informations about authorized
#' traveling direction on lines. if NULL, then all lines can be used in both
#' directions. Must be the name of a column otherwise. The values of the
#' column must be "FT" (From - To), "TF" (To - From) or "Both".
#' @return A SpatialLinesDataFrame
#' @importFrom sp coordinates Line Lines SpatialLines SpatialLinesDataFrame
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
lines_direction <- function(lines,field){
    listlines <- lapply(1:nrow(lines),function(i){
        line <- lines[i,]
        if(line[[field]]=="TF"){
            coords <- coordinates(line)[[1]][[1]]
            new_line <- Lines(list(Line(coords[nrow(coords):1,])), ID = i)
        }else{
            coords <- coordinates(line)[[1]][[1]]
            new_line <- Lines(list(Line(coords)), ID = i)
        }
        return(new_line)
    })
    spLines <- SpatialLines(listlines)
    df <- SpatialLinesDataFrame(spLines,lines@data,match.ID = FALSE)
    raster::crs(df) <- raster::crs(lines)
    return(df)
}


#' @title Add vertices to a single line
#'
#' @description Add vertices (SpatialPoints) to a single line (SpatialLines),
#'   may fail if the lines geometries are self intersecting.
#'
#' @param line The SpatialLine to modify
#' @param points The SpatialPoints to add to as vertex to the lines
#' @param i The index of the line (for recombining after all lines)
#' @param mindist The minimum distance between one point and the extremity of
#'   the line to add the point as a vertex.
#' @return A matrix of coordinates to build the new line
#' @importFrom rgeos gProject
#' @importFrom sp coordinates SpatialPoints Line Lines
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
add_vertices <- function(line, points, i, mindist) {
    # extract coordinates
    line_coords <- coordinates(line)[[1]][[1]]
    original_distances <- sapply(1:nrow(line_coords),function(i){
        if (i==0){
            return(0)
        }else{
            return(sqrt(sum((line_coords[i,]-line_coords[i-1,])**2)))
        }
    })
    line_coords <- cbind(line_coords, cumsum(original_distances))
    tot_lengths <- max(line_coords[,3])
    # calculate lengths
    pt_coords <- coordinates(points)
    lengths <- gProject(line, points)
    pt_coords <- cbind(pt_coords, lengths)
    pt_coords <- pt_coords[(lengths>mindist & lengths<(tot_lengths-mindist)),]
    all_coords <- rbind(line_coords,pt_coords)
    # reorder the coordinate matrix
    ord_coords <- all_coords[order(all_coords[,3]), ]
    return(ord_coords[,1:2])
}



#' @title Add vertices to a SpatialLinesDataFrame
#'
#' @description Add vertices (SpatialPoints) to their nearest lines
#'   (SpatialLines), may fail if the lines geometries are self intersecting.
#'
#' @param lines The SpatialLinesDataframe to modify
#' @param points The SpatialPoints to add to as vertex to the lines
#' @param nearest_lines_idx For each point, the index of the nearest line
#' @param mindist The minimum distance between one point and the extremity of
#'   the line to add the point as a vertex.
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @importFrom sp coordinates SpatialPoints SpatialLinesDataFrame Line
#'   SpatialLines
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
add_vertices_lines <- function(lines, points, nearest_lines_idx, mindist) {
    #pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    new_lines_list <- lapply(1:nrow(lines), function(i) {
        #setTxtProgressBar(pb, i)
        line <- lines[i, ]
        testpts <- nearest_lines_idx == i
        if (any(testpts)) {
            okpts <- subset(points,testpts)
            newline <- add_vertices(line, okpts, i, mindist)
            return(newline)
        } else {
            sline <- coordinates(line)[[1]][[1]]
            return(sline)
        }

    })
    final_lines <- do.call(raster::spLines,new_lines_list)
    final_lines <- SpatialLinesDataFrame(final_lines,
                                         lines@data,match.ID = FALSE)
    raster::crs(final_lines) <- raster::crs(lines)
    return(final_lines)
}


#' @title Cut lines into lixels
#'
#' @description Cut a SpatialLines object into lixels with a specified minimal
#'   distance may fail if the lines geometries are self intersecting.
#'
#' @param lines The SpatialLinesDataframe to modify
#' @param lx_length The length of a lixel
#' @param mindist The minimum length of a lixel. After cut, if the length of the
#'   final lixel is shorter than the minimum distance, then it is added to the
#'   previous lixel. if NULL, then mindist = maxdist/10. Note that the segments
#'   that are already shorter than the minimum distance are not modified.
#' @param verbose A Boolean indicating if a progress bar should be displayed
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom  sp coordinates Line Lines SpatialLines SpatialLinesDataFrame
#'   SpatialPoints
#' @importFrom rgeos gLength gInterpolate
#' @export
#' @examples
#' \dontrun{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' lixels <- lixelize_lines(mtl_network,150,50)
#' }
lixelize_lines <- function(lines, lx_length, mindist = NULL, verbose = FALSE) {
    if (is.null(mindist)) {
        mindist <- lx_length/10
    }
    if(verbose){
        pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    }
    newlixels <- lapply(1:nrow(lines), function(i) {
        if(verbose){
            setTxtProgressBar(pb, i)
        }
        line <- lines[i, ]
        tot_length <- gLength(line)
        if (tot_length < lx_length+mindist) {
            coords <- coordinates(line)
            return(list(coords[[1]][[1]]))
        } else {
            # producing the points to snapp on
            distances <- seq(lx_length, tot_length, lx_length)
            if ((tot_length - distances[[length(distances)]]) < mindist) {
                distances <- distances[1:(length(distances) - 1)]
            }
            points <- t(sapply(distances, function(d) {
                return(coordinates(gInterpolate(line, d)))
            }))
            points <- data.frame(x = points[, 1], y = points[, 2],
                                 distance = distances,
                                 type = "cut")
            # extracting the original coordinates
            coords <- SpatialPoints(coordinates(line))
            xy <- coordinates(coords)
            points2 <- data.frame(x = xy[, 1],
                                  y = xy[, 2],
                                  distance = gProject(line,coords),
                                  type = "base")
            if(points2$distance[nrow(points2)]==0){
                points2$distance[nrow(points2)] <- tot_length
            }
            # merging both and sorting
            allpts <- rbind(points, points2)
            allpts <- allpts[order(allpts$distance), ]
            # and now splitting this motherfucker
            indices <- c(0, which(allpts$type == "cut"), nrow(allpts))
            lixels <- lapply(1:(length(indices) - 1), function(j) {
                pts <- allpts[indices[[j]]:indices[[j + 1]], ]
                return(as.matrix(pts[, 1:2]))
            })
            return(lixels)
        }
    })
    oids <- lapply(1:nrow(lines), function(i) {
        return(rep(i, length(newlixels[[i]])))
    })
    oids <- do.call("c", oids)

    new_lines <- do.call(raster::spLines,unlist(newlixels,recursive = FALSE))
    df <- lines@data[oids, ]
    if(class(df) != "dataframe"){
        df <- data.frame("oid" = df)
    }
    new_splines <- SpatialLinesDataFrame(new_lines, df, match.ID = FALSE)
    raster::crs(new_splines) <- raster::crs(lines)
    return(new_splines)

}


#'@title Cut lines into lixels (multicore)
#'
#'@description Cut a SpatialLines object into lixels with a specified minimal
#'  distance may fail if the lines geometries are self intersecting with
#'  multicore support.
#'
#'@param lines The SpatialLinesDataframe to modify
#'@param lx_length The length of a lixel
#'@param mindist The minimum length of a lixel. After cut, if the length of the
#'  final lixel is shorter than the minimum distance, then it is added to the
#'  previous lixel. If NULL, then mindist = maxdist/10
#'@param verbose A Boolean indicating if a progress bar must be displayed
#'@param chunk_size The size of a chunk used for multiprocessing. Default is
#'  100.
#'@return An object of class SpatialLinesDataFrame (package sp)
#'@export
#'@importFrom utils capture.output
#' @examples
#' \dontrun{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' future::plan(future::multisession(workers=2))
#' lixels <- lixelize_lines.mc(mtl_network,150,50)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")){
#' future::plan(future::sequential)
#' }
#'}
lixelize_lines.mc <- function(lines, lx_length, mindist = NULL, verbose = TRUE, chunk_size = 100) {
    chunks <- split(1:nrow(lines), rep(1:ceiling(nrow(lines) / chunk_size),
                each = chunk_size, length.out = nrow(lines)))
    chunks <- lapply(chunks,function(x){return(lines[x,])})
    # step2 : starting the function
    iseq <- 1:length(chunks)
    if (verbose) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                chunk_lines <- chunks[[i]]
                invisible(capture.output(new_lines <- lixelize_lines(chunk_lines,
                    lx_length, mindist)))
                return(new_lines)
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, function(i) {
            chunk_lines <- chunks[[i]]
            invisible(capture.output(new_lines <- lixelize_lines(chunk_lines,
                lx_length, mindist)))
            return(new_lines)
        })
    }

    final_lines <- do.call("rbind", values)
    return(final_lines)
}


#' @title LineString to simple Line
#'
#' @description Split the polylines of a SpatialLinesDataFrame object in simple
#'   lines.
#'
#' @param lines The SpatialLinesDataframe to modify
#' @return An object of class SpatialLinesDataFrame (package sp)
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom sp coordinates SpatialLinesDataFrame SpatialLines Lines Line
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
simple_lines <- function(lines) {
    ## extracting the coordinates of the lines
    allcoords <- coordinates(lines)
    counts <- sapply(1:length(allcoords), function(i) {
        return(nrow(allcoords[[i]][[1]]) - 1)
    })
    oids <- sapply(1:nrow(lines), function(i) {
        return(rep(i, counts[[i]]))
    })
    if(nrow(lines)>1){
        oids <- do.call("c", oids)
    }

    ## using the coordinates to create newlines
    new_lines <- lapply(1:length(allcoords), function(i) {
        coords <- allcoords[[i]][[1]]
        segment_lines <- lapply(1:(nrow(coords) - 1), function(i) {
            mat <- coords[i:(i + 1), ]
            return(mat)
        })
        return(segment_lines)
    })
    data <- lines@data[oids, ]

    final_lines <- do.call(raster::spLines,unlist(new_lines,recursive = FALSE))
    final_lines <- SpatialLinesDataFrame(final_lines,data, match.ID = FALSE)
    raster::crs(final_lines) <- raster::crs(lines)
    return(final_lines)
}

#' @title Center points of lines
#'
#' @description Generate a SpatialPointsDataFrame with line center points. The
#'   points are located at center of the line based on the length of the line.
#'
#' @param lines The SpatialLinesDataframe to use
#' @return An object of class SpatialPointsDataFrame (package sp)
#' @importFrom sp coordinates SpatialPointsDataFrame SpatialPoints
#' @importFrom rgeos gInterpolate
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' \dontrun{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' centers <- lines_center(mtl_network)
#' }
lines_center <- function(lines) {
    centerpoints <- maptools::SpatialLinesMidPoints(lines)
    return(centerpoints)
}

#' @title Add center vertex to lines
#'
#' @description Add to each line of a SpatialLinesDataFrame an additional vertex at its
#' center.
#'
#' @param lines The SpatialLinesDataframe to use
#' @return An object of class SpatialLinesDataframe (package sp)
#' @importFrom sp coordinates SpatialPointsDataFrame SpatialPoints
#' @importFrom rgeos gInterpolate
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
add_center_lines <- function(lines) {
    pb <- txtProgressBar(min = 0, max = nrow(lines), style = 3)
    listline <- sapply(1:nrow(lines), function(i) {
        setTxtProgressBar(pb, i)
        line <- lines[i, ]
        pt <- gInterpolate(line, gLength(line) / 2)
        newline <- add_vertices(line, pt, i)
        return(newline)
    })
    final_lines <- SpatialLinesDataFrame(SpatialLines(listline), lines@data,
        match.ID = FALSE)
    raster::crs(final_lines) <- raster::crs(lines)
    return(final_lines)
}


#'@title Add center vertex to lines (multicore)
#'
#'@description Add to each line of a SpatialLinesDataFrame an additional vertex
#'  at its center with multicore support.
#'
#'@param lines The SpatialLinesDataframe to use
#'@param show_progress A Boolean indicating if a progress bar must be displayed
#'@param chunk_size The size of a chunk used for multiprocessing. Default is
#'  100.
#'@return An object of class SpatialLinesDataframe (package sp)
#'@importFrom utils capture.output
#'@keywords internal
#' @examples
#' #This is an internal function, no example provided
#' \dontshow{
#'    ## R CMD check: make sure any open connections are closed afterward
#'    if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#'}
add_center_lines.mc <- function(lines, show_progress = TRUE, chunk_size = 100) {
    # step1 : splitting the data into chunks
    chunks <- split(1:nrow(lines), rep(1:ceiling(nrow(lines) / chunk_size), each = chunk_size,
        length.out = nrow(lines)))
    chunks <- lapply(chunks, function(x){return(lines[x,])})
    # step2 : starting the function
    iseq <- 1:length(chunks)
    if (show_progress) {
        progressr::with_progress({
            p <- progressr::progressor(along = iseq)
            values <- future.apply::future_lapply(iseq, function(i) {
                p(sprintf("i=%g", i))
                chunk_lines <- chunks[[i]]
                invisible(capture.output(new_lines <- add_center_lines(chunk_lines)))
                return(new_lines)
            })
        })
    } else {
        values <- future.apply::future_lapply(iseq, function(i) {
            chunk_lines <- chunks[[i]]
            invisible(capture.output(new_lines <- add_center_lines(chunk_lines)))
            return(new_lines)
        })
    }

    final_lines <- do.call("rbind", values)
    return(final_lines)
}


#' @title Points along lines
#'
#' @description Generate points along the lines of a SpatialLinesDataFrame.
#'
#' @param lines The SpatialLinesDataframe to use
#' @param dist The distance between the points along the lines
#' @return An object of class SpatialLinesDataframe (package sp)
#' @importFrom utils capture.output
#' @export
#' @examples
#' \dontrun{
#' networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
#' mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
#' new_pts <- lines_points_along(mtl_network,50)
#' }
lines_points_along <- function(lines,dist){
    lenghts <- gLength(lines, byid = TRUE)
    list_pts <- lapply(1:nrow(lines),function(i){
        line <- lines[i,]
        line_lenght <- lenghts[i]
        distances <- seq(0,line_lenght,dist)
        pts <- gInterpolate(line,distances)
        return(pts)
    })
    oids <- lapply(1:length(list_pts),function(i){rep(i,length(list_pts[[i]]))})
    oids <- do.call("c",oids)
    all_pts <- do.call(rbind,list_pts)
    data <- lines@data[oids,]
    all_pts <- sp::SpatialPointsDataFrame(all_pts,data)
    return(all_pts)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### functions on polygons ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @title Points along polygon boundary
#'
#' @description Generate a SpatialPointsDataFrame by placing points along the border of
#' polygons of a SpatialPolygonsDataFrame.
#'
#' @param polygons A SpatialPolygonsDataFrame
#' @param dist The distance between the points
#' @importFrom rgeos gBoundary
#' @return A SpatialPolygonsDataFrame representing the grid
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
surrounding_points <- function(polygons,dist){
    #extracting the boundaries and their lengths
    boundaries <- gBoundary(polygons, byid = TRUE)
    df <- sp::SpatialLinesDataFrame(boundaries,polygons@data)
    all_pts <- lines_points_along(df,dist)
    return(df)
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### gridding function ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Spatial grid
#'
#' @description Generate a grid of a specified shape in the bbox of a Spatial object.
#'
#' @param grid_shape A numeric vector of length 2 indicating the number
#' of rows and the numbers of columns of the grid
#' @param spatial A list of SpatialDataFrames objects (package sp)
#' @return A SpatialPolygonsDataFrame representing the grid
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
build_grid <- function(grid_shape, spatial) {
    boxes <- lapply(spatial,sp::bbox)
    boxes <- do.call(cbind,boxes)
    v1 <- as.numeric(c(min(boxes[1,]),max(boxes[1,])))
    v2 <- as.numeric(c(min(boxes[2,]),max(boxes[2,])))
    box <- rbind(v1,v2)
    if(prod(grid_shape)==1){
        #ext <- raster::extent(spatial)
        ext <- raster::extent(box[1,1],box[1,2],box[2,1],box[2,2])
        poly <- methods::as(ext, "SpatialPolygons")
        raster::crs(poly) <- raster::crs(spatial)
        return(poly)
    }else{
        ## step1 : creating the grid
        x <- seq(from = box[1, 1], to = box[1, 2], length.out = grid_shape[[1]])
        y <- seq(from = box[2, 1], to = box[2, 2], length.out = grid_shape[[2]])
        xy <- expand.grid(x = x, y = y)
        grid.pts <- sp::SpatialPointsDataFrame(coords = xy, data = xy)
        raster::crs(grid.pts) <- raster::crs(spatial)
        sp::gridded(grid.pts) <- TRUE
        grid <- methods::as(grid.pts, "SpatialPolygons")
        return(grid)
    }
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### snapping function ####
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Nearest point on segment
#'
#' @description Find the nearest projected point on a segment (from maptools)
#'
#' @param s The coordinates of the segment
#' @param p The coordinates of the point
#'
#' @return A numeric vector with the coordinates of the projected point
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
nearestPointOnSegment <- function(s, p){
    # Adapted from http://pastebin.com/n9rUuGRh
    ap <- c(p[1] - s[1,1], p[2] - s[1,2])
    ab <- c(s[2,1] - s[1,1], s[2,2] - s[1,2])
    t <- sum(ap*ab) / sum(ab*ab)
    t <- ifelse(t<0,0,ifelse(t>1,1,t))
    # when start and end of segment are identical t is NA
    t <- ifelse(is.na(t), 0, t)
    x <- s[1,1] + ab[1] * t
    y <- s[1,2] + ab[2] * t
    # Return nearest point and distance
    result <- c(x, y, sqrt((x-p[1])^2 + (y-p[2])^2))
    names(result) <- c("X","Y","distance")
    return(result)
}


#' @title Nearest point on Line
#'
#' @description Find the nearest projected point on a LineString (from maptools)
#'
#' @param coordsLine The coordinates of the line (matrix)
#' @param coordsPoint The coordinates of the point
#'
#' @return A numeric vector with the coordinates of the projected point
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
nearestPointOnLine <- function(coordsLine, coordsPoint){
    nearest_points <- vapply(2:nrow(coordsLine), function(x){
        nearestPointOnSegment(coordsLine[(x-1):x,], coordsPoint)
    }, FUN.VALUE=c(0,0,0))

    # Return coordinates of the nearest point on this line
    return(nearest_points[1:2, which.min(nearest_points[3,])])

}


#' @title Nearest feature
#'
#' @description Find the nearest feature from set X in set Y
#'
#' @param x A SpatialDataFrame
#' @param y A SpatialDataFrame
#'
#' @return A numeric vector with the index of the nearest features
#' @keywords internal
#' @examples
#' #This is an internal function, no example provided
nearest <- function(x,y){
    sf1 <- sf::st_as_sf(x)
    sf2 <- sf::st_as_sf(y)
    sf::st_crs(sf1) <- sf::st_crs(sf2)
    best <- sf::st_nearest_feature(sf1,sf2)
    return(best)
}


#' @title Snap points to lines
#'
#' @description Snap points to their nearest lines (edited from maptools)
#'
#' @param points A SpatialPointsDataFrame
#' @param lines A SpatialLinesDataFrame
#' @param idField The name of the column to use as index for the lines
#'
#' @return A SpatialPointsDataFrame with the projected geometries
#' @keywords internal
#' @importFrom methods slot
#' @examples
#' #This is an internal function, no example provided
snapPointsToLines2 <- function(points, lines ,idField = NA) {

    nearest_line_index <- nearest(points,lines)
    coordsLines <- sp::coordinates(lines)
    coordsPoints <- sp::coordinates(points)


    # Get coordinates of nearest points lying on nearest lines
    mNewCoords <- vapply(1:length(points),function(x){
            return(nearestPointOnLine(coordsLines[[nearest_line_index[x]]][[1]],
                               coordsPoints[x,]))}, FUN.VALUE=c(0,0))

    # Recover lines' Ids (If no id field has been specified, take the sp-lines id)
    if (!is.na(idField)) {
        nearest_line_id <- lines@data[,idField][nearest_line_index]
    }  else {
        nearest_line_id <- sapply(slot(lines, "lines"), function(i){slot(i, "ID")})[nearest_line_index]
    }
    # Create data frame and sp points
    df <- points@data
    df$nearest_line_id <- nearest_line_id

    return(SpatialPointsDataFrame(coords=t(mNewCoords),
                           data=df,
                           proj4string=raster::crs(points)))

}



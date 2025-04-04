#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### NETWORK DIGGLE CORRECTION FACTOR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# this function is not usefull anymore. We keep it only for debuging purpose
# corrfactor_simple <- function(graph,events,edges,bws){
#   tree_edges <- build_quadtree(edges)
#   buffers <- st_buffer(events, dist = bws)
#   #iterons sur chacun des evenements
#   dfs <- lapply(1:nrow(events),function(i){
#     e <- events[i,]
#     y <- e$vertex_id
#     bw <- bws[[i]]
#     ## step1 selecting the edges inside of the radius
#     #buff <- gBuffer(e,width=bw)
#     buff <- buffers[,i]
#     ok_edges <- spatial_request(buff,tree_edges,edges)
#     ## Step3 for each edge, find its two vertices
#     vertices <- ends(graph,ok_edges$edge_id,names = FALSE)
#     ## step4 calculate the the distance between the start node and each edge
#     ## vertex
#     un_vertices <- unique(c(vertices[,1],vertices[,2]))
#     dist1 <- as.numeric(distances(graph,y,to=un_vertices,mode="out"))
#
#     dist_table <- data.frame("vertex"=un_vertices,
#                              "distance" = dist1)
#     ## step5 aggregate all the data
#     df_edges <- data.frame("edge_id" = ok_edges$edge_id,
#                            "weight" = ok_edges$weight,
#                            "node1" = vertices[,1],
#                            "node2" = vertices[,2]
#     )
#     A <- data.table(df_edges)
#     B <- data.table(dist_table)
#     df_edges$d1 <- A[B, on = c("node1" = "vertex"),
#                      names(B) := mget(paste0("i.", names(B)))]$distance
#     df_edges$d2 <- A[B, on = c("node2" = "vertex"),
#                      names(B) := mget(paste0("i.", names(B)))]$distance
#
#     df_edges <- subset(df_edges,df_edges$d1<bw & df_edges$d2<bw)
#
#     df_edges$ecart <- abs(df_edges$d1 - df_edges$d2)
#     df_edges$lower <- with(df_edges, pmin(d1, d2))
#     df_edges$upper <- with(df_edges, pmax(d1, d2))
#     ## creation de la partie 1
#     part1 <- df_edges[c("lower","weight","edge_id")]
#     part1$distances <- part1$lower
#     part1$edge_size <- df_edges$weight - df_edges$ecart
#     part1$lower <- NULL
#     ## creation de la partie 2
#     part2 <- df_edges[c("upper","weight","edge_id")]
#     part2$distances <- part2$upper
#     part2$edge_size <- df_edges$ecart
#     part2$upper <- NULL
#     totdf <- rbind(part1,part2)
#     totdf$alpha <- 1
#     return(subset(totdf,totdf$edge_size>0))
#   })
#   return(dfs)
# }



#' @title Split boundary of polygon
#'
#' @description A function to cut the boundary of the study area into chunks.
#'
#' @param polygon The polygon representing the study area
#' @param bw The maximum bandwidth
#' @keywords internal
#' @importFrom sf st_boundary st_sf
#' @return A feature collection of linestrings
split_border <- function(polygon,bw){
  boundaries <- st_boundary(polygon)
  boundaries <- st_sf(id = 1, geom = boundaries)
  lixels <- lixelize_lines(boundaries,lx_length = bw*4, mindist = bw)
  return(lixels)
}



# correction_factor <- function(study_area, events, lines, method, bws, kernel_name, tol, digits, max_depth, sparse){
#
#   #merging the polygons (just in case)
#   #study_area <- gUnaryUnion(study_area)
#   study_area <- st_union(study_area)
#
#   kernel_func <- select_kernel(kernel_name)
#   events$goid <- 1:nrow(events)
#   events$bws <- bws
#   bw <- max(bws)
#
#   # step 0 calculate the distances between the points and the border
#   boundaries <- st_boundary(study_area)
#   dists <- as.numeric(st_distance(events, boundaries))
#   ok_events <- subset(events, dists < bw)
#   # step 1 create the border elements
#   chunks <- split_border(study_area,bw)
#   chunks$oid <- 1:nrow(chunks)
#   # step 2  associate each events to it closest chunk
#   snapped <- snapPointsToLines2(ok_events,chunks)
#   ok_events$nearest_line_id <- as.numeric(snapped$nearest_line_id)
#   # step 2 select the elements in each chunks
#   tree_lines <- build_quadtree(lines)
#
#   selections <- lapply(1:nrow(chunks),function(i){
#     part <- chunks[i,]
#     sel_events <- subset(ok_events,ok_events$nearest_line_id == part$oid)
#
#     #si on a aucun evenement ici, on passe
#     if(nrow(sel_events)==0){
#       return(NULL)
#     }
#     poly <- st_convex_hull(sel_events)
#     buff <- st_buffer(poly, dist = bw)
#     sel_lines <- spatial_request(buff,tree_lines,lines)
#
#     sel_lines$oid <- 1:nrow(sel_lines)
#     # splitting the lines at the border intersection
#     inter <- st_intersection(st_geometry(sel_lines),st_geometry(boundaries))
#     if(length(inter) == 0){
#       lines2 <- sel_lines
#     }else{
#       inter <- st_sf(ptid = 1:length(inter), geometry = inter)
#       idx <- nearest_lines(inter, sel_lines, snap_dist = bw, max_iter = 10)
#       lines2 <- add_vertices_lines(sel_lines,inter,idx,tol)
#     }
#
#     #and finaly split the lines at the events
#     sel_events <- snapPointsToLines2(sel_events,lines2,idField = "oid")
#     new_lines <- add_vertices_lines(lines2,sel_events,
#                                     sel_events$nearest_line_id,tol)
#
#     new_lines <- simple_lines(new_lines)
#     new_lines$length <- as.numeric(st_length(new_lines))
#     new_lines <- subset(new_lines,new_lines$length>0)
#     new_lines$oid <- 1:nrow(new_lines)
#     new_lines <- new_lines[c("length","oid")]
#     return(list("sel_lines" = new_lines,
#                 "sel_events" = sel_events))
#   })
#   selections <- selections[lengths(selections) != 0]
#   # step 3 iterate over the selections and calculate the corrections factor
#   values <- lapply(selections,function(sel){
#     sel_lines <- sel$sel_lines
#     sel_events <- sel$sel_events
#
#     #building the local graph
#     graph_result <- build_graph(sel_lines, digits = digits,
#                                 line_weight = "length")
#     graph <- graph_result$graph
#     nodes <- graph_result$spvertices
#     edges <- graph_result$spedges
#
#     edges$is_inside <- st_within(edges,study_area, sparse = FALSE)[,1]
#
#     sel_events$vertex_id <- closest_points(sel_events, nodes)
#
#     neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
#
#     neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})
#     # lets obtain the potential values of each line
#     if(method=="continuous"){
#       if(sparse){
#         dfs <- corrfactor_continuous_sparse(neighbour_list,
#                                                           sel_events$vertex_id,
#                                                           graph_result$linelist,
#                                                           bws, max_depth)
#       }else{
#         dfs <- corrfactor_continuous(neighbour_list,
#                                                    sel_events$vertex_id,
#                                                    graph_result$linelist,
#                                                    bws, max_depth)
#       }
#     }
#     if(method=="discontinuous"){
#       if(sparse){
#         dfs <- corrfactor_discontinuous_sparse(neighbour_list,
#                                                              sel_events$vertex_id,
#                                                              graph_result$linelist,
#                                                              bws, max_depth)
#       }else{
#         dfs <- corrfactor_discontinuous(neighbour_list,
#                                                       sel_events$vertex_id,
#                                                       graph_result$linelist,
#                                                       bws, max_depth)
#       }
#     }
#     if(method=="simple"){
#       #dfs <- corrfactor_simple(graph,sel_events,edges,sel_events$bws)
#       dfs <- corrfactor_discontinuous_sparse(neighbour_list,
#                                              sel_events$vertex_id,
#                                              graph_result$linelist,
#                                              bws, max_depth)
#       dfs <- lapply(dfs, function(x){
#         x$alpha <- 1
#         return(x)
#       })
#     }
#
#     # and finaly calculate the correction factor by integrals
#     # more exactly, it calculate the mass of the event OUTSIDE
#     # the study area
#     corrfactor <- sapply(1:length(dfs),function(j){
#       df <- dfs[[j]]
#       bw <- sel_events$bws[[j]]
#       contribs <- sapply(1:nrow(df),function(i){
#         row <- df[i,]
#         start <- row$distances
#         end <- row$edge_size + start
#         if(bw-start<tol){
#           return(0)
#         }else{
#           val <- cubintegrate(kernel_func,lower=start,upper=end,
#                               bw=bw, relTol = 1e-15)
#           return(val$integral * row$alpha)
#         }
#       })
#       df$contribs <- contribs
#       df2 <- merge(df,st_drop_geometry(edges[c("edge_id","is_inside")]),by="edge_id")
#       outside <- subset(df2,df2$is_inside == FALSE)
#       return(sum(outside$contribs))
#     })
#     #corrfactor <- ifelse(corrfactor==0,1,corrfactor)
#     # and we get here the correction factor by getting the inverse of
#     # the mass inside (1-outside)
#     corrfactor <- 1/(1-corrfactor)
#     df <- data.frame("corrfactor" = corrfactor,
#                      "goid" = sel_events$goid)
#     return(df)
#   })
#   comb_values <- do.call(rbind,values)
#   events_df <- data.table(st_drop_geometry(events))
#   B <- data.table(comb_values)
#   events_df[B,  on = c("goid"),
#             names(B) := mget(paste0("i.", names(B)))]
#
#   events_df$corrfactor <- ifelse(is.na(events_df$corrfactor), 1,
#                                  events_df$corrfactor)
#   return(events_df$corrfactor)
# }




#' @title Border correction for NKDE
#'
#' @description Function to calculate the border correction factor.
#'
#' @param study_area A feature collection of polygons or a polygon, the limit
#'   of the study area.
#' @param events A feature collection of points representing the events on the
#'   network.
#' @param lines The lines used to create the network
#' @param method The method to use when calculating the NKDE, must be one of
#'   simple / discontinuous / continuous (see details for more information)
#' @param bws The kernel bandwidth (in meters) for each event
#' @param kernel_name The name of the kernel to use
#' @param tol When adding the events and the sampling points to the network, the
#'   minimum distance between these points and the lines extremities. When
#'   points are closer, they are added at the extremity of the lines.
#' @param digits The number of digits to keep in the spatial coordinates. It
#'   ensures that topology is good when building the network. Default is 3
#' @param max_depth When using the continuous and discontinuous methods, the
#'   calculation time and memory use can go wild  if the network has a lot of
#'   small edges (area with a lot of intersections and a lot of events). To
#'   avoid it, it is possible to set here a maximum depth. Considering that the
#'   kernel is divided at intersections, a value of 8 should yield good
#'   estimates. A larger value can be used without problem for the discontinuous
#'   method. For the continuous method, a larger value will strongly impact
#'   calculation speed.
#' @param sparse A boolean indicating if sparse or regular matrix should be
#'   used by the Rcpp functions. Regular matrices are faster, but require more
#'   memory and could lead to error, in particular with multiprocessing. Sparse
#'   matrices are slower, but require much less memory.
#' @importFrom sf st_union st_boundary st_distance st_convex_hull st_buffer st_length st_intersection st_within st_crs<-
#' @importFrom stats integrate
#' @importFrom cubature cubintegrate
#' @return A numeric vector with the correction factor values for each event
#' @keywords internal
#' @examples
#' #no example provided, this is an internal function
correction_factor <- function(study_area, events, lines, method, bws, kernel_name, tol, digits, max_depth, sparse){

  #merging the polygons (just in case)
  #study_area <- gUnaryUnion(study_area)
  study_area <- st_union(study_area)

  kernel_func <- select_kernel(kernel_name)
  events$goid <- 1:nrow(events)
  events$bws <- bws
  bw <- max(bws)

  # step 0 calculate the distances between the points and the border
  boundaries <- st_boundary(study_area)
  dists <- as.numeric(st_distance(events, boundaries))
  ok_events <- subset(events, dists < bw)

  # step 1 create the border elements
  chunks <- split_border(study_area,bw)
  chunks$oid <- as.character(1:nrow(chunks))

  # step 2  associate each events to it closest chunk
  snapped <- snapPointsToLines2(ok_events,chunks)
  ok_events$nearest_line_id <- as.character(snapped$nearest_line_id)

  inter_events <- split(ok_events,ok_events$nearest_line_id)

  # step 2 select the elements in each chunks

  polys <- do.call(rbind,lapply(inter_events, function(x){
    sf::st_convex_hull(st_union(x))
  }))
  polys <- st_buffer(st_as_sf(st_as_sfc(polys)), bw)
  st_crs(polys) <- st_crs(lines)
  polys$chunk_id <- names(inter_events)


  inter_lines <- sf::st_join(lines, polys)
  inter_lines <- split(inter_lines,inter_lines$chunk_id)

  selections <- lapply(1:nrow(chunks),function(i){

    chunk_id <- chunks$oid[[i]]
    sel_events <- inter_events[[chunk_id]]

    #si on a aucun evenement ici, on passe
    if(is.null(sel_events)){
      return(NULL)
    }
    sel_lines <- inter_lines[[chunk_id]]

    sel_lines$oid <- 1:nrow(sel_lines)
    # splitting the lines at the border intersection
    inter <- st_intersection(st_geometry(sel_lines),st_geometry(boundaries))
    if(length(inter) == 0){
      lines2 <- sel_lines
    }else{
      inter <- st_sf(ptid = 1:length(inter), geometry = inter)
      idx <- nearest_lines(inter, sel_lines, snap_dist = bw, max_iter = 10)
      lines2 <- add_vertices_lines(sel_lines,inter,idx,tol)
    }

    #and finaly split the lines at the events
    sel_events <- snapPointsToLines2(sel_events,lines2,idField = "oid")
    new_lines <- add_vertices_lines(lines2,sel_events,
                                    sel_events$nearest_line_id,tol)

    new_lines <- simple_lines(new_lines)
    new_lines$length <- as.numeric(st_length(new_lines))
    new_lines <- subset(new_lines,new_lines$length>0)
    new_lines$oid <- 1:nrow(new_lines)
    new_lines <- new_lines[c("length","oid")]
    return(list("sel_lines" = new_lines,
                "sel_events" = sel_events))
  })

  selections <- selections[lengths(selections) != 0]
  # step 3 iterate over the selections and calculate the corrections factor
  values <- lapply(selections,function(sel){
    sel_lines <- sel$sel_lines
    sel_events <- sel$sel_events

    #building the local graph
    graph_result <- build_graph(sel_lines, digits = digits,
                                line_weight = "length")
    graph <- graph_result$graph
    nodes <- graph_result$spvertices
    edges <- graph_result$spedges

    edges$is_inside <- st_within(edges,study_area, sparse = FALSE)[,1]

    sel_events$vertex_id <- closest_points(sel_events, nodes)

    neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")

    neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})
    # lets obtain the potential values of each line
    if(method=="continuous"){
      if(sparse){
        dfs <- corrfactor_continuous_sparse(neighbour_list,
                                            sel_events$vertex_id,
                                            graph_result$linelist,
                                            bws, max_depth)
      }else{
        dfs <- corrfactor_continuous(neighbour_list,
                                     sel_events$vertex_id,
                                     graph_result$linelist,
                                     bws, max_depth)
      }
    }
    if(method=="discontinuous"){
      if(sparse){
        dfs <- corrfactor_discontinuous_sparse(neighbour_list,
                                               sel_events$vertex_id,
                                               graph_result$linelist,
                                               bws, max_depth)
      }else{
        dfs <- corrfactor_discontinuous(neighbour_list,
                                        sel_events$vertex_id,
                                        graph_result$linelist,
                                        bws, max_depth)
      }
    }
    if(method=="simple"){
      #dfs <- corrfactor_simple(graph,sel_events,edges,sel_events$bws)
      dfs <- corrfactor_discontinuous_sparse(neighbour_list,
                                             sel_events$vertex_id,
                                             graph_result$linelist,
                                             bws, max_depth)
      dfs <- lapply(dfs, function(x){
        x$alpha <- 1
        return(x)
      })
    }

    # and finaly calculate the correction factor by integrals
    # more exactly, it calculate the mass of the event OUTSIDE
    # the study area
    corrfactor <- sapply(1:length(dfs),function(j){
      df <- dfs[[j]]
      bw <- sel_events$bws[[j]]
      contribs <- sapply(1:nrow(df),function(i){
        row <- df[i,]
        start <- row$distances
        end <- row$edge_size + start
        if(bw-start<tol){
          return(0)
        }else{
          val <- cubintegrate(kernel_func,lower=start,upper=end,
                              bw=bw, relTol = 1e-15)
          return(val$integral * row$alpha)
        }
      })
      df$contribs <- contribs
      df2 <- merge(df,st_drop_geometry(edges[c("edge_id","is_inside")]),by="edge_id")
      outside <- subset(df2,df2$is_inside == FALSE)
      return(sum(outside$contribs))
    })
    #corrfactor <- ifelse(corrfactor==0,1,corrfactor)
    # and we get here the correction factor by getting the inverse of
    # the mass inside (1-outside)
    corrfactor <- 1/(1-corrfactor)
    df <- data.frame("corrfactor" = corrfactor,
                     "goid" = sel_events$goid)
    return(df)
  })
  comb_values <- do.call(rbind,values)
  events_df <- data.table(st_drop_geometry(events))
  B <- data.table(comb_values)
  events_df[B,  on = c("goid"),
            names(B) := mget(paste0("i.", names(B)))]

  events_df$corrfactor <- ifelse(is.na(events_df$corrfactor), 1,
                                 events_df$corrfactor)
  return(events_df$corrfactor)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### TIME DIGGLE CORRECTION FACTOR ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Time extent correction for NKDE
#'
#' @description Function to calculate the time extent correction factor in tnkde.
#'
#' @param events_time A numeric vector representing when the events occurred
#' @param samples_time A numeric vector representing when the densities will
#' be sampled
#' @param bws_time A numeric vector with the temporal bandwidths
#' @param kernel_name The name of the kernel to use
#' @param time_limits A vector with the upper and lower limit of the time period studied
#' @importFrom cubature cubintegrate
#' @return A numeric vector with the correction factor values for each event
#' @keywords internal
#' @examples
#' #no example provided, this is an internal function
correction_factor_time <- function(events_time, samples_time, bws_time, kernel_name, time_limits = NULL){
  kernel_func <- select_kernel(kernel_name)

  # step1 : determine the study period
  if(is.null(time_limits)){
    start_time <- min(samples_time)
    end_time <- max(samples_time)
  }else{
    start_time <- min(time_limits)
    end_time <- max(time_limits)
  }



  # step2 : calculating for each event if it is above or under the limits
  low_diff <- events_time-bws_time - start_time
  low_diff <- ifelse(low_diff < 0, abs(low_diff), 0)

  up_diff <- end_time - events_time
  up_diff <- ifelse(up_diff < 0, 0, up_diff)

  # calculating the part of the density outside the area (lower bound)
  get_integral1 <- function(bw, low_diff){
    cubintegrate(kernel_func,lower=(bw-low_diff),upper=bw,
                 bw=bw, relTol = 1e-15)$integral
  }
  get_inegral_vec <- Vectorize(get_integral1, vectorize.args = c("low_diff","bw"))
  out_lower <- get_inegral_vec(low_diff = low_diff, bw = bws_time)

  # calculating the part of the density outside the area (upper bound)
  get_integral2 <- function(bw, up_diff){
    cubintegrate(kernel_func,lower=up_diff,upper=bw,
                 bw=bw, relTol = 1e-15)$integral
  }
  get_inegral_vec <- Vectorize(get_integral2, vectorize.args = c("up_diff","bw"))
  out_upper <- get_inegral_vec(up_diff = up_diff, bw = bws_time)

  in_density <- 1 - (out_upper + out_lower)
  return(1/in_density)

  #outer_part <- out_lower + out_upper
  #return(outer_part)

}

# PREVIOUS VERSION KEPT FOR DEBUGING
# bw_cvl_calc <- function(bw_range,bw_step,lines, events, w, kernel_name, method,  diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){
#
#   ## step0 basic checks
#   samples <- events
#   div <- "bw"
#
#   if(verbose){
#     print("checking inputs ...")
#   }
#   if((kernel_name %in% c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", 'epanechnikov','uniform'))==FALSE){
#     stop('the name of the kernel function must be one of c("triangle", "gaussian", "scaled gaussian", "tricube", "cosine" ,"triweight", "quartic", "epanechnikov" ,"uniform")')
#   }
#
#   if(method %in% c("simple","continuous","discontinuous") == FALSE){
#     stop('the method must be one of c("simple","continuous","discontinuous"')
#   }
#   if(method == "continuous" & kernel_name == "gaussian"){
#     stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandwidth, please consider using the quartic kernel instead")
#   }
#
#   if(min(bw_range)<=0){
#     stop("the bandwidth for the kernel must be superior to 0")
#   }
#
#   if(bw_step<=0){
#     stop("the step between two bandwidths must be greater than 0")
#   }
#
#   if(diggle_correction & is.null(study_area)){
#     stop("the study_area must be defined if the Diggle correction factor is used")
#   }
#   if(check){
#     check_geometries(lines,samples,events, study_area)
#   }
#
#
#   ## step1 : preparing the data
#   if(verbose){
#     print("prior data preparation ...")
#   }
#   data <- prepare_data(samples, lines, events,w,digits,tol,agg)
#   lines <- data$lines
#   Wl <- rgeos::gLength(lines) / 1000
#   samples <- data$events
#   events <- data$events
#
#   ## step2  creating the grid
#   grid <- build_grid(grid_shape,list(lines,samples,events))
#
#   ## calculating the correction factor for each bw
#   all_bws <- seq(min(bw_range),max(bw_range),bw_step)
#   if(verbose){
#     print("Calculating the correction factor if required")
#   }
#   for (bw in all_bws){
#     if(diggle_correction){
#       bws <- rep(bw,nrow(events))
#       corr_factor <- correction_factor(study_area,events,lines,method,bws, kernel_name, tol, digits, max_depth, sparse)
#     }else{
#       corr_factor<- rep(1,nrow(events))
#     }
#
#     events[[paste("weight_",bw,sep="")]] <- events$weight * corr_factor
#     samples[[paste("weight_",bw,sep="")]] <- samples$weight * corr_factor
#
#   }
#   max_bw <- max(bw_range)
#
#   ## step3 splitting the dataset with each rectangle
#   selections <- split_by_grid(grid,samples,events,lines,max_bw, tol, digits, split_all = FALSE)
#
#   ## sub sampling the quadra if required
#   if (sub_sample < 1){
#     nb <- ceiling(length(selections) * sub_sample)
#     selections <- selections[sample(1:length(selections),size = nb,replace = FALSE)]
#   }
#
#   ## step 4 calculating the CVl values
#   if(verbose){
#     print("start calculating the CVl values ...")
#   }
#
#   n_quadra <- length(selections)
#
#   if (verbose){
#     pb <- txtProgressBar(min = 0, max = n_quadra, style = 3)
#   }
#   dfs <- lapply(1:n_quadra,function(i){
#     sel <- selections[[i]]
#
#     values <- nkde_worker_bw_sel_cvl(sel$lines, sel$events,
#                                  sel$samples, kernel_name, all_bws,
#                                  method, div, digits,
#                                  tol,sparse, max_depth, verbose)
#
#     if(verbose){
#       setTxtProgressBar(pb, i)
#     }
#
#     ## on combine les resultats de chaque bw dans une matrice
#     df <- data.frame(do.call(cbind,values))
#     names(df) <- paste("k",1:ncol(df),sep="")
#     df$goid <- sel$samples$goid
#     return(df)
#   })
#
#   ## step5  combining the results for each quadra
#   tot_df <- do.call(rbind,dfs)
#   tot_df <- tot_df[order(tot_df$goid),]
#
#
#   ## calculer les valeurs de CVl
#   cv_scores <- sapply(1:length(all_bws), function(i){
#     bw <- all_bws[[i]]
#     kvalues <- tot_df[,i]
#     score <- (Wl - sum(1/kvalues))**2
#     return(score)
#   })
#
#   finaldf <- data.frame(
#     "bw" = all_bws,
#     "cvl_scores" = cv_scores
#   )
#
#   return(finaldf)
# }


# nkde_worker_bw_sel_cvl <- function(lines, events, samples, kernel_name, bws, method, div, digits, tol, sparse, max_depth, verbose = FALSE){
#
#   # if we do not have event in that space, just return 0 values
#   if(nrow(events)==0){
#     values <- lapply(bws,function(i){rep(0,nrow(samples))})
#     return(values)
#   }
#
#   ## step1 creating the graph
#   graph_result <- build_graph(lines,digits = digits,line_weight = "length")
#   graph <- graph_result$graph
#   nodes <- graph_result$spvertices
#   edges <- graph_result$spedges
#
#   ## step2 finding for each event, its corresponding node
#   ## NOTE : there will be less samples than events most of the time
#   ## because of the avoidance of island effects.
#   events$vertex_id <- closest_points(events, nodes)
#   samples$vertex_id <- closest_points(samples, nodes)
#   events$oid <- 1:nrow(events)
#
#   ## step3 starting the calculations !
#   neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
#   neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})
#
#   ## we calculate the nkde values for each bw provided
#   kernel_values <- lapply(bws, function(bw){
#
#     repbws <- rep(bw,nrow(events))
#
#     if(method == "simple"){
#       values <- spNetwork::get_loo_values_simple(neighbour_list, samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
#                                                  events$vertex_id, events[[paste("weight_",bw,sep="")]],
#                                                  repbws, kernel_name, graph_result$linelist, max_depth)
#     }else if(method=="continuous"){
#       ##and finally calculating the values
#       # NOTE : Values is a dataframe with two columns :
#       # the kvalue calculated for each event when all the events are considered
#       # the kvalue specific at each event (when all other events are not considered)
#       values <- spNetwork::get_loo_values_continuous(neighbour_list,  samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
#                                                      events$vertex_id, events[[paste("weight_",bw,sep="")]],
#                                                      repbws, kernel_name, graph_result$linelist, max_depth)
#       values <- values$sum_k
#
#     }else if(method == "discontinuous"){
#       values <- spNetwork::get_loo_values_discontinuous(neighbour_list, samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
#                                                         events$vertex_id, events[[paste("weight_",bw,sep="")]], repbws,
#                                                         kernel_name, graph_result$linelist, max_depth)
#     }
#
#     ## step7 adjusting the kernel values !
#     # dividing by bw is crucial, other wise, larger BW are always better !
#     return(values * (1/bw))
#   })
#
#   ## at that point, we have a list of numeric vectors for each bw
#   return(kernel_values)
# }



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### EN DEVELOPPEMENT ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Bandwidth selection by Cronie and Van Lieshout's Criterion
#'
#' @description Calculate for multiple bandwidth the Cronie and Van Lieshout's Criterion to
#' select an appropriate bandwidth in a data-driven approach.
#'
#' @details The Cronie and Van Lieshout's Criterion \insertCite{cronie2018non}{spNetwork}
#' find the optimal bandwidth by minimizing the difference between the size of the observation
#' window and the sum of the reciprocal of the estimated kernel density at the
#' events locations. In the network case, the size of the study area is the sum
#' of the length of each line in the network. Thus, it is important to only
#' use the necessary parts of the network.
#'
#' @references{
#'     \insertAllCited{}
#' }
#'
#' @template bw_selection-args
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @template nkde_geoms-args
#' @template sparse-arg
#' @template grid_shape-arg
#' @template bw_selection_adapt-args
#' @param sub_sample A float between 0 and 1 indicating the percentage of quadra
#' to keep in the calculus. For large datasets, it may be useful to limit the
#' bandwidth evaluation and thus reduce calculation time.
#' @param verbose A Boolean, indicating if the function should print messages
#' about the process.
#' @param zero_strat A string indicating what to do when density is 0 when calculating LOO density estimate for an isolated event.
#' "min_double" (default) replace the 0 value by the minimum double possible on the machine. "remove" will remove them from the final
#' score. The first approach penalizes more strongly the small bandwidths.
#' @template check-arg
#' @return A dataframe with two columns, one for the bandwidths and the second for
#' the Cronie and Van Lieshout's Criterion.
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#' cv_scores <- bw_cvl_calc(seq(200,400,50),
#'                                mtl_network, bike_accidents,
#'                                rep(1,nrow(bike_accidents)),
#'                                "quartic", "discontinuous",
#'                                diggle_correction = FALSE, study_area = NULL,
#'                                max_depth = 8,
#'                                digits=2, tol=0.1, agg=5,
#'                                sparse=TRUE, grid_shape=c(1,1),
#'                                sub_sample = 1, verbose=TRUE, check=TRUE)
#'}
bw_cvl_calc <- function(bws = NULL,
                        lines, events, w, kernel_name, method,
                        diggle_correction = FALSE, study_area = NULL,
                        adaptive = FALSE, trim_bws = NULL, mat_bws = NULL,
                        max_depth = 15,
                        digits=5, tol=0.1, agg=NULL, sparse=TRUE,
                        zero_strat = "min_double",
                        grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  div <- "bw"
  events$weight <- w
  events$wid <- 1:nrow(events)

  wl <- sum(as.numeric(st_length(lines)))

  if(verbose){
    print("checking inputs ...")
  }

  passed <- bw_checks(check,lines,samples,events,
                      kernel_name, method, bws_net = bws, bws_time = NULL,
                      adaptive = adaptive, trim_net_bws = trim_bws,
                      diggle_correction = diggle_correction, study_area = study_area)

  if(zero_strat %in% c("min_double", "remove") == FALSE){
    stop("zero_strat argument must be one of c('min_double', 'remove')")
  }

  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }

  data <- prepare_data(samples, lines, events, w , digits,tol,agg)
  lines <- data$lines
  events_loc <- data$events

  #idx <- FNN::knnx.index(st_coordinates(events_loc),st_coordinates(events), k = 1)
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## calculating the correction factor for each bw
  ## they must be calculate for the location of the events and then stored in a matrix
  # all_bws <- seq(min(bw_range),max(bw_range),bw_step)
  all_bws <- bws

  if(verbose){
    print("Calculating the correction factor if required")
  }

  ## we will construct here a matrix with all the bandwidths
  ## if we are not in the adaptive mode, then all the columns will have
  ## unique values
  if(adaptive == FALSE){
    mat_bws <- sapply(all_bws, function(x){
      rep(x,nrow(events))
    })
  }else{
    if(is.null(mat_bws)){
      mat_bws <- adaptive_bw(grid, events, lines, all_bws, trim_bws, method,
                             kernel_name, max_depth, tol, digits, sparse, verbose)
    }else{

      # in the case where all the bandwidths were provided by the user
      # we only have to extract the relevant informations from this matrix
      all_bws <- colnames(mat_bws)
      if(is.null(all_bws)){
        all_bws <- 1:ncol(mat_bws)
      }

      if(nrow(mat_bws) != nrow(events)){
        stop("The number of rows in mat_bws must be the same as the number of rows in events")
      }
    }
  }

  events_weight <- apply(mat_bws, MARGIN = 2, FUN = function(bws){

    if(diggle_correction){
      corr_factor <- correction_factor(study_area,events_loc,lines,method, bws, kernel_name, tol, digits, max_depth, sparse)
      corr_factor <- corr_factor[events$goid] * events$weight

    }else{
      corr_factor <- rep(1,nrow(events)) * events$weight
    }
    return(corr_factor)
  })


  max_bw <- max(mat_bws)

  ## step3 splitting the dataset with each rectangle
  # NB : here we select the events in the quadra (samples) and the events locations in the buffer (events_loc)
  selections <- split_by_grid(grid, events, events_loc, lines,max_bw, tol, digits, split_all = FALSE)

  ## sub sampling the quadra if required
  if (sub_sample < 1){
    nb <- ceiling(length(selections) * sub_sample)
    selections <- selections[sample(1:length(selections),size = nb,replace = FALSE)]
  }

  ## step 4 calculating the CV values
  if(verbose){
    print("start calculating the CV values ...")
  }

  n_quadra <- length(selections)

  if (verbose){
    pb <- txtProgressBar(min = 0, max = n_quadra, style = 3)
  }
  dfs <- lapply(1:n_quadra,function(i){

    sel <- selections[[i]]

    # the events_loc must cover the quadra and the bw
    sel_events_loc <- sel$events

    # idem for all the events
    sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

    # but I also need to know on which events I must calculate the densities (in the quadra)
    quad_events <- sel$samples
    sel_weights <- events_weight[sel_events$wid,]

    # I extract here the bws required
    sel_bws <- mat_bws[sel_events$wid,]

    values <- nkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                  kernel_name, sel_bws,
                                  method, div, digits,
                                  tol,sparse, max_depth, verbose, cvl = TRUE)
    # NOTE : values is a matrix

    return(values)

  })


  # removing NULL elements in list
  dfs[sapply(dfs, is.null)] <- NULL

  # all the elements are matrices, we must combine them by row
  all_loo_scores <- do.call(rbind, dfs)

  if(zero_strat == "min_double"){
    all_loo_scores <- ifelse(all_loo_scores == 0, .Machine$double.xmin, all_loo_scores)
  }else{
    all_loo_scores <- ifelse(all_loo_scores == 0, NA, all_loo_scores)
  }

  part2 <- colSums(all_loo_scores ** (-1), na.rm = TRUE)

  cv_scores <- (part2 - wl)**2
  bw_scores <- data.frame(
    "bw" = all_bws,
    "cvl_scores" = cv_scores
  )

  return(bw_scores)
}


#' @title Bandwidth selection by Cronie and Van Lieshout's Criterion (multicore version)
#'
#' @description Calculate for multiple bandwidths the Cronie and Van Lieshout's Criterion to
#' select an appropriate bandwidth in a data-driven approach. A plan from the package future can be used
#' to split the work across several cores. The different cells generated in accordance with the
#' argument grid_shape are used for the parallelization. So if only one cell is
#' generated (`grid_shape = c(1,1)`), the function will use only one core. The progress bar
#' displays the progression for the cells.
#'
#'
#' @details For more details, see help(bw_cvl_calc)
#'
#' @template bw_selection-args
#' @template nkde_params-arg
#' @template diggle_corr-arg
#' @template nkde_geoms-args
#' @template sparse-arg
#' @template grid_shape-arg
#' @template bw_selection_adapt-args
#' @param adaptive A boolean indicating if an adaptive bandwidth must be used.
#' If adaptive = TRUE, the local bandwidth are derived from the global bandwidths
#' calculated from bw_range and bw_step.
#' @param trim_bws A vector indicating the maximum value an adaptive bandwidth can
#' reach. Higher values will be trimmed. It must have the same length as `bws`.
#' @param sub_sample A float between 0 and 1 indicating the percentage of quadra
#' to keep in the calculus. For large datasets, it may be useful to limit the
#' bandwidth evaluation and thus reduce calculation time.
#' @param verbose A Boolean, indicating if the function should print messages
#' about the process.
#' @template check-arg
#' @param zero_strat A string indicating what to do when density is 0 when calculating LOO density estimate for an isolated event.
#' "min_double" (default) replace the 0 value by the minimum double possible on the machine. "remove" will remove them from the final
#' score. The first approach penalizes more strongly the small bandwidths.
#' @return A dataframe with two columns, one for the bandwidths and the second for
#' the Cronie and Van Lieshout's Criterion.
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multisession(workers=1))
#' cv_scores <- bw_cvl_calc.mc(seq(200,400,50),
#'                                mtl_network, bike_accidents,
#'                                rep(1,nrow(bike_accidents)),
#'                                "quartic", "discontinuous",
#'                                diggle_correction = FALSE, study_area = NULL,
#'                                max_depth = 8,
#'                                digits=2, tol=0.1, agg=5,
#'                                sparse=TRUE, grid_shape=c(1,1),
#'                                sub_sample = 1, verbose=TRUE, check=TRUE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
bw_cvl_calc.mc <- function(bws = NULL, lines, events, w, kernel_name, method,
                           diggle_correction = FALSE, study_area = NULL,
                           adaptive = FALSE, trim_bws = NULL, mat_bws = NULL,
                           max_depth = 15,
                           digits=5, tol=0.1, agg=NULL, sparse=TRUE,
                           zero_strat = "min_double",
                           grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  div <- "bw"
  events$weight <- w
  events$wid <- 1:nrow(events)

  wl <- sum(as.numeric(st_length(lines)))

  if(verbose){
    print("checking inputs ...")
  }

  passed <- bw_checks(check,lines,samples,events,
                      kernel_name, method, bws_net = bws, bws_time = NULL,
                      adaptive = adaptive, trim_net_bws = trim_bws,
                      diggle_correction = diggle_correction, study_area = study_area)


  ## step1 : preparing the data
  if(verbose){
    print("prior data preparation ...")
  }

  data <- prepare_data(samples, lines, events, w , digits,tol,agg)
  lines <- data$lines
  events_loc <- data$events

  #idx <- FNN::knnx.index(st_coordinates(events_loc),st_coordinates(events), k = 1)
  idx <- closest_points(events, events_loc)
  events$goid <- events_loc$goid[idx]

  ## step2  creating the grid
  grid <- build_grid(grid_shape,list(lines,samples,events))

  ## calculating the correction factor for each bw
  ## they must be calculate for the location of the events and then stored in a matrix
  # all_bws <- seq(min(bw_range),max(bw_range),bw_step)
  all_bws <- bws

  if(verbose){
    print("Calculating the correction factor if required")
  }

  ## we will construct here a matrix with all the bandwidths
  ## if we are not in the adaptive mode, then all the columns will have
  ## unique values
  if(adaptive == FALSE){
    mat_bws <- sapply(all_bws, function(x){
      rep(x,nrow(events))
    })
  }else{
    if(is.null(mat_bws)){
      mat_bws <- adaptive_bw(grid, events, lines, all_bws, trim_bws, method,
                             kernel_name, max_depth, tol, digits, sparse, verbose)
    }else{

      # in the case where all the bandwidths were provided by the user
      # we only have to extract the relevant informations from this matrix
      all_bws <- colnames(mat_bws)
      if(is.null(all_bws)){
        all_bws <- 1:ncol(mat_bws)
      }

      if(nrow(mat_bws) != nrow(events)){
        stop("The number of rows in mat_bws must be the same as the number of rows in events")
      }
    }
  }

  events_weight <- apply(mat_bws, MARGIN = 2, FUN = function(bws){

    if(diggle_correction){
      corr_factor <- correction_factor(study_area,events_loc,lines,method, bws, kernel_name, tol, digits, max_depth, sparse)
      corr_factor <- corr_factor[events$goid] * events$weight

    }else{
      corr_factor <- rep(1,nrow(events)) * events$weight
    }
    return(corr_factor)
  })

  max_bw <- max(mat_bws)
  ## step3 splitting the dataset with each rectangle
  # NB : here we select the events in the quadra (samples) and the events locations in the buffer (events_loc)
  selections <- split_by_grid(grid, events, events_loc, lines,max_bw, tol, digits, split_all = FALSE)

  ## sub sampling the quadra if required
  if (sub_sample < 1){
    nb <- ceiling(length(selections) * sub_sample)
    selections <- selections[sample(1:length(selections),size = nb,replace = FALSE)]
  }

  ## step 4 calculating the CV values
  if(verbose){
    print("start calculating the CV values ...")
  }

  n_quadra <- length(selections)


  if(verbose){
    progressr::with_progress({
      p <- progressr::progressor(along = selections)
      dfs <- future.apply::future_lapply(selections, function(sel) {

        sel_events_loc <- sel$events

        # idem for all the events
        sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

        # but I also need to know on which events I must calculate the densities (in the quadra)
        quad_events <- sel$samples
        sel_weights <- events_weight[sel_events$wid,]

        # I extract here the bws required
        sel_bws <- mat_bws[sel_events$wid,]

        values <- nkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                     kernel_name, sel_bws,
                                     method, div, digits,
                                     tol,sparse, max_depth, verbose, cvl = TRUE)

        p(sprintf("i=%g", sel$index))

        return(values)

      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      sel_events_loc <- sel$events

      # idem for all the events
      sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

      # but I also need to know on which events I must calculate the densities (in the quadra)
      quad_events <- sel$samples
      sel_weights <- events_weight[sel_events$wid,]

      # I extract here the bws required
      sel_bws <- mat_bws[sel_events$wid,]

      values <- nkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                   kernel_name, sel_bws,
                                   method, div, digits,
                                   tol,sparse, max_depth, verbose, cvl = TRUE)
      return(values)
    })
  }

  # removing NULL elements in list
  dfs[sapply(dfs, is.null)] <- NULL

  # all the elements are matrices, we must combine them by row
  all_loo_scores <- do.call(rbind, dfs)

  if(zero_strat == "min_double"){
    all_loo_scores <- ifelse(all_loo_scores == 0, .Machine$double.xmin, all_loo_scores)
  }else{
    all_loo_scores <- ifelse(all_loo_scores == 0, NA, all_loo_scores)
  }

  part2 <- colSums(all_loo_scores ** (-1), na.rm = TRUE)

  cv_scores <- (part2 - wl)**2

  finaldf <- data.frame(
    "bw" = all_bws,
    "cvl_scores" = cv_scores
  )

  # # removing NULL elements in list
  # dfs[sapply(dfs, is.null)] <- NULL
  #
  # add <- function(x) Reduce("+", x)
  #
  # cv_scores <- (add(dfs) - wl)**2
  #
  #
  # finaldf <- data.frame(
  #   "bw" = all_bws,
  #   "cvl_scores" = cv_scores
  # )

  return(finaldf)
}


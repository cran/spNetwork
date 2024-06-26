# THis is the older verison of the function (kept for debuging)
# bw_cv_likelihood_calc <- function(bw_range,bw_step,lines, events, w, kernel_name, method,  diggle_correction = FALSE, study_area = NULL, max_depth = 15, digits=5, tol=0.1, agg=NULL, sparse=TRUE, grid_shape=c(1,1), sub_sample=1, verbose=TRUE, check=TRUE){
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
#     stop("using the continuous NKDE and the gaussian kernel function can yield negative values for densities because the gaussian kernel does not integrate to 1 within the bandiwdth, please consider using the quartic kernel instead")
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
#   ## step 4 calculating the CV values
#   if(verbose){
#     print("start calculating the CV values ...")
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
#     values <- nkde_worker_bw_sel(sel$lines, sel$events,
#                           sel$samples, kernel_name, all_bws,
#                           method, div, digits,
#                           tol,sparse, max_depth, verbose)
#
#     if(verbose){
#       setTxtProgressBar(pb, i)
#     }
#
#     # si on est en mode continu, on obtient une liste de dataframes
#     # sinon une liste de vecteurs numeriques
#     if(method == "continuous"){
#       dfk <- data.frame(do.call(cbind,lapply(values,function(v){v$kvalues})))
#       dfloo <- data.frame(do.call(cbind,lapply(values,function(v){v$loovalues})))
#       dfk$goid <- sel$samples$goid
#       dfloo$goid <- sel$samples$goid
#       return(list(dfk,dfloo))
#
#     }else{
#       ## on combine les resultats de chaque bw dans une matrice
#       df <- data.frame(do.call(cbind,values))
#       names(df) <- paste("k",1:ncol(df),sep="")
#       df$goid <- sel$samples$goid
#       return(df)
#     }
#
#
#   })
#
#   ## step5  combining the results for each quadra
#   if(method == "continuous"){
#     tot_df <- do.call(rbind,lapply(dfs,function(df){df[[1]]}))
#     tot_df <- tot_df[order(tot_df$goid),]
#     tot_loos <- do.call(rbind,lapply(dfs,function(df){df[[2]]}))
#     tot_loos <- tot_loos[order(tot_loos$goid),]
#   }else{
#     tot_df <- do.call(rbind,dfs)
#     tot_df <- tot_df[order(tot_df$goid),]
#   }
#
#
#   ## calculer les valeurs de CV
#   kern_fun <- select_kernel(kernel_name)
#   cv_scores <- sapply(1:length(all_bws), function(i){
#     bw <- all_bws[[i]]
#     kvalues <- tot_df[,i]
#     if (method %in% c("simple","discontinuous")){
#       correction <- kern_fun(0,bw) * (1/bw)
#     } else{
#       correction <- tot_loos[,i]
#     }
#     cv_values <- kvalues - correction
#     cv_values[cv_values==0] <- .Machine$double.xmin
#     score <- sum(log(cv_values))
#     return(score)
#   })
#
#   finaldf <- data.frame(
#     "bw" = all_bws,
#     "cv_scores" = cv_scores
#   )
#
#   return(finaldf)
# }


# nkde_worker_bw_sel <- function(lines, events, samples, kernel_name, bws, method, div, digits, tol, sparse, max_depth, verbose = FALSE){
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
#
#     }else if(method == "discontinuous"){
#       values <- spNetwork::get_loo_values_discontinuous(neighbour_list, samples$vertex_id, samples[[paste("weight_",bw,sep="")]],
#                                                         events$vertex_id, events[[paste("weight_",bw,sep="")]], repbws,
#                                                         kernel_name, graph_result$linelist, max_depth)
#     }
#
#     ## step7 adjusting the kernel values !
#     # dividing by bw is crucial, other wise, larger BW are always better !
#     if(method == "continuous"){
#       df <- data.frame(
#         kvalues = values$sum_k * (1/bw),
#         loovalues = values$loo * (1/bw)
#       )
#       return(df)
#     }else {
#       return(values * (1/bw))
#     }
#   })
#
#   ## at that point, we have a list of numeric vectors or a list of dataframes, one for each bw
#   return(kernel_values)
# }


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### CV LIKELIHOOD ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Bandwidth selection by likelihood cross validation
#'
#' @description Calculate for multiple bandwidth the cross validation likelihood to
#' select an appropriate bandwidth in a data-driven approach
#'
#' @details  The function calculates the likelihood cross validation score for several
#' bandwidths in order to find the most appropriate one. The general idea is to find the
#' bandwidth that would produce the most similar results if one event was removed from
#' the dataset (leave one out cross validation). We use here the shortcut formula as
#' described by the package spatstat \insertCite{spatstatpkg}{spNetwork}.
#'
#'  \eqn{LCV(h) = \sum_i \log\hat\lambda_{-i}(x_i)}
#'
#' Where the sum is taken for all events \eqn{x_i} and where \eqn{\hat\lambda_{-i}(x_i)} is the leave-one-out kernel
#' estimate at \eqn{x_i} for a bandwidth h. A higher value indicates a better bandwidth.
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
#' the cross validation score (the lower the better).
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#' cv_scores <- bw_cv_likelihood_calc(seq(200,800,50),
#'                                mtl_network, bike_accidents,
#'                                rep(1,nrow(bike_accidents)),
#'                                "quartic", "simple",
#'                                diggle_correction = FALSE, study_area = NULL,
#'                                max_depth = 8,
#'                                digits=2, tol=0.1, agg=5,
#'                                sparse=TRUE, grid_shape=c(1,1),
#'                                sub_sample = 1, verbose=TRUE, check=TRUE)
#' }
bw_cv_likelihood_calc <- function(bws = NULL,
                                  lines, events, w, kernel_name, method,
                                  diggle_correction = FALSE, study_area = NULL,
                                  adaptive = FALSE, trim_bws = NULL, mat_bws = NULL,
                                  max_depth = 15, digits=5, tol=0.1, agg=NULL,
                                  sparse=TRUE, grid_shape=c(1,1), sub_sample=1,
                                  zero_strat = "min_double",
                                  verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  div <- "bw"
  events$weight <- w
  events$wid <- 1:nrow(events)

  if(verbose){
    print("checking inputs ...")
  }

  passed <- bw_checks(check,lines,samples,events,
           kernel_name, method, bws_net = bws,
           adaptive = adaptive, trim_net_bws = trim_bws, arr_bws_net = mat_bws,
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
  ## they must be calculated for the location of the events and then stored in a matrix
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
      dim(mat_bws) <- c(nrow(events),length(all_bws))
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
  if(verbose){
    print("Splittig the dataset within the grid ...")
  }
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
    print(i)
    sel <- selections[[i]]

    # the events_loc must cover the quadra and the bw
    sel_events_loc <- sel$events

    # idem for all the events
    sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

    # but I also need to know on which events I must calculate the densities (in the quadra)
    quad_events <- sel$samples
    sel_weights <- events_weight[sel_events$wid,,drop=FALSE]

    # I extract here the bws required
    sel_bws <- mat_bws[sel_events$wid,,drop=FALSE]

    values <- nkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                  kernel_name, sel_bws,
                                  method, div, digits,
                                  tol,sparse, max_depth, verbose)

    if(verbose){
      setTxtProgressBar(pb, i)
    }

    return(values)

  })

  # removing NULL elements in list
  dfs[sapply(dfs, is.null)] <- NULL

  # all the elements are matrices, we must combine them by row
  all_loo_scores <- do.call(rbind, dfs)


  # and we can calculate now the scores
  if(zero_strat == "min_double"){
    all_loo_scores <- ifelse(all_loo_scores <= 0, .Machine$double.xmin, all_loo_scores)
    cv_scores <- colSums(log(all_loo_scores)) / nrow(all_loo_scores)
  }else{
    binary_mat <- all_loo_scores <= 0
    all_loo_scores <- ifelse(binary_mat, 1, all_loo_scores)
    cv_scores <- colSums(log(all_loo_scores)) / colSums(binary_mat == FALSE)
  }

  # add <- function(x) Reduce("+", x)
  #
  # cv_scores <- add(dfs)


  finaldf <- data.frame(
    "bw" = all_bws,
    "cv_scores" = cv_scores
  )

  return(finaldf)
}



#' @title Bandwidth selection by likelihood cross validation (multicore)
#'
#' @description Calculate for multiple bandwidth the cross validation likelihood to
#' select an appropriate bandwidth in a data-driven approach
#'
#' @details  See the function bw_cv_likelihood_calc for more details. The calculation is split
#' according to the parameter grid_shape. If `grid_shape = c(1,1)`, then parallel processing cannot be used.
#'
#' @template bw_selection-args
#' @template diggle_corr-arg
#' @template nkde_params-arg
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
#' the cross validation score (the lower the better).
#' @export
#' @examples
#' \donttest{
#' data(mtl_network)
#' data(bike_accidents)
#' future::plan(future::multisession(workers=1))
#' cv_scores <- bw_cv_likelihood_calc.mc(seq(200,800,50),
#'                                mtl_network, bike_accidents,
#'                                rep(1,nrow(bike_accidents)),
#'                                "quartic", "simple",
#'                                diggle_correction = FALSE, study_area = NULL,
#'                                max_depth = 8,
#'                                digits=2, tol=0.1, agg=5,
#'                                sparse=TRUE, grid_shape=c(1,1),
#'                                sub_sample = 1, verbose=TRUE, check=TRUE)
#' ## make sure any open connections are closed afterward
#' if (!inherits(future::plan(), "sequential")) future::plan(future::sequential)
#' }
bw_cv_likelihood_calc.mc <- function(bws, lines, events, w, kernel_name, method,
                                     diggle_correction = FALSE, study_area = NULL,
                                     adaptive = FALSE, trim_bws = NULL, mat_bws = NULL,
                                     max_depth = 15, digits=5, tol=0.1, agg=NULL,
                                     sparse=TRUE, grid_shape=c(1,1), sub_sample=1,
                                     zero_strat = "min_double",
                                     verbose=TRUE, check=TRUE){

  ## step0 basic checks
  samples <- events
  div <- "bw"
  events$weight <- w
  events$wid <- 1:nrow(events)

  if(verbose){
    print("checking inputs ...")
  }

  passed <- bw_checks(check,lines,samples,events,
                      kernel_name, method, bws_net = bws, bws_time = NULL,
                      adaptive = adaptive, trim_net_bws = trim_bws, arr_bws_net = mat_bws,
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
      mat_bws <- adaptive_bw.mc(grid = grid,
                                events = events,
                                lines = lines,
                                bw = all_bws,
                                trim_bw = trim_bws,
                                method = method,
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

        # the events_loc must cover the quadra and the bw
        sel_events_loc <- sel$events

        # idem for all the events
        sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

        # but I also need to know on which events I must calculate the densities (in the quadra)
        quad_events <- sel$samples
        sel_weights <- events_weight[sel_events$wid,,drop=FALSE]

        # I extract here the bws required
        sel_bws <- mat_bws[sel_events$wid,,drop=FALSE]

        values <- nkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                     kernel_name, sel_bws,
                                     method, div, digits,
                                     tol,sparse, max_depth, verbose)

        p(sprintf("i=%g", sel$index))

        return(values)

      })
    })
  }else{
    dfs <- future.apply::future_lapply(selections, function(sel) {

      # the events_loc must cover the quadra and the bw
      sel_events_loc <- sel$events

      # idem for all the events
      sel_events <- subset(events, events$goid %in% sel_events_loc$goid)

      # but I also need to know on which events I must calculate the densities (in the quadra)
      quad_events <- sel$samples
      sel_weights <- events_weight[sel_events$wid,,drop=FALSE]

      # I extract here the bws required
      sel_bws <- mat_bws[sel_events$wid,,drop=FALSE]

      values <- nkde_worker_bw_sel(sel$lines, quad_events, sel_events_loc, sel_events, sel_weights,
                                   kernel_name, sel_bws,
                                   method, div, digits,
                                   tol,sparse, max_depth, verbose)

      return(values)
    })
  }


  # # removing NULL elements in list
  # dfs[sapply(dfs, is.null)] <- NULL
  #
  # add <- function(x) Reduce("+", x)
  #
  # cv_scores <- add(dfs)
  #
  #
  # finaldf <- data.frame(
  #   "bw" = all_bws,
  #   "cv_scores" = cv_scores
  # )

  # removing NULL elements in list
  dfs[sapply(dfs, is.null)] <- NULL

  # all the elements are matrices, we must combine them by row
  all_loo_scores <- do.call(rbind, dfs)

  # and we can calculate now the scores
  if(zero_strat == "min_double"){
    all_loo_scores <- ifelse(all_loo_scores <= 0, .Machine$double.xmin, all_loo_scores)
    cv_scores <- colSums(log(all_loo_scores)) / nrow(all_loo_scores)
  }else{
    binary_mat <- all_loo_scores <= 0
    all_loo_scores <- ifelse(binary_mat, 1, all_loo_scores)
    cv_scores <- colSums(log(all_loo_scores)) / colSums(binary_mat == FALSE)
  }

  # add <- function(x) Reduce("+", x)
  #
  # cv_scores <- add(dfs)


  finaldf <- data.frame(
    "bw" = all_bws,
    "cv_scores" = cv_scores
  )

  return(finaldf)
}




#' @title Bandwidth selection by likelihood cross validation worker function
#'
#' @description worker function for calculating for multiple bandwidth the cross validation likelihood to
#' select an appropriate bandwidth in a data-driven approach
#'
#' @param lines A feature collection of linestrings representing the underlying network
#' @param quad_events a feature collection of points indicating for which events the densities must be calculated
#' @param events_loc A feature collection of points representing the location of the events
#' @param events A feature collection of points representing the events. Multiple events can share
#' the same location. They are linked by the goid column
#' @param w A numeric matrix with the weight of the events for each bandwdith
#' @param kernel_name The name of the kernel to use (string)
#' @param bws_net A numeric matrix with the network bandwidths for each event
#' @param method The type of NKDE to use (string)
#' @param zero_strat A string indicating what to do when density is 0 when calculating LOO density estimate for an isolated event.
#' "min_double" (default) replace the 0 value by the minimum double possible on the machine. "remove" will remove them from the final
#' score. The first approach penalizes more strongly the small bandwidths.
#' @param digits The number of digits to retain from the spatial coordinates. It
#'   ensures that topology is good when building the network. Default is 3. Too high a
#'   precision (high number of digits) might break some connections
#' @param tol A float indicating the minimum distance between the events and the
#'   lines' extremities when adding the point to the network. When points are
#'   closer, they are added at the extremity of the lines.
#' @template sparse-arg
#' @param verbose A boolean
#' @param cvl A boolean indicating if the cvl method (TRUE) or the loo (FALSE) method must be used
#' @keywords internal
#' @examples
#' # no example provided, this is an internal function
nkde_worker_bw_sel <- function(lines, quad_events, events_loc, events, w,
                               kernel_name, bws_net, method, div,
                               digits, tol, sparse, max_depth,
                               zero_strat = "min_double",
                               verbose = FALSE, cvl = FALSE){


  # if we do not have event in that space, just return NULL
  if(nrow(events)==0){
    return(NULL)
  }

  ## step1 creating the graph
  graph_result <- build_graph(lines,digits = digits, line_weight = "length")
  graph <- graph_result$graph
  nodes <- graph_result$spvertices
  edges <- graph_result$spedges

  ## step2 finding for each event, its corresponding node
  ## NOTE : there will be less samples than events most of the time
  events_loc$vertex_id <- closest_points(events_loc, nodes)

  events_loc2 <- events_loc
  events_loc2$geometry <- NULL
  events2 <- events
  events2$geometry <- NULL
  quad_events2 <- quad_events
  quad_events2$geometry <- NULL

  #first a join for all the events in the bw
  vertex_id <- NULL #avoid a note
  i.vertex_id<- NULL #avoid a note
  setDT(events2)[events_loc2, on = "goid", vertex_id := i.vertex_id]

  #and a second join for the quad_events
  setDT(quad_events2)[events_loc2, on = "goid", vertex_id := i.vertex_id]

  ## step3 starting the calculations !
  neighbour_list <- adjacent_vertices(graph,nodes$id,mode="out")
  neighbour_list <- lapply(neighbour_list,function(x){return (as.numeric(x))})

  if(nrow(events2) == 1){
    w <- matrix(w, ncol = ncol(bws_net))
  }

  # do not forget to add spNetwork:: below
  kernel_values <- nkde_get_loo_values(method,
                                        neighbour_list,
                                        quad_events2$vertex_id, quad_events2$wid,
                                        events2$vertex_id, events2$wid, w,
                                        bws_net,
                                        kernel_name, graph_result$linelist, max_depth,
                                        cvl)
  # kernel values are a matrix

  return(kernel_values)
}

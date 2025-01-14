## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
backup_option <- options()
base_wd <- getwd()
library(ggplot2)

## ----include=FALSE------------------------------------------------------------
load(system.file("extdata", "results_vignette_kfunc.rda",
                           package = "spNetwork", mustWork = TRUE))

## ----fig.show='hold', fig.align = 'center', warning=FALSE, message=FALSE------
library(spNetwork)
library(tmap)

data(main_network_mtl)
data(mtl_libraries)
data(mtl_theatres)

tm_shape(main_network_mtl) + 
  tm_lines("black") + 
  tm_shape(mtl_libraries) + 
  tm_dots(col = "red", size = 0.2) +
  tm_shape(mtl_theatres) + 
  tm_dots(col = "blue", size = 0.2)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  kfun_theatre <- kfunctions(lines = main_network_mtl,
#                             points = mtl_theatres,
#                             start = 0,
#                             end = 5000,
#                             step = 50,
#                             width = 1000,
#                             nsim = 50,
#                             resolution = 50,
#                             verbose = FALSE,
#                             conf_int = 0.05,
#                             digits = 2,
#                             tol = 0.1,
#                             agg = NULL,
#                             return_sims = FALSE,
#                             calc_g_func = TRUE)
#  kfun_theatre$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
plot_df <- kfun_theatre$values

ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical K-function")

## ----eval = FALSE-------------------------------------------------------------
#  kfun_theatre$plotg

## ----echo=FALSE, fig.show='hold', fig.align = 'center'------------------------
ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
                fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
    labs(x = "distances",
         y = "empirical G-function")

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  kfun_biblio <- kfunctions(main_network_mtl,
#                            mtl_libraries,
#                            start = 0,
#                            end = 5000,
#                            step = 50,
#                            width = 1000,
#                            nsim = 50,
#                            resolution = 50,
#                            verbose = FALSE,
#                            conf_int = 0.05,
#                            digits = 2,
#                            tol = 0.1,
#                            agg = NULL,
#                            return_sims = FALSE,
#                            calc_g_func = TRUE)
#  
#  kfun_biblio$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
plot_df <- kfun_biblio$values

ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical K-function")

## ----eval = FALSE-------------------------------------------------------------
#  kfun_biblio$plotg

## ----echo=FALSE, fig.show='hold', fig.align = 'center'------------------------
ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_g", ymax = "upper_g"),
                fill = grDevices::rgb(0.1,0.1,0.1),alpha=0.4, )+
    geom_path(aes_string(x = "distances", y = "lower_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_g"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_g"), col="blue")+
    labs(x = "distances",
         y = "empirical G-function")

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  cross_biblio_theatre <- cross_kfunctions(
#    lines = main_network_mtl,
#    pointsA = mtl_libraries,
#    pointsB = mtl_theatres,
#    start = 0,
#    end = 5000,
#    step = 50,
#    width = 1000,
#    nsim = 50,
#    conf_int = 0.05,
#    digits = 2,
#    tol = 0.01,
#    resolution = NULL,
#    agg = NULL,
#    return_sims = FALSE,
#    calc_g_func = TRUE,
#    verbose = FALSE)
#  
#  cross_biblio_theatre$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
plot_df <- cross_biblio_theatre$values

ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical K-function")

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  cross_theatre_biblio <- cross_kfunctions(
#    lines = main_network_mtl,
#    pointsA = mtl_theatres,
#    pointsB = mtl_libraries,
#    start = 0,
#    end = 5000,
#    step = 50,
#    width = 1000,
#    nsim = 50,
#    conf_int = 0.05,
#    digits = 2,
#    tol = 0.01,
#    resolution = NULL,
#    agg = NULL,
#    return_sims = FALSE,
#    calc_g_func = TRUE,
#    verbose = FALSE)
#  
#  cross_theatre_biblio$plotk

## ----echo=FALSE, fig.align="center", fig.show='hold', message=FALSE, warning=FALSE----
plot_df <- cross_theatre_biblio$values

ggplot(plot_df)+
    geom_ribbon(aes_string(x = "distances", ymin="lower_k",ymax = "upper_k"),
                fill = grDevices::rgb(0.1,0.1,0.1), alpha=0.4)+
    geom_path(aes_string(x = "distances", y = "lower_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "upper_k"), col="black",
              linetype="dashed")+
    geom_path(aes_string(x = "distances", y = "obs_k"), col="blue")+
    labs(x = "distances",
         y = "empirical K-function")

## ----include = FALSE----------------------------------------------------------
# reset all the user parameters
options(backup_option)
setwd(base_wd)
oldpar <- par(mfrow = c(1,2))
par(oldpar)

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  
#  future::plan(future::multisession, workers = 1)
#  
#  cross_biblio_theatre_mc <- cross_kfunctions.mc(
#    lines = main_network_mtl,
#    pointsA = mtl_libraries,
#    pointsB = mtl_theatres,
#    start = 0,
#    end = 5000,
#    step = 50,
#    width = 1000,
#    nsim = 50,
#    conf_int = 0.05,
#    digits = 2,
#    tol = 0.01,
#    resolution = NULL,
#    agg = NULL,
#    return_sims = FALSE,
#    calc_g_func = TRUE,
#    verbose = TRUE,
#    grid_shape = c(2,2)
#    )

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE--------

data("bike_accidents")
data("mtl_network")

tm_shape(mtl_network) + 
  tm_lines("black") + 
  tm_shape(bike_accidents) + 
  tm_dots(col = "red", size = 0.2)


## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, eval = FALSE----
#  
#  # converting the Date field to a numeric field
#  # (counting days from first accident)
#  
#  bike_accidents$Time <- as.POSIXct(bike_accidents$Date, format = "%Y/%m/%d")
#  start <- as.POSIXct("2016/01/01", format = "%Y/%m/%d")
#  bike_accidents$Time <- difftime(bike_accidents$Time, start, units = "days")
#  bike_accidents$Time <- as.numeric(bike_accidents$Time)
#  
#  
#  net_time_values <- k_nt_functions(
#    lines = mtl_network,
#    points = bike_accidents,
#    points_time = bike_accidents$Time,
#    start_net = 0,
#    end_net = 1500,
#    step_net = 50,
#    width_net = 200,
#    start_time = 0,
#    end_time = max(bike_accidents$Time),
#    step_time = 7,
#    width_time = 14,
#    nsim = 500,
#    conf_int = 0.05,
#    digits = 2,
#    tol = 0.01,
#    resolution = NULL,
#    agg = NULL,
#    verbose = TRUE,
#    calc_g_func = TRUE
#  )

## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, out.width='100%', fig.dpi = 100----

df_plot <- lapply(1:nrow(net_time_values$obs_k), function(i){
  
  df <- data.frame(
    k_obs = net_time_values$obs_k[i,],
    k_low = net_time_values$lower_k[i,],
    k_upp = net_time_values$upper_k[i,],
    net_dist = net_time_values$distances_net[[i]],
    time_dist = net_time_values$distances_time
  )
  return(df)
  
})

df_plot <- do.call(rbind, df_plot)

df_plot$sign <- ifelse(df_plot$k_obs > df_plot$k_upp, 'concentrated', 'random')
df_plot$sign <- ifelse(df_plot$k_obs < df_plot$k_low, 'dispersed', df_plot$sign)

ggplot() + 
  geom_raster(
    data = subset(df_plot, df_plot$sign != 'random'),
    mapping = aes(x = time_dist, y = net_dist, fill = k_obs  )) + 
  geom_raster(
    data = subset(df_plot, df_plot$sign == 'random'),
    mapping = aes(x = time_dist, y = net_dist), fill = 'grey')



## ----fig.align="center", fig.show='hold', message=FALSE, warning=FALSE, out.width='100%', fig.dpi = 100----

df_plot <- subset(df_plot, df_plot$net_dist > 0 & df_plot$time_dist > 0)
df_plot$diff <- (df_plot$k_obs - df_plot$k_upp) / df_plot$k_obs

ggplot() + 
  geom_raster(
    data = subset(df_plot, df_plot$sign != 'random'),
    mapping = aes(x = time_dist, y = net_dist, fill = diff)) + 
  geom_raster(
    data = subset(df_plot, df_plot$sign == 'random'),
    mapping = aes(x = time_dist, y = net_dist), fill = 'grey')



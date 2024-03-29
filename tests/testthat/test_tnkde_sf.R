context("testing the kernel functions")
library(sf)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Testing the kernels with a SPARSE matrix ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## -- TEST FOR THE SIMPLE TNKDE -- ##

test_that("Testing the simple tnkde with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1),
                      id = c(1))
  event <- st_as_sf(event, coords = c("x","y"))
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                      y=c(0.1), id = 1)
  sp_point <- st_as_sf(sp_point, coords = c("x","y"))
  sample_times <- c(1,2,3,4,5)


  # real distance is 6, and let us say that the bw_net is 10 and bw_times is 5
  net_density <- (1/10) * quartic_kernel(6,10)
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

  #let us calculate the value with our function
  obs_value <- tnkde(
       lines = all_lines,
       events = event,
       w = c(1),
       time_field = "time",
       samples_loc = sp_point,
       samples_time = sample_times,
       check = F,
       kernel_name = "quartic",
       bw_net  = 10,
       bw_time = 5,
       adaptive = F,
       method = "simple",
       div = "bw",
       agg = 0.01,
       verbose = TRUE,
       tol = 0.0001,
       digits = 2
       )
  comp <- sum(round(abs(obs_value - result_mat),9))
  expect_equal(comp, 0)
})


## -- TEST FOR THE DISCONTINUOUS TNKDE -- ##

test_that("Testing the discontinuous tnkde with a simple case", {
  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 5.1, 0.1 0.1)",
    "LINESTRING (-5.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -5.1, 0.1 0.1)",
    "LINESTRING (5.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  event <- st_as_sf(event, coords = c("x","y"))
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(3.1),
                         y=c(0.1))
  sp_point <- st_as_sf(sp_point, coords = c("x","y"))
  sample_times <- c(1,2,3,4,5)


  # real distance is 6, and let us say that the bw_net is 10 and bw_times is 5
  net_density <- (1/10) * quartic_kernel(6,10) * 1/3
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

  #let us calculate the value with our function
  obs_value <- tnkde(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 10,
    bw_time = 5,
    adaptive = F,
    method = "discontinuous",
    div = "bw",
    agg = 0.01,
    verbose = F,
    tol = 0.0001,
    digits = 2
  )

  obs_value2 <- tnkde(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 10,
    bw_time = 5,
    adaptive = F,
    method = "discontinuous",
    div = "bw",
    agg = 0.01,
    verbose = TRUE,
    tol = 0.0001,
    digits = 2, sparse = FALSE
  )

  test1 <- sum(abs(obs_value2 - obs_value)) == 0
  test2 <- sum(round(abs(obs_value - result_mat),9)) == 0

  expect_true(test1 & test2)
})



## -- TEST FOR THE CONTINUOUS TNKDE -- ##

test_that("Testing the continuous tnkde with a simple case", {

  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 10.1, 0.1 0.1)",
    "LINESTRING (-10.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -10.1, 0.1 0.1)",
    "LINESTRING (10.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  event <- st_as_sf(event, coords = c("x","y"))
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(0.1),
                         y=c(1.1))
  sp_point <- st_as_sf(sp_point, coords = c("x","y"))
  sample_times <- c(1,2,3,4,5)

  # real distance is 2, and let us say that the network bw is 5 and bw_times is 5
  net_density <- (1/5) * (quartic_kernel(2,5) - ((1/2) * quartic_kernel(4,5)))
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

  #let us calculate the value with our function
  obs_value <- tnkde(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 5,
    bw_time = 5,
    adaptive = F,
    method = "continuous",
    div = "bw",
    agg = 0.01,
    verbose = F,
    tol = 0.0001,
    digits = 2
  )

  obs_value2 <- tnkde(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 5,
    bw_time = 5,
    adaptive = F,
    method = "continuous",
    div = "bw",
    agg = 0.01,
    verbose = TRUE,
    sparse = FALSE,
    tol = 0.0001,
    digits = 2
  )


  test1 <- sum(abs(obs_value2 - obs_value)) == 0
  test2 <- sum(round(abs(obs_value - result_mat),9)) == 0

  expect_true(test1 & test2)

})



## -- TEST FOR THE CONTINUOUS TNKDE -- ##

test_that("Testing the multicore version of the function TNKDE", {

  ## creating the simple situation
  # start with de definition of some lines
  wkt_lines <- c(
    "LINESTRING (0.1 10.1, 0.1 0.1)",
    "LINESTRING (-10.1 0.1, 0.1 0.1)",
    "LINESTRING (0.1 -10.1, 0.1 0.1)",
    "LINESTRING (10.1 0.1, 0.1 0.1)")

  linesdf <- data.frame(wkt = wkt_lines,
                        id = paste("l",1:length(wkt_lines),sep=""))

  all_lines <- st_as_sf(linesdf, wkt = "wkt")

  # definition of one event
  event <- data.frame(x=c(0.1),
                      y=c(3.1))
  event <- st_as_sf(event, coords = c("x","y"))
  event$time <- 5

  # definition of one sampling point
  sp_point <- data.frame(x=c(0.1),
                         y=c(1.1))
  sp_point <- st_as_sf(sp_point, coords = c("x","y"))

  sample_times <- c(1,2,3,4,5)

  # real distance is 2, and let us say that the network bw is 5 and bw_times is 5
  net_density <- (1/5) * (quartic_kernel(2,5) - ((1/2) * quartic_kernel(4,5)))
  time_densities <- (1/5) * quartic_kernel(abs(sample_times-event$time),5)

  result_mat <-sapply(time_densities, function(x){
    x*net_density
  })

  #let us calculate the value with our function
  obs_value <- tnkde(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 5,
    bw_time = 5,
    adaptive = F,
    method = "continuous",
    div = "bw",
    agg = 0.01,
    verbose = TRUE,
    tol = 0.0001,
    digits = 2
  )

  future::plan(future::multisession(workers=1))
  obs_value2 <- tnkde.mc(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 5,
    bw_time = 5,
    adaptive = F,
    method = "continuous",
    div = "bw",
    agg = 0.01,
    verbose = TRUE,
    tol = 0.0001,
    digits = 2, grid_shape = c(3,3)
  )

  future::plan(future::multisession(workers=1))
  obs_value3 <- tnkde.mc(
    lines = all_lines,
    events = event,
    w = c(1),
    time_field = "time",
    samples_loc = sp_point,
    samples_time = sample_times,
    check = F,
    kernel_name = "quartic",
    bw_net  = 5,
    bw_time = 5,
    adaptive = F,
    method = "continuous",
    div = "bw",
    agg = 0.01,
    verbose = FALSE, sparse = FALSE,
    tol = 0.0001,
    digits = 2, grid_shape = c(3,3)
  )

  test0 <- sum(abs(obs_value2 - obs_value3)) == 0
  test1 <- sum(abs(obs_value2 - obs_value)) == 0
  test2 <- sum(round(abs(obs_value - result_mat),9)) == 0
  expect_true(test1 & test2 & test0)

})




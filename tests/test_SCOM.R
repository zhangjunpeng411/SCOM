library(testthat)
library(SCOM)
library(igraph)

# Network similarity
net1 <- list(as_data_frame(sample_k_regular(10, 2)), as_data_frame(sample_k_regular(10, 2))) 
net2 <- list(as_data_frame(sample_k_regular(10, 3)), as_data_frame(sample_k_regular(10, 3))) 
Sim.net <- Sim.ceRNet(net1, net2)

test_that("Test Sim.ceRNet", {
  expect_equal(Sim.ceRNet(net1, net2), Sim.net)
})

# Hub similarity
hub1 <- list(c("ncRNA1", "ncRNA2", "ncRNA3", "ncRNA4"), c("ncRNA1", "ncRNA2", "ncRNA4", "ncRNA5")) 
hub2 <- list(c("ncRNA1", "ncRNA2", "ncRNA5", "ncRNA6"), c("ncRNA1", "ncRNA3", "ncRNA4", "ncRNA6"))
Sim_hub <- Sim.hub(hub1, hub2) 

test_that("Test Sim.hub", {
    expect_equal(Sim.hub(hub1, hub2), Sim_hub)
})

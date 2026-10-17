#### Sampling error in online social network data -- a study exploring the potential of data donations ####

## author: Anonymized for review

## notice: as we cannot share the data used for this study due to legal reasons this is not a reproducible script; any personally identifiable information was removed;
##         for the sake of sparsity the code is not fully repeated if only individual parameters were changed e.g. the number of layers to sample or the bias in seed selection,
##         in this case see remarks in comments in code.

## The following files are used in this script
# graph_full_inner - igraph object, this is the giant component of the Twitter follower network as of Sept 2021 containing 911,331 vertices and vertex attributes gender (3 categories), party (3 categories) and ethnicity (6 categories)
# covid_states_donation - data.frame, contains all eligible data donations from the Covid States Study
# covid_states_donation_user_details - data.frame, details on donated accounts

# filenames have been anonymized

#### Table of Content ####

## Step 1: calculation of the full network statistics and attribute distributions
## Step 2: Initialize donated accounts as seeds
## Step 3: BFS samples with donated seeds
## Step 4: BFS samples with random seeds
## Step 5: BFS with gender bias in seed selection
### Step 5a: oversampling females
### Step 5b: Oversampling males
## Step 6: BFS with party bias in seed selection
### Step 6a: oversampling democrates
### Step 6b: oversampling republicans
## Step 7: BFS with ethnicity/race bias in seed selection
### Step 7a: oversampling Caucasians
### Step 7b: oversampling Asian
## Step 8: BFS with popularity bias in seed selection
## Step 9: Random Walk sampling with j = 0.15
## Step 10: calculate attribute representations for all samples
## Step 11: Reweighting attributes
## Step 12: Reweighting statistics: density, local clc and ethnic assortativity
## Step 13: Reweighting statistics: global clc
## Step 14: NMSE unweighted estimates
## Step 15: NMSE reweighted estimates
## Step 16: Indegree inclusion share correlation for Appendix
## Step 17: Plots main paper: Attributes
## Step 18: Plots main paper: Statistics
## Step 19: Tables for Appendix: Descriptives
## Step 20: Tables for Appendix: Attributes 
## Step 21: Tables for Appendix: Network Statistics
## Step 21: Tables for Appendix: NMSE


##########################################################################################
##########################################################################################
##########################################################################################

library(tidyverse)
library(igraph)
library(stringr)
library(parallel)
#for plots:
library(gridExtra)
library(grid)
library(cowplot)
library(ggh4x)
library(ggnewscale)
library(patchwork)
#for tables
library(knitr)
library(kableExtra)
library(stringr)


options(arrow.skip_nul = TRUE,
        scipen = 999)

## Step 1: calculation of the full network statistics and attribute distributions

#group ethnic
V(graph_full_inner)$group <- NA
V(graph_full_inner)$group[V(graph_full_inner)$race == "Caucasian"] <- 1
V(graph_full_inner)$group[V(graph_full_inner)$race == "Asian"]     <- 2
V(graph_full_inner)$group[V(graph_full_inner)$race == "African-American"] <- 3
V(graph_full_inner)$group[V(graph_full_inner)$race == "Hispanic"]     <- 4
V(graph_full_inner)$group[V(graph_full_inner)$race == "Other"] <- 5
V(graph_full_inner)$group[V(graph_full_inner)$race == "Unknown"]     <- 6

#group gender
V(graph_full_inner)$group_gender <- NA
V(graph_full_inner)$group_gender[V(graph_full_inner)$gender == "Female"] <- 1
V(graph_full_inner)$group_gender[V(graph_full_inner)$gender == "Male"] <- 2
V(graph_full_inner)$group_gender[V(graph_full_inner)$gender == "Unknown"] <- 3

#group party affiliation
V(graph_full_inner)$group_party <- NA
V(graph_full_inner)$group_party[V(graph_full_inner)$party == "Democrat"] <- 1
V(graph_full_inner)$group_party[V(graph_full_inner)$party == "Republican"] <- 2
V(graph_full_inner)$group_party[V(graph_full_inner)$party == "Independent"] <- 3

compute_network_summary_whole <- data.frame(
  sample_number = "whole network",
  num_nodes = vcount(graph_full_inner),
  num_edges = ecount(graph_full_inner),
  avg_indegree = mean(degree(graph_full_inner, mode = "in")),
  avg_outdegree = mean(degree(graph_full_inner, mode = "out")),
  density = edge_density(graph_full_inner, loops = FALSE),
  clustering_coeff = transitivity(graph_full_inner, type = "localaverage"),
  num_components = components(graph_full_inner)$no, # 1 comp. is correct
  largest_component_size = max(components(graph_full_inner)$csize), 
  assortativity_ethnic = assortativity_nominal(graph_full_inner, V(graph_full_inner)$group),
  assortativity_gender = assortativity_nominal(graph_full_inner, V(graph_full_inner)$group_gender),
  assortativity_party = assortativity_nominal(graph_full_inner, V(graph_full_inner)$group_party),
  clustering_coeff_global_vector= transitivity(graph_full_inner, type = "global")
)

#compute_network_summary_whole
write.csv(compute_network_summary_whole, file="whole_network_descriptives_in_voterfile.csv")

attributes_full_df <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(attributes_full_df) <- c("sample_num","sample_method","sample_depth","sample_size","female","male","gender_unknown",
                                  "caucasian","african_american","asian","hispanic","race_other","race_unknown",
                                  "democrat","republican","independent")
#get statistics 
sample_size <- vcount(graph_full_inner)
female <- sum(V(graph_full_inner)$gender == "Female")
male <- sum(V(graph_full_inner)$gender == "Male")
gender_unknown <- sum(V(graph_full_inner)$gender == "Unknown")
caucasian <- sum(V(graph_full_inner)$race == "Caucasian")
african_american <- sum(V(graph_full_inner)$race == "African-American")
asian <- sum(V(graph_full_inner)$race == "Asian")
hispanic <- sum(V(graph_full_inner)$race == "Hispanic")
race_other <- sum(V(graph_full_inner)$race == "Other")
race_unknown <- sum(V(graph_full_inner)$race == "Unknown")
democrat <- sum(V(graph_full_inner)$party == "Democrat")
republican <- sum(V(graph_full_inner)$party == "Republican")
independent <- sum(V(graph_full_inner)$party == "Independent")

attributes_full_df <- rbind(attributes_full_df, c(0,"full","full",sample_size,female,male,gender_unknown,caucasian,african_american,asian,
                                                  hispanic,race_other,race_unknown,democrat,republican,independent))

colnames(attributes_full_df) <- c("sample_num","sample_method","sample_depth","sample_size","female","male","gender_unknown",
                                  "caucasian","african_american","asian","hispanic","race_other","race_unknown",
                                  "democrat","republican","independent")

write.csv(attributes_full_df, file="whole_network_attributes.csv")


## Step 2: Initialize donated accounts as seeds
mat_inner_senders <- V(graph_full_inner)$name # account list from full network considered

covid_states_overlap <- covid_states_donation %>% 
  mutate(in_inner_matrix = if_else(user_id %in% mat_inner_senders,1,0),
         user_id_str = as.character(user_id))

covid_states_overlap_inner <- covid_states_overlap %>% filter(in_inner_matrix == 1)

#find number of followers for donnated handles
accounts_donnated <- covid_states_overlap %>% filter(in_inner_matrix == 1) %>% select(uid = user_id)
user_details_overlap <- covid_states_overlap_inner %>% left_join(., covid_states_donation_user_details, by = "user_id_str")

write.csv(user_details_overlap, "overlap_covidstates_ids_inner_graph.csv")

## Step 3: BFS samples with donated seeds

seed_data <- user_details_overlap # read.csv("overlap_covidstates_ids_inner_graph.csv")
seed_nodes <- intersect(V(graph_full_inner)$userid, seed_data$user_id) 
length(seed_nodes) #564

V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Functions

bfs_sample_single <- function(graph, seed_user_ids, k) {
  library(igraph)
  #seed_nodes <- seed_nodes
  seed_nodes <- match(seed_user_ids, V(graph)$userid)
  sampled_nodes <- seed_nodes
  
  for (seed in seed_nodes) {
    reached_nodes <- unlist(neighborhood(graph,k,seed, mode = "out")) 
    sampled_nodes <- unique(c(sampled_nodes, reached_nodes))
  }
  
  subgraph <- induced_subgraph(graph, sampled_nodes)
  
  return(list(seed_nodes = seed_nodes, sampled_nodes = sampled_nodes, subgraph = subgraph))
}


#Initialization
N = length(seed_nodes) # number of seed nodes per sample
k = 1 # layers to be sampled 1 and 2
S = 1 # number of samples to draw - in case of donation always 1

#Sampling
results_bfs <- bfs_sample_single(graph_full_inner, seed_user_ids=seed_nodes, k = k)


#Remove edges in last layer
edge_in_layst_layer_remover_single <- function(k,result, user_id_graph_full_inner){
  library(igraph)
  subgraph <- result$subgraph
  seed_user_ids <- user_id_graph_full_inner[result$seed_nodes]
  
  #seeds
  V(subgraph)$is_seed <- V(subgraph)$userid %in% seed_user_ids
  table(V(subgraph)$is_seed)
  
  #layer 1: neigbors of seeds
  neighbor_verts <- unlist(neighborhood(subgraph, order = 1, nodes = seed_user_ids, mode = "out"))
  V(subgraph)$is_layer1 <- FALSE
  V(subgraph)$is_layer1[neighbor_verts] <- TRUE
  V(subgraph)$is_layer1[V(subgraph)$is_seed == TRUE] <- FALSE
  table(V(subgraph)$is_seed, V(subgraph)$is_layer1)
  
  #layer2: neighbors of neigbors of seeds
  V(subgraph)$is_layer2 <- FALSE
  V(subgraph)$is_layer2[V(subgraph)$is_seed == FALSE & V(subgraph)$is_layer1 == FALSE] <- TRUE
  table(V(subgraph)$is_layer1, V(subgraph)$is_layer2)
  
  if (k == 1) {
    
    last_layer <- which(V(subgraph)$is_layer1)
    
  } else { if (k == 2) {
    last_layer <- which(V(subgraph)$is_layer2)
    
  } else {cat("Error: k is neither 1 or 2")}}
  
  #identify edges within last layer
  edge_ends <- ends(subgraph, E(subgraph), names = FALSE)
  
  E_last_last <- E(subgraph)[
    edge_ends[,1] %in% last_layer &
      edge_ends[,2] %in% last_layer
  ]
  
  #remove edges within last layer
  subgraph_cleaned <- delete_edges(subgraph, E_last_last)
  
  #share edges preseved
  share_preserved  <- ecount(subgraph_cleaned) / ecount(subgraph)
  
  return(list(sample_number = 1, seed_nodes = result$seed_nodes, sampled_nodes = result$sampled_nodes, subgraph = subgraph_cleaned, 
              edge_share_preserved = share_preserved))
}

#apply to donated samples: this is repeated for k1 and k2
user_id_graph_full_inner <- V(graph_full_inner)$userid
results_bfs <- edge_in_layst_layer_remover_single(k,result,user_id_graph_full_inner)


#Summary statistics
compute_network_summary_single <- function(results) {
  graph <- results$subgraph
  num_nodes <- vcount(graph)
  num_edges <- ecount(graph)
  comp <- components(graph)
  
  summary_df <- data.frame(
    sample_number = 1,
    num_nodes = num_nodes,
    num_edges = num_edges,
    avg_degree = mean(degree(graph, mode = "in")),
    max_outdegree = max(degree(graph, mode = "out")),
    min_outdegree = min(degree(graph, mode = "out")),
    max_indegree = max(degree(graph, mode = "in")),
    min_indegree = min(degree(graph, mode = "in")),
    density = edge_density(graph, loops = FALSE),
    clustering_coeff = transitivity(graph, type = "localaverage"),
    num_components = comp$no,
    largest_component_size = max(comp$csize),
    share_largest_comp = max(comp$csize) / 911331,
    num_seeds = length(results$seed_nodes)
  )
  
  return(summary_df)
}

network_descriptives_df <- compute_network_summary_single(results_bfs)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("S",S,"N",N,"k",k,"network_sample_descriptives_BFS",".csv"))
#save results
save(results_bfs,network_descriptives_df, file=paste0("donated","_S",S,"N",N,"k",k,"network_sample_results_BFS",".rdata"))


## Step 4: BFS samples with random seeds

V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)) #save original IDs before sampling creates new IDs

#Functions

bfs_sample <- function(graph, N, k, sample_number) {
  library(igraph)
  seed_nodes <- as.numeric(sample(V(graph), N, replace = FALSE))
  sampled_nodes <- seed_nodes
  
  for (seed in seed_nodes) {
    reached_nodes <- unlist(neighborhood(graph,k,seed, mode = "out")) 
    sampled_nodes <- unique(c(sampled_nodes, reached_nodes))
  }
  
  subgraph <- induced_subgraph(graph, sampled_nodes)
  
  return(list(sample_number = sample_number, seed_nodes = seed_nodes, sampled_nodes = sampled_nodes, subgraph = subgraph))
}

#Function to perform multiple BFS samples in parallel
bfs_multiple_samples <- function(graph, N, k, S, num_cores) {
  cl <- makeCluster(num_cores)
  clusterExport(cl, list("graph", "N", "k", "bfs_sample"), envir = environment())
  results <- parLapply(cl, 1:S, function(s) bfs_sample(graph, N, k, s))
  stopCluster(cl)
  return(results)
}

#Initialization
run <- 1 # 1 to 10 - sampling is split up into 10 runs because 100 at a time is too memory intense
N = 564 # number of seed nodes per sample
k = 1 # layers to be sampled: 1 and 2
S = 10 # number of samples to draw
C = 10 #number of cores to use # note: depending on the resources available on the shared server this was varied between 3 and 10

results_bfs <- bfs_multiple_samples(graph_full_inner, N = N, k = k, S = S, num_cores = C)

user_id_graph_full_inner <- V(graph_full_inner)$userid

edge_in_layst_layer_remover <- function(k,result, user_id_graph_full_inner, sample_number){
  library(igraph)
  subgraph <- result[[sample_number]]$subgraph
  seed_user_ids <- user_id_graph_full_inner[result[[sample_number]]$seed_nodes]
  
  #seeds
  V(subgraph)$is_seed <- V(subgraph)$userid %in% seed_user_ids
  table(V(subgraph)$is_seed)
  
  #layer 1: neigbors of seeds
  neighbor_verts <- unlist(neighborhood(subgraph, order = 1, nodes = seed_user_ids, mode = "out"))
  V(subgraph)$is_layer1 <- FALSE
  V(subgraph)$is_layer1[neighbor_verts] <- TRUE
  V(subgraph)$is_layer1[V(subgraph)$is_seed == TRUE] <- FALSE
  table(V(subgraph)$is_seed, V(subgraph)$is_layer1)
  
  #layer2: neighbors of neigbors of seeds
  V(subgraph)$is_layer2 <- FALSE
  V(subgraph)$is_layer2[V(subgraph)$is_seed == FALSE & V(subgraph)$is_layer1 == FALSE] <- TRUE
  table(V(subgraph)$is_layer1, V(subgraph)$is_layer2)
  
  if (k == 1) {
    
    last_layer <- which(V(subgraph)$is_layer1)
    
  } else { if (k == 2) {
    last_layer <- which(V(subgraph)$is_layer2)
    
  } else {cat("Error: k is neither 1 or 2")}}
  
  #identify edges within last layer
  
  edge_ends <- ends(subgraph, E(subgraph), names = FALSE)
  
  E_last_last <- E(subgraph)[
    edge_ends[,1] %in% last_layer &
      edge_ends[,2] %in% last_layer
  ]
  
  #remove edges within last layer
  subgraph_cleaned <- delete_edges(subgraph, E_last_last)
  
  #share edges preseved
  share_preserved  <- ecount(subgraph_cleaned) / ecount(subgraph)
  
  return(list(sample_number = sample_number, seed_nodes = result[[sample_number]]$seed_nodes, sampled_nodes = result[[sample_number]]$sampled_nodes, subgraph = subgraph_cleaned, 
              edge_share_preserved = share_preserved))
}

# for all samples
edge_in_layst_layer_remover_multiple_samples <- function(result, user_id_graph_full_inner, k, S, num_cores) {
  cl <- makeCluster(num_cores)
  clusterExport(cl, list("result", "user_id_graph_full_inner", "k", "edge_in_layst_layer_remover"), envir = environment())
  results <- parLapply(cl, 1:S, function(s) edge_in_layst_layer_remover(k, result, user_id_graph_full_inner, s))
  stopCluster(cl)
  return(results)
}

# apply to samples: this is repeated for all runs of k1 and k2
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)


#Summary statistics
compute_network_summary <- function(results, num_cores = C) {
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, library(igraph))  #igraph needs to be reloaded in workers of cluster
  clusterExport(cl, list("results", "vcount", "ecount", "degree", "edge_density", "components", "transitivity"), envir = environment())
  
  summary_list <- parLapply(cl, results, function(res) {
    graph <- res$subgraph
    num_nodes <- vcount(graph)
    num_edges <- ecount(graph)
    comp <- components(graph)
    
    data.frame(
      sample_number = res$sample_number,
      num_nodes = num_nodes,
      num_edges = num_edges,
      avg_indegree = mean(degree(graph, mode = "in")),
      avg_outdegree = mean(degree(graph, mode = "out")),
      density = edge_density(graph, loops = FALSE),
      clustering_coeff = transitivity(graph, type = "localaverage"),  
      num_components = comp$no,
      largest_component_size = max(comp$csize),
      share_largest_comp = max(comp$csize) / 911331
    )
  })
  
  stopCluster(cl)
  return(rbindlist(summary_list))
}

network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("BFS_unbiased/","S",S,"N",N,"k",k,"network_sample_descriptives_BFS_run",run,".csv"))
#save results
save(results_bfs,network_descriptives_df, file=paste0("BFS_unbiased/","_S",S,"N",N,"k",k,"network_sample_results_BFS_run",run,".rdata"))


##Step 5: BFS with gender bias in seed selection

## Step 5a: oversampling females

table(V(graph_full_inner)$gender) # absolute num vertices by gender

table(V(graph_full_inner)$gender) / vcount(graph_full_inner) #percentage by gender
table(V(graph_full_inner)$gender) / vcount(graph_full_inner) * 564 #expected random seeds by gender; M = 564

#keep sampling probability for the gender == unknown constant and skew the number of male/female nodes
#calculate how many seeds to sample for each gender
expected_sample_counts <- function(
    n_female, n_male, n_neutral,
    n_total,
    female_bias = 1.05,  # e.g., 1.05 means females are 5% more likely # other values used: 1.10 for 10% more and 2 for 100% more
    sample_size = 564
){
  p_neutral <- n_neutral / n_total
  M_neutral <- round(sample_size*p_neutral)
  sample_size_remain <- sample_size - M_neutral
  
  weighted_total = n_female*female_bias+n_male
  p <- sample_size_remain / weighted_total
  
  M_female <- round(n_female*female_bias*p)
  M_male <- round(n_male*p)
  
  if (sum(M_female,M_male,M_neutral) == sample_size) {
    return(c(M_female,M_male,M_neutral))
  } else {
    if (sum(M_female,M_male,M_neutral) > sample_size) {
      M_neutral <- M_neutral -1
    } else {
      M_neutral <- M_neutral +1
    }
    return(c(M_female,M_male,M_neutral))
  }
}

#calculate number of nodes to sample
M_by_gender <- expected_sample_counts(sum(V(graph_full_inner)$gender == "Female"),
                                      sum(V(graph_full_inner)$gender == "Male"),
                                      sum(V(graph_full_inner)$gender == "Unknown"), vcount(graph_full_inner))


#now sample seeds by gender
# Filter female and male vertices
female_nodes <- V(graph_full_inner)[gender == "Female"]
male_nodes   <- V(graph_full_inner)[gender == "Male"]
unknown_nodes   <- V(graph_full_inner)[gender == "Unknown"]

# Sample
non_random_seed_sampling <- function(){
  sampled_females <- sample(female_nodes, M_by_gender[1]) #sample by gender
  sampled_males   <- sample(male_nodes, M_by_gender[2])
  sampled_unknown   <- sample(unknown_nodes, M_by_gender[3])
  
  sampled_females <- V(graph_full_inner)$userid[sampled_females] # Extract the 'userid' attribute
  sampled_males <- V(graph_full_inner)$userid[sampled_males]
  sampled_unknown <- V(graph_full_inner)$userid[sampled_unknown]
  
  sampled_nodes <- c(sampled_females,sampled_males,sampled_unknown)
  
  return(sampled_nodes)
  
}

sampling_time <- Sys.time()
# list of samples
seeds_list <- replicate(10, non_random_seed_sampling(), simplify = FALSE)
length(seeds_list[[1]]) #length of each of teh 100 sampled vectors

V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Functions

# Function to perform breadth-first sampling
bfs_sample <- function(graph, seed_user_ids, k, sample_number) {
  library(igraph)
  seed_nodes <- match(seed_user_ids, V(graph)$userid)
  sampled_nodes <- seed_nodes
  
  for (seed in seed_nodes) {
    reached_nodes <- unlist(neighborhood(graph,k,seed, mode = "out"))  # Stores nodes within k layers
    sampled_nodes <- unique(c(sampled_nodes, reached_nodes))
  }
  
  subgraph <- induced_subgraph(graph, sampled_nodes)
  
  return(list(sample_number = sample_number, seed_nodes = seed_nodes, sampled_nodes = sampled_nodes, subgraph = subgraph))
}


# Function to perform multiple BFS samples in parallel with specific non_random seeds
bfs_multiple_samples_bias <- function(graph, k, seed_user_ids_list, num_cores) {
  cl <- makeCluster(num_cores)
  clusterExport(cl, list("graph", "k", "bfs_sample"), envir = environment())
  results <- parLapply(cl, seq_along(seed_user_ids_list), function(i) {
    seeds <- seed_user_ids_list[[i]]
    bfs_sample(graph, seeds, k, sample_number = i)
  })
  stopCluster(cl)
  return(results)
}

#Initialization
k = 1 # layers to be sampled: 1 or 2
#S not needed as defined by length of seed_user_ids_list, always 10 in the application
C = 10 #number of cores to use
run <- 1

#apply 
results_bfs <- bfs_multiple_samples_bias(graph_full_inner, k = k, seed_user_ids_list = seeds_list, num_cores = C)
S = length(results_bfs)
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)

network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("women5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv")) # equivalently women10percentmore. women100percentmore
#save results
save(results_bfs,network_descriptives_df, file=paste0("women5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata")) # equivalently women10percentmore. women100percentmore

## Step 5b: Oversampling males

expected_sample_counts <- function(
    n_female, n_male, n_neutral,
    n_total,
    male_bias = 1.05,  # e.g., 1.05 means males are 5% more likely to donate than females; # other values used: 1.10 for 10% more and 2 for 100% more
    sample_size = 564
){
  p_neutral <- n_neutral / n_total
  M_neutral <- round(sample_size*p_neutral)
  sample_size_remain <- sample_size - M_neutral
  
  weighted_total = n_female+n_male*male_bias
  p <- sample_size_remain / weighted_total
  
  M_female <- round(n_female*p)
  M_male <- round(n_male*p*male_bias)
  
  if (sum(M_female,M_male,M_neutral) == sample_size) {
    return(c(M_female,M_male,M_neutral))
  } else {
    if (sum(M_female,M_male,M_neutral) > sample_size) {
      M_neutral <- M_neutral -1
    } else {
      M_neutral <- M_neutral +1
    }
    return(c(M_female,M_male,M_neutral))
  }
}

#calculate number of nodes to sample
M_by_gender <- expected_sample_counts(sum(V(graph_full_inner)$gender == "Female"),
                                      sum(V(graph_full_inner)$gender == "Male"),
                                      sum(V(graph_full_inner)$gender == "Unknown"), vcount(graph_full_inner))


#now sample seeds by gender
# Filter female and male vertices
female_nodes <- V(graph_full_inner)[gender == "Female"]
male_nodes   <- V(graph_full_inner)[gender == "Male"]
unknown_nodes   <- V(graph_full_inner)[gender == "Unknown"]

# Sample
seeds_list <- replicate(10, non_random_seed_sampling(), simplify = FALSE)
V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Initialization
k = 2 # layers to be sampled 1 or 2
#S not needed as defined by length of seed_user_ids_list, always 10 in the application
C = 10
run <- 1

#apply
results_bfs <- bfs_multiple_samples_bias(graph_full_inner, k = k, seed_user_ids_list = seeds_list, num_cores = C)
S = length(results_bfs)
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)

network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("men5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv")) # equivalently men10percentmore. men100percentmore
#save results
save(results_bfs,network_descriptives_df, file=paste0("men5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata")) # equivalently men10percentmore. men100percentmore


##Step 6: BFS with party bias in seed selection

## Step 6a: oversampling democrates
table(V(graph_full_inner)$party) # absolute num vertices by party

table(V(graph_full_inner)$party) / vcount(graph_full_inner) #percentage by party
table(V(graph_full_inner)$party) / vcount(graph_full_inner) * 564 #expected random seeds by party; M = 564

#keep sampling probability for independent constant and skew the number of dem/rep nodes
#calculate how many seeds to sample for each party
expected_sample_counts <- function(
    n_dem, n_rep, n_indep,
    n_total,
    democrat_bias = 1.05,  # e.g., 1.05 means democrats are 5% more likely to donate than republicans; # other values used: 1.10 for 10% more and 2 for 100% more
    sample_size = 564
){
  p_indep <- n_indep / n_total
  M_indep <- round(sample_size*p_indep)
  sample_size_remain <- sample_size - M_indep
  
  weighted_total = n_dem*democrat_bias+n_rep
  p <- sample_size_remain / weighted_total
  
  M_dem <- round(n_dem*democrat_bias*p)
  M_rep <- round(n_rep*p)
  
  if (sum(M_dem,M_rep,M_indep) == sample_size) {
    return(c(M_dem,M_rep,M_indep))
  } else {
    if (sum(M_dem,M_rep,M_indep) > sample_size) {
      M_indep <- M_indep -1
    } else {
      M_indep <- M_indep +1
    }
    return(c(M_dem,M_rep,M_indep))
  }
}

#calculate number of nodes to sample
M_by_party <- expected_sample_counts(sum(V(graph_full_inner)$party == "Democrat"),
                                     sum(V(graph_full_inner)$party == "Republican"),
                                     sum(V(graph_full_inner)$party == "Independent"), vcount(graph_full_inner))

dem_nodes <- V(graph_full_inner)[party == "Democrat"]
rep_nodes   <- V(graph_full_inner)[party == "Republican"]
indep_nodes   <- V(graph_full_inner)[party == "Independent"]

# Sample
non_random_seed_sampling <- function(){
  sampled_dem <- sample(dem_nodes, M_by_party[1]) #sample by party
  sampled_rep   <- sample(rep_nodes, M_by_party[2])
  sampled_indep   <- sample(indep_nodes, M_by_party[3])
  
  sampled_dem <- V(graph_full_inner)$userid[sampled_dem] # Extract the 'userid' attribute
  sampled_rep <- V(graph_full_inner)$userid[sampled_rep]
  sampled_indep <- V(graph_full_inner)$userid[sampled_indep]
  
  sampled_nodes <- c(sampled_dem,sampled_rep,sampled_indep)
  
  return(sampled_nodes)
  
}

sampling_time <- Sys.time()
# list of samples 
seeds_list <- replicate(10, non_random_seed_sampling(), simplify = FALSE)
V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Initialization
k = 1 # layers to be sampled: 1 or 2
C = 10 #number of cores to use
run <- 1

#apply
results_bfs <- bfs_multiple_samples_bias(graph_full_inner, k = k, seed_user_ids_list = seeds_list, num_cores = C)
S = length(results_bfs)
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)

network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("democrat05percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv")) # equivalently democrat10percentmore, democrat100percentmore
#save results
save(results_bfs,network_descriptives_df, file=paste0("democrat05percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata")) # equivalently democrat10percentmore, democrat100percentmore


## Step 6b: oversampling republicans

expected_sample_counts <- function(
    n_dem, n_rep, n_indep,
    n_total,
    republican_bias = 1.05,  # e.g., 1.05 means republicans are 5% more likely to donate than democrats; # other values used: 1.10 for 10% more and 2 for 100% more
    sample_size = 564
){
  p_indep <- n_indep / n_total
  M_indep <- round(sample_size*p_indep)
  sample_size_remain <- sample_size - M_indep
  
  weighted_total = n_dem+n_rep*republican_bias
  p <- sample_size_remain / weighted_total
  
  M_dem <- round(n_dem*p)
  M_rep <- round(n_rep*republican_bias*p)
  
  if (sum(M_dem,M_rep,M_indep) == sample_size) {
    return(c(M_dem,M_rep,M_indep))
  } else {
    if (sum(M_dem,M_rep,M_indep) > sample_size) {
      M_indep <- M_indep -1
    } else {
      M_indep <- M_indep +1
    }
    return(c(M_dem,M_rep,M_indep))
  }
}

#calculate number of nodes to sample: females 5% more likely to donate
M_by_party <- expected_sample_counts(sum(V(graph_full_inner)$party == "Democrat"),
                                     sum(V(graph_full_inner)$party == "Republican"),
                                     sum(V(graph_full_inner)$party == "Independent"), vcount(graph_full_inner))


#now sample seeds by gender
# Filter female and male vertices
dem_nodes <- V(graph_full_inner)[party == "Democrat"]
rep_nodes   <- V(graph_full_inner)[party == "Republican"]
indep_nodes   <- V(graph_full_inner)[party == "Independent"]

# Sample 520 females and 480 males
#set.seed(429)  # for reproducibility
non_random_seed_sampling <- function(){
  sampled_dem <- sample(dem_nodes, M_by_party[1]) #sample by party
  sampled_rep   <- sample(rep_nodes, M_by_party[2])
  sampled_indep   <- sample(indep_nodes, M_by_party[3])
  
  sampled_dem <- V(graph_full_inner)$userid[sampled_dem] # Extract the 'userid' attribute
  sampled_rep <- V(graph_full_inner)$userid[sampled_rep]
  sampled_indep <- V(graph_full_inner)$userid[sampled_indep]
  
  sampled_nodes <- c(sampled_dem,sampled_rep,sampled_indep)
  
  return(sampled_nodes)
  
}

sampling_time <- Sys.time()
# list of samples with females being 5% more likely to donate
seeds_list <- replicate(10, non_random_seed_sampling(), simplify = FALSE)
length(seeds_list[[1]]) #length of each of the 100 sampled vectors

V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Initialization
#N = 564 # number of seed nodes per sample
k = 1 # layers to be sampled: 1 and 2
C = 10 #number of cores to use
run <- 1

#apply
results_bfs <- bfs_multiple_samples_bias(graph_full_inner, k = k, seed_user_ids_list = seeds_list, num_cores = C)
S = length(results_bfs)
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)
network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("republican5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv")) # equivalently republican10percentmore, republican100percentmore
#save results
save(results_bfs,network_descriptives_df, file=paste0("republican5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata")) # equivalently republican10percentmore, republican100percentmore

## Step 7: BFS with ethnicity/race bias in seed selection

## Step 7a: oversampling Caucasians
table(V(graph_full_inner)$race) # absolute num vertices by race

table(V(graph_full_inner)$race) / vcount(graph_full_inner) #percentage by race
table(V(graph_full_inner)$race) / vcount(graph_full_inner) * 564 #expected random seeds by race; M = 564

#keep sampling probability for the other ethnicities identical 
#calculate how many seeds to sample for each ethnicity
expected_sample_counts <- function(
    n_cauc, n_asia, n_afam, n_hisp, n_other, n_unknown,
    n_total,
    caucasian_bias = 1.05,  # e.g., 1.05 means caucasians are 5% more likely to donate than other ethicities ; # other values used: 1.10 for 10% more and 2 for 100% more
    sample_size = 564
){
  sample_size_remain <- sample_size 
  
  weighted_total = n_cauc*caucasian_bias+n_asia+n_afam+n_hisp+n_other+n_unknown
  p <- sample_size_remain / weighted_total
  
  M_cauc <- round(n_cauc*caucasian_bias*p)
  M_asia <- round(n_asia*p)
  M_afam <- round(n_afam*p)
  M_hisp <- round(n_hisp*p)
  M_other <- round(n_other*p)
  M_unknown <- round(n_unknown*p)
  
  while (sum(M_cauc, M_asia, M_afam, M_hisp, M_other, M_unknown) != sample_size) {
    
    if (sum(M_cauc, M_asia, M_afam, M_hisp, M_other, M_unknown) > sample_size) {
      M_unknown <- M_unknown - 1
    } else {
      M_unknown <- M_unknown + 1
    }
    
  } # the number of unknown ethnicity nodes is increased or reduced until there are as many nodes as there should be as seeds
  
  return(c(M_cauc, M_asia, M_afam, M_hisp, M_other, M_unknown))
}

#calculate number of nodes to sample
M_by_ethnicity <- expected_sample_counts(sum(V(graph_full_inner)$race == "Caucasian"),
                                         sum(V(graph_full_inner)$race == "Asian"),
                                         sum(V(graph_full_inner)$race == "African-American"),
                                         sum(V(graph_full_inner)$race == "Hispanic"),
                                         sum(V(graph_full_inner)$race == "Other"),
                                         sum(V(graph_full_inner)$race == "Unknown"), vcount(graph_full_inner))


#now sample seeds by ethnicity
# Filter vertices
cauc_nodes <- V(graph_full_inner)[race == "Caucasian"]
asia_nodes   <- V(graph_full_inner)[race == "Asian"]
afam_nodes <- V(graph_full_inner)[race == "African-American"]
hisp_nodes   <- V(graph_full_inner)[race == "Hispanic"]
other_nodes <- V(graph_full_inner)[race == "Other"]
unknown_nodes   <- V(graph_full_inner)[race == "Unknown"]

non_random_seed_sampling <- function(){
  sampled_cauc <- sample(cauc_nodes, M_by_ethnicity[1]) #sample by ethnicity
  sampled_asia   <- sample(asia_nodes, M_by_ethnicity[2])
  sampled_afam   <- sample(afam_nodes, M_by_ethnicity[3])
  sampled_hisp <- sample(hisp_nodes, M_by_ethnicity[4]) 
  sampled_other   <- sample(other_nodes, M_by_ethnicity[5])
  sampled_unknown   <- sample(unknown_nodes, M_by_ethnicity[6])
  
  sampled_cauc <- V(graph_full_inner)$userid[sampled_cauc] # Extract the 'userid' attribute
  sampled_asia <- V(graph_full_inner)$userid[sampled_asia]
  sampled_afam <- V(graph_full_inner)$userid[sampled_afam]
  sampled_hisp <- V(graph_full_inner)$userid[sampled_hisp]
  sampled_other <- V(graph_full_inner)$userid[sampled_other]
  sampled_unknown <- V(graph_full_inner)$userid[sampled_unknown]
  
  sampled_nodes <- c(sampled_cauc,sampled_asia,sampled_afam,sampled_hisp,sampled_other,sampled_unknown)
  
  return(sampled_nodes)
  
}

seeds_list <- replicate(10, non_random_seed_sampling(), simplify = FALSE)
V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Initialization
k = 1 # layers to be sampled: 1 and 2
C = 10 #number of cores to use
run <- 1

#apply
results_bfs <- bfs_multiple_samples_bias(graph_full_inner, k = k, seed_user_ids_list = seeds_list, num_cores = C)
S = length(results_bfs)
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)
network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("caucasian5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv"))
#save results
save(results_bfs,network_descriptives_df, file=paste0("caucasian5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata"))


## Step 7b: oversampling Asian

#keep sampling probability for the other ethnicities identical 
#calculate how many seeds to sample for each ethnicity
expected_sample_counts <- function(
    n_cauc, n_asia, n_afam, n_hisp, n_other, n_unknown,
    n_total,
    asian_bias = 1.05,  # e.g., 1.05 means asians are 5% more likely to donate than other ethnicities; # other values used: 1.10 for 10% more and 2 for 100% more
    sample_size = 564
){  sample_size_remain <- sample_size 

weighted_total = n_cauc+n_asia*asian_bias+n_afam+n_hisp+n_other+n_unknown
p <- sample_size_remain / weighted_total

M_cauc <- round(n_cauc*p)
M_asia <- round(n_asia*asian_bias*p)
M_afam <- round(n_afam*p)
M_hisp <- round(n_hisp*p)
M_other <- round(n_other*p)
M_unknown <- round(n_unknown*p)

while (sum(M_cauc, M_asia, M_afam, M_hisp, M_other, M_unknown) != sample_size) {
  
  if (sum(M_cauc, M_asia, M_afam, M_hisp, M_other, M_unknown) > sample_size) {
    M_unknown <- M_unknown - 1
  } else {
    M_unknown <- M_unknown + 1
  }
  
} # the number of unknown ethnicity nodes is increased or reduced until there are as many nodes as there should be sampled

return(c(M_cauc, M_asia, M_afam, M_hisp, M_other, M_unknown))
}

#calculate number of nodes to sample
M_by_ethnicity <- expected_sample_counts(sum(V(graph_full_inner)$race == "Caucasian"),
                                         sum(V(graph_full_inner)$race == "Asian"),
                                         sum(V(graph_full_inner)$race == "African-American"),
                                         sum(V(graph_full_inner)$race == "Hispanic"),
                                         sum(V(graph_full_inner)$race == "Other"),
                                         sum(V(graph_full_inner)$race == "Unknown"), vcount(graph_full_inner))


#now sample seeds by ethnicity
# Filter vertices
cauc_nodes <- V(graph_full_inner)[race == "Caucasian"]
asia_nodes   <- V(graph_full_inner)[race == "Asian"]
afam_nodes <- V(graph_full_inner)[race == "African-American"]
hisp_nodes   <- V(graph_full_inner)[race == "Hispanic"]
other_nodes <- V(graph_full_inner)[race == "Other"]
unknown_nodes   <- V(graph_full_inner)[race == "Unknown"]

non_random_seed_sampling <- function(){
  sampled_cauc <- sample(cauc_nodes, M_by_ethnicity[1]) #sample by ethnicity
  sampled_asia   <- sample(asia_nodes, M_by_ethnicity[2])
  sampled_afam   <- sample(afam_nodes, M_by_ethnicity[3])
  sampled_hisp <- sample(hisp_nodes, M_by_ethnicity[4]) 
  sampled_other   <- sample(other_nodes, M_by_ethnicity[5])
  sampled_unknown   <- sample(unknown_nodes, M_by_ethnicity[6])
  
  sampled_cauc <- V(graph_full_inner)$userid[sampled_cauc] # Extract the 'userid' attribute
  sampled_asia <- V(graph_full_inner)$userid[sampled_asia]
  sampled_afam <- V(graph_full_inner)$userid[sampled_afam]
  sampled_hisp <- V(graph_full_inner)$userid[sampled_hisp]
  sampled_other <- V(graph_full_inner)$userid[sampled_other]
  sampled_unknown <- V(graph_full_inner)$userid[sampled_unknown]
  
  sampled_nodes <- c(sampled_cauc,sampled_asia,sampled_afam,sampled_hisp,sampled_other,sampled_unknown)
  
  return(sampled_nodes)
  
}

seeds_list <- replicate(10, non_random_seed_sampling(), simplify = FALSE)
V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)$userid) #save original IDs before sampling creates new IDs

#Initialization
k = 1 # layers to be sampled: 1 and 2
C = 10 #number of cores to use
run <- 1

results_bfs <- bfs_multiple_samples_bias(graph_full_inner, k = k, seed_user_ids_list = seeds_list, num_cores = C)
S = length(results_bfs)
results_bfs <- edge_in_layst_layer_remover_multiple_samples(results_bfs,user_id_graph_full_inner, k, S, C)
network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("asian5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv"))
#save results
save(results_bfs,network_descriptives_df, file=paste0("asian5percentmore/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata"))

## Step 8: BFS with popularity bias in seed selection

#get 95% intervall for degree for bottom and top coding
quantile(degree(graph_full_inner, mode = "in"), probs = c(0.025, 0.975))

#degrees_full_graph <- degree(graph_full_inner, mode = "in")  # Get out-degrees of nodes
V(graph_full_inner)$indegree <- degree(graph_full_inner, mode = "in")
V(graph_full_inner)$indegree <- pmin(pmax(V(graph_full_inner)$indegree, 1), 125)


V(graph_full_inner)$prob_bydegree_norm <- V(graph_full_inner)$indegree /sum(V(graph_full_inner)$indegree)


bfs_sample_biasdegree <- function(graph, N, k, sample_number) {
  library(igraph)
  seed_nodes <- as.numeric(sample(V(graph), N, replace = FALSE, prob = V(graph)$prob_bydegree_norm))
  sampled_nodes <- seed_nodes
  
  for (seed in seed_nodes) {
    reached_nodes <- unlist(neighborhood(graph,k,seed, mode = "out"))  # Stores nodes within k layers
    sampled_nodes <- unique(c(sampled_nodes, reached_nodes))
  }
  
  subgraph <- induced_subgraph(graph, sampled_nodes)
  
  return(list(sample_number = sample_number, seed_nodes = seed_nodes, sampled_nodes = sampled_nodes, subgraph = subgraph))
}

# Function to perform multiple BFS samples in parallel
bfs_multiple_samples_biasdegree <- function(graph, N, k, S, num_cores) {
  cl <- makeCluster(num_cores)
  clusterExport(cl, list("graph", "N", "k", "bfs_sample_biasdegree"), envir = environment())
  results <- parLapply(cl, 1:S, function(s) bfs_sample_biasdegree(graph, N, k, s))
  stopCluster(cl)
  return(results)
}
#Initialization
N = 564 # number of seed nodes per sample
k = 1 # layers to be sampled #2 #3
S = 10 # number of samples to draw
C = 10 #number of cores to use
run <- 1

#apply
results_bfs <- bfs_multiple_samples_biasdegree(graph_full_inner, N = N, k = k, S = S, num_cores = C)
network_descriptives_df <- compute_network_summary(results_bfs, num_cores = C)

#save network descriptives
write.csv(network_descriptives_df, file=paste0("outputBFSbiasdegree/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_descriptives_BFS_","run",run,".csv"))
#save results
save(results_bfs,network_descriptives_df, file=paste0("outputBFSbiasdegree/k",k,"/S",length(results_bfs),"N",564,"k",k,"network_sample_results_BFS_","run",run,".rdata"))


## Step 9: Random Walk sampling with j = 0.15

V(graph_full_inner)$id_str <- as.character(V(graph_full_inner)) #save original IDs before sampling creates new IDs

# Metropolis-Hastings Random Walk
metropolis_hastings_rw <- function(graph, N, j = 0.15, sample_number) {
  require(igraph)
  
  num_nodes <- vcount(graph)
  nodes <- as.integer(V(graph))  # Convert node set to integer indices
  degrees <- degree(graph, mode = "out")  # compute degrees
  
  visited_set <- logical(num_nodes)  #track unique nodes visited
  sampled_nodes <- integer(N)  # Store unique nodes
  full_walk <- integer(N * 2)  # Store full walk
  
  # Start at a random node
  current <- sample(nodes, 1)
  visited_set[current] <- TRUE
  sampled_nodes[1] <- current
  full_walk[1] <- current
  count <- 1  # Unique node counter
  step_count <- 1  # Track total steps taken
  
  # Move tracking
  random_jumps <- 0
  forced_jumps <- 0
  valid_moves <- 0
  
  while (count < N) {
    if (runif(1) < j) {
      # Random jump to any node
      current <- sample(nodes, 1)
      random_jumps <- random_jumps + 1
    } else {
      # Attempt to move to a neighbor
      neighbors <- as.integer(neighbors(graph, current, mode = "out"))
      if (length(neighbors) > 0) {
        proposed <- sample(neighbors, 1)
        alpha <- min(1, degrees[current] / degrees[proposed])
        if (runif(1) < alpha) {
          current <- proposed
          valid_moves <- valid_moves + 1
        }
      } else {
        # No outgoing edges, force a random jump
        current <- sample(nodes, 1)
        forced_jumps <- forced_jumps + 1
      }
    }
    
    # Store the move
    step_count <- step_count + 1
    full_walk[step_count] <- current
    
    # If new unique node, add to sample
    if (!visited_set[current]) {
      count <- count + 1
      sampled_nodes[count] <- current
      visited_set[current] <- TRUE
    }
  }
  
  subgraph <- induced_subgraph(graph, sampled_nodes)
  total_steps <- random_jumps + forced_jumps + valid_moves
  jump_percentage <- ifelse(total_steps > 0, (random_jumps + forced_jumps) / total_steps, 0)
  
  return(list(
    sample_number = sample_number,
    sampled_nodes = sampled_nodes,
    full_walk = full_walk[1:step_count],  
    subgraph = subgraph,
    jump_percentage = jump_percentage
  ))
}

# Parallelized Sampling Function
parallel_mhrw <- function(graph, N, j = 0.15, S, C) {
  cl <- makeCluster(C)
  
  clusterEvalQ(cl, library(igraph))  # Load igraph in parallel workers
  
  clusterExport(cl, varlist = c("graph", "N", "j", "metropolis_hastings_rw"), envir = environment())
  
  results <- parLapply(cl, 1:S, function(s) metropolis_hastings_rw(graph, N, j, s))
  
  stopCluster(cl)
  return(results)
}


#Initialization
run <- 1
N = 10373   #number of nodes per sample (final network size) -- BFS unbiased seeds averages: 10373, 295426
j = 0.15    # jump probability (from literature)
S = 10 # number of samples to draw
C = 10 #number of cores to use

results_rw <- parallel_mhrw(graph_full_inner, N = N, j = j, S = S, C = C)

#changed way to compute network summary due to runtime
network_descriptives_list <- mclapply(results_rw, function(res) {
  graph <- res$subgraph
  num_nodes <- igraph::vcount(graph)
  num_edges <- igraph::ecount(graph)
  comp <- igraph::components(graph, mode = "weak")
  jump_percentage <- res$jump_percentage
  
  data.frame(
    sample_number = res$sample_number,
    num_nodes = num_nodes,
    num_edges = num_edges,
    avg_indegree = mean(igraph::degree(graph, mode = "in")),
    avg_outdegree = mean(igraph::degree(graph, mode = "out")),
    density = igraph::edge_density(graph, loops = FALSE),
    clustering_coeff = igraph::transitivity(graph, type = "localaverage"),
    num_components = comp$no,
    largest_component_size = max(comp$csize),
    share_largest_comp = max(comp$csize) / 911331,
    jump_percentage = jump_percentage,
    assortativity_ethnic = assortativity_nominal(graph, V(graph)$group),
    assortativity_gender = assortativity_nominal(graph, V(graph)$group_gender),
    assortativity_party = assortativity_nominal(graph, V(graph)$group_party),
    clustering_coeff_global_vector= transitivity(graph, type = "global")
  )
}, mc.cores = C)
network_descriptives_df_rw <-
  rbind(network_descriptives_list[[1]],network_descriptives_list[[2]], network_descriptives_list[[2]],
        network_descriptives_list[[4]],network_descriptives_list[[5]],network_descriptives_list[[6]],
        network_descriptives_list[[7]],network_descriptives_list[[8]],network_descriptives_list[[9]],
        network_descriptives_list[[10]])


#save network descriptives
write.csv(network_descriptives_df_rw, file=paste0("outputRW/S",S,"N",N,"j",j,"network_sample_descriptives_RW_run",run,".csv"))
#save results
save(results_rw,network_descriptives_df_rw, file=paste0("outputRW/S",S,"N",N,"j",j,"network_sample_results_RW_run",run,".rdata"))


## Step 10: calculate attribute representations for all samples

attributes_full_df <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(attributes_full_df) <- c("sample_num","sample_method","sample_depth","sample_size","female","male","gender_unknown",
                                  "caucasian","african_american","asian","hispanic","race_other","race_unknown",
                                  "democrat","republican","independent")

#get statistics 
sample_size <- vcount(graph_full_inner)
female <- sum(V(graph_full_inner)$gender == "Female")
male <- sum(V(graph_full_inner)$gender == "Male")
gender_unknown <- sum(V(graph_full_inner)$gender == "Unknown")
caucasian <- sum(V(graph_full_inner)$race == "Caucasian")
african_american <- sum(V(graph_full_inner)$race == "African-American")
asian <- sum(V(graph_full_inner)$race == "Asian")
hispanic <- sum(V(graph_full_inner)$race == "Hispanic")
race_other <- sum(V(graph_full_inner)$race == "Other")
race_unknown <- sum(V(graph_full_inner)$race == "Unknown")
democrat <- sum(V(graph_full_inner)$party == "Democrat")
republican <- sum(V(graph_full_inner)$party == "Republican")
independent <- sum(V(graph_full_inner)$party == "Independent")

attributes_full_df <- rbind(attributes_full_df, c(0,"full","full",sample_size,female,male,gender_unknown,caucasian,african_american,asian,
                                                  hispanic,race_other,race_unknown,democrat,republican,independent))

colnames(attributes_full_df) <- c("sample_num","sample_method","sample_depth","sample_size","female","male","gender_unknown",
                                  "caucasian","african_american","asian","hispanic","race_other","race_unknown",
                                  "democrat","republican","independent")



# #### no bias ####

# The following function was used to create a dataframe with attribute representations for each sample samples by each sampling procedure

add_attributes <- function(result){
  
  attributes_full_df <- data.frame(matrix(ncol = 16, nrow = 0))
  colnames(attributes_full_df) <- c("sample_num","sample_method","sample_depth","sample_size","female","male","gender_unknown",
                                    "caucasian","african_american","asian","hispanic","race_other","race_unknown",
                                    "democrat","republican","independent")
  
  
  
  for (i in 1:length(results_bfs)){
    
    sample_num <- i
    sample_method <- "BFS" #this was "RW" for random walk samples
    sample_depth <-NA
    sample_size <- vcount(results_bfs[[i]]$subgraph)
    
    female <- sum(V(results_bfs[[i]]$subgraph)$gender == "Female")
    male <- sum(V(results_bfs[[i]]$subgraph)$gender == "Male")
    gender_unknown <- sum(V(results_bfs[[i]]$subgraph)$gender == "Unknown")
    caucasian <- sum(V(results_bfs[[i]]$subgraph)$race == "Caucasian")
    african_american <- sum(V(results_bfs[[i]]$subgraph)$race == "African-American")
    asian <- sum(V(results_bfs[[i]]$subgraph)$race == "Asian")
    hispanic <- sum(V(results_bfs[[i]]$subgraph)$race == "Hispanic")
    race_other <- sum(V(results_bfs[[i]]$subgraph)$race == "Other")
    race_unknown <- sum(V(results_bfs[[i]]$subgraph)$race == "Unknown")
    democrat <- sum(V(results_bfs[[i]]$subgraph)$party == "Democrat")
    republican <- sum(V(results_bfs[[i]]$subgraph)$party == "Republican")
    independent <- sum(V(results_bfs[[i]]$subgraph)$party == "Independent")
    
    attribute_df <- data.frame(sample_num,sample_method,sample_depth,sample_size,female,male,gender_unknown,
                               caucasian,african_american,asian,hispanic,race_other,race_unknown,
                               democrat, republican, independent)
    
    attributes_full_df <- rbind(attributes_full_df, attribute_df)
    
    
  }
  
  return(attributes_full_df)
}

load("FILENAME.rdata") #load a file containing all 100 samples (repeated for all BFS with any bias and for RW)

new_statistics <- add_attributes(results_bfs)
new_statistics$sample_depth <- "k1" # "k2" for larger samples

write.csv(new_statistics, file = "FILENAME_attributes.csv")


## Step 11: Reweighting attributes

edge_in_last_layer_detector <- function(k,result, user_id_graph_full_inner, sample_number){
  
  if (method == "donated"){
    subgraph <- result$subgraph
    
    seed_user_ids <- user_id_graph_full_inner[result$seed_nodes]
  }else{
    subgraph <- result[[sample_number]]$subgraph
    
    seed_user_ids <- user_id_graph_full_inner[result[[sample_number]]$seed_nodes]
  }
  #seeds
  V(subgraph)$is_seed <- V(subgraph)$userid %in% seed_user_ids
  
  #layer 1: neighbors of seeds
  neighbor_verts <- unlist(neighborhood(subgraph, order = 1, nodes = seed_user_ids, mode = "out"))
  V(subgraph)$is_layer1 <- FALSE
  V(subgraph)$is_layer1[neighbor_verts] <- TRUE
  V(subgraph)$is_layer1[V(subgraph)$is_seed == TRUE] <- FALSE
  
  #layer2: neighbors of neighbors of seeds
  V(subgraph)$is_layer2 <- FALSE
  V(subgraph)$is_layer2[V(subgraph)$is_seed == FALSE & V(subgraph)$is_layer1 == FALSE] <- TRUE
  
  if (k == 1) {
    
    last_layer <- V(subgraph)$is_layer1
    
  } else { if (k == 2) {
    last_layer <- V(subgraph)$is_layer2
    
  } else {cat("Error: k is neither 1 or 2")}}
  
  return(last_layer)
}


edge_layer_detector <- function(k, result, user_id_graph_full_inner, sample_number){
  
  if (method == "donated"){
    subgraph <- result$subgraph
  } else {
    subgraph <- result[[sample_number]]$subgraph
  }
  
  
  seed_user_ids <- user_id_graph_full_inner[result[[sample_number]]$seed_nodes]
  
  #seeds
  V(subgraph)$is_seed <- V(subgraph)$userid %in% seed_user_ids
  
  #layer 1: neigbors of seeds
  neighbor_verts <- unlist(neighborhood(subgraph, order = 1, nodes = seed_user_ids, mode = "out"))
  V(subgraph)$is_layer1 <- FALSE
  V(subgraph)$is_layer1[neighbor_verts] <- TRUE
  V(subgraph)$is_layer1[V(subgraph)$is_seed == TRUE] <- FALSE
  
  #layer2: neighbors of neigbors of seeds
  V(subgraph)$is_layer2 <- FALSE
  V(subgraph)$is_layer2[V(subgraph)$is_seed == FALSE & V(subgraph)$is_layer1 == FALSE] <- TRUE
  
  #layer attribute
  V(subgraph)$layer <- NA
  V(subgraph)$layer[V(subgraph)$is_seed == TRUE] <- 0
  V(subgraph)$layer[V(subgraph)$is_layer1 == TRUE] <- 1
  V(subgraph)$layer[V(subgraph)$is_layer2 == TRUE] <- 2
  
  return(subgraph)
}

#get real indegrees for each vertex in complete graph
ind_full <- degree(graph_full_inner, mode = "in")
names(ind_full) <- V(graph_full_inner)$userid
degree_full_df <- data.frame(useridchr = names(ind_full),
                             real_indegree = ind_full)

user_id_graph_full_inner <- V(graph_full_inner)$userid
rm(graph_full_inner)

#load attribute data from voter files
load("FILENAME.rdata") #load voter data

#combine degree from full graph with attribute df
voter_data_attributes_reduced <- voter_data_attributes_reduced %>%
  filter(useridchr %in% degree_full_df$useridchr) %>%
  left_join(., degree_full_df, by = "useridchr") %>%
  mutate(
    #gender
    female =  case_match(gender, "Female" ~ 1, .default = 0),
    male =  case_match(gender, "Male" ~ 1, .default = 0),
    nogender =  case_match(gender, "Unknown" ~ 1, .default = 0),
    #party
    democrat = case_match(party,"Democrat" ~1, .default = 0),
    republican = case_match(party,"Republican" ~1, .default = 0),
    independent = case_match(party,"Independent" ~1, .default = 0),
    #race
    caucasian = case_match(race,"Caucasian" ~1, .default = 0),
    asian = case_match(race,"Asian" ~1, .default = 0),
    hispanic = case_match(race,"Hispanic" ~1, .default = 0),
    afram = case_match(race,"African-American" ~1, .default = 0),
    other_race = case_match(race,"Other" ~1, .default = 0),
    unknown_race = case_match(race,"Unknown" ~1, .default = 0),
    #real indegree
    weight_real_indegree = if_else(real_indegree == 0,1,real_indegree))
head(voter_data_attributes_reduced)


# define hajek weight
hajek_share <- function(y,w){ #y is the estimate and w is the weight here indegree
  (sum(y/w, na.rm = TRUE)/sum(1/w, na.rm = TRUE))
}

#prepare df to save results
attributes_BFS_df_weight_attributes <- data.frame(matrix(ncol = 40, nrow = 0))
colnames(attributes_BFS_df_weight_attributes) <- c("sample_num","sample_method","sample_depth","sample_size",
                                                      "female_naive_mean","female_hajek_sampindegree","female_hajek_realindegree",
                                                      "male_naive_mean","male_hajek_sampindegree","male_hajek_realindegree",
                                                      "nogender_naive_mean","nogender_hajek_sampindegree","nogender_hajek_realindegree",
                                                      "democrat_naive_mean","democrat_hajek_sampindegree","democrat_hajek_realindegree",
                                                      "republican_naive_mean","republican_hajek_sampindegree","republican_hajek_realindegree",
                                                      "independent_naive_mean","independent_hajek_sampindegree","independent_hajek_realindegree",
                                                      "caucasian_naive_mean","caucasian_hajek_sampindegree","caucasian_hajek_realindegree",
                                                      "asian_naive_mean","asian_hajek_sampindegree","asian_hajek_realindegree",
                                                      "hispanic_naive_mean","hispanic_hajek_sampindegree","hispanic_hajek_realindegree",
                                                      "afram_naive_mean","afram_hajek_sampindegree","afram_hajek_realindegree",
                                                      "other_race_naive_mean","other_race_hajek_sampindegree","other_race_hajek_realindegree",
                                                      "unknown_race_naive_mean","unknown_race_hajek_sampindegree","unknown_race_hajek_realindegree"
)

for (method in c("BFS","donated","BFS_biasdegree_linear","RW0.15")){ 
  
  for (k in c(1,2)) {
    
    if (k == 1){
      ntotal = 10373
    } else{
      ntotal = 295426
    }
    
    #load sampled data
    if (method == "BFS" & k == 1) {
      load("FILENAME.rdata") # unbiased BFS sample k = 1
      runs <- 1
    } else {
      if (method == "BFS_biasdegree_linear"  & k == 1){
        load("FILENAME.rdata") #biasdegree sample k = 1
        runs <- 1
      } else {
        if (method == "donated") {
          load(paste0("FILENAME.rdata"))
          runs <- 1
        } else {
          runs <- c(1:10)
        }
      }
    }
    
    for (run in runs) {
      
      if (method == "BFS" & k == 2){
        load(paste0("FILENAME.rdata"))
      } else {
        if (method == "BFS_biasdegree_linear"  & k == 2){
          load(paste0("FILENAME.rdata"))
        } else {
            if (method == "RW0.15"){
              load(paste0("FILENAME.rdata"))
              # this is rw
              results_bfs <- results_rw #assign name for loop below
            }
          }
        
      }
      
      
      if (method == "donated"){
        num_results <- 1
      } else{
        num_results <- length(results_bfs)
      }
      
      
      for (sample in 1:num_results){
        
        if (method == "donated"){
          subgraph <- results_bfs$subgraph
          subset_voterfile <- voter_data_attributes_reduced %>%
            filter(useridchr %in% V(subgraph)$useridchr) %>%
            arrange(match(useridchr, V(subgraph)$useridchr)) %>%
            mutate(indegree_sample = degree(subgraph, mode = "in"),
                   weight_indegree_sample = if_else(indegree_sample == 0,1,indegree_sample),
                   is_seed = V(subgraph)$userid %in% user_id_graph_full_inner[results_bfs$seed_nodes],
                   last_layer = edge_in_last_layer_detector(k,results_bfs, user_id_graph_full_inner,sample))
        } else{
          subgraph <- results_bfs[[sample]]$subgraph
          subset_voterfile <- voter_data_attributes_reduced %>%
            filter(useridchr %in% V(subgraph)$useridchr) %>%
            arrange(match(useridchr, V(subgraph)$useridchr)) %>%
            mutate(indegree_sample = degree(subgraph, mode = "in"),
                   weight_indegree_sample = if_else(indegree_sample == 0,1,indegree_sample),
                   is_seed = V(subgraph)$userid %in% user_id_graph_full_inner[results_bfs[[sample]]$seed_nodes],
                   last_layer = edge_in_last_layer_detector(k,results_bfs, user_id_graph_full_inner,sample))
        }
        #get statistics 
        sample_size <- vcount(subgraph)
        
        #female ratio
        female_naive_mean <- length(subset_voterfile$gender[subset_voterfile$gender == "Female"]) / sample_size
        female_hajek_sampindegree <- hajek_share(subset_voterfile$female,subset_voterfile$weight_indegree_sample)
        female_hajek_realindegree <- hajek_share(subset_voterfile$female,subset_voterfile$weight_real_indegree)
        
        #male ratio
        male_naive_mean <- length(subset_voterfile$gender[subset_voterfile$gender == "Male"]) / sample_size
        male_hajek_sampindegree <- hajek_share(subset_voterfile$male,subset_voterfile$weight_indegree_sample)
        male_hajek_realindegree <- hajek_share(subset_voterfile$male,subset_voterfile$weight_real_indegree)
        
        #gender unknown ratio
        nogender_naive_mean <- length(subset_voterfile$gender[subset_voterfile$gender == "Unknown"]) / sample_size
        nogender_hajek_sampindegree <- hajek_share(subset_voterfile$nogender,subset_voterfile$weight_indegree_sample)
        nogender_hajek_realindegree <- hajek_share(subset_voterfile$nogender,subset_voterfile$weight_real_indegree)
        
        #democrat ratio
        democrat_naive_mean <- length(subset_voterfile$party[subset_voterfile$party == "Democrat"]) / sample_size
        democrat_hajek_sampindegree <- hajek_share(subset_voterfile$democrat,subset_voterfile$weight_indegree_sample)
        democrat_hajek_realindegree <- hajek_share(subset_voterfile$democrat,subset_voterfile$weight_real_indegree)
        
        #republican ratio
        republican_naive_mean <- length(subset_voterfile$party[subset_voterfile$party == "Republican"]) / sample_size
        republican_hajek_sampindegree <- hajek_share(subset_voterfile$republican,subset_voterfile$weight_indegree_sample)
        republican_hajek_realindegree <- hajek_share(subset_voterfile$republican,subset_voterfile$weight_real_indegree)
        
        #independent ratio
        independent_naive_mean <- length(subset_voterfile$party[subset_voterfile$party == "Independent"]) / sample_size
        independent_hajek_sampindegree <- hajek_share(subset_voterfile$independent,subset_voterfile$weight_indegree_sample)
        independent_hajek_realindegree <- hajek_share(subset_voterfile$independent,subset_voterfile$weight_real_indegree)
        
        #caucasian ratio
        caucasian_naive_mean <- length(subset_voterfile$race[subset_voterfile$race == "Caucasian"]) / sample_size
        caucasian_hajek_sampindegree <- hajek_share(subset_voterfile$caucasian,subset_voterfile$weight_indegree_sample)
        caucasian_hajek_realindegree <- hajek_share(subset_voterfile$caucasian,subset_voterfile$weight_real_indegree)
        
        #asian ratio
        asian_naive_mean <- length(subset_voterfile$race[subset_voterfile$race == "Asian"]) / sample_size
        asian_hajek_sampindegree <- hajek_share(subset_voterfile$asian,subset_voterfile$weight_indegree_sample)
        asian_hajek_realindegree <- hajek_share(subset_voterfile$asian,subset_voterfile$weight_real_indegree)
        
        #hispanic ratio
        hispanic_naive_mean <- length(subset_voterfile$race[subset_voterfile$race == "Hispanic"]) / sample_size
        hispanic_hajek_sampindegree <- hajek_share(subset_voterfile$hispanic,subset_voterfile$weight_indegree_sample)
        hispanic_hajek_sampindegree <- hajek_share(subset_voterfile$hispanic,subset_voterfile$weight_real_indegree)
        
        #African-American ratio
        afram_naive_mean <- length(subset_voterfile$race[subset_voterfile$race == "African-American"]) / sample_size
        afram_hajek_sampindegree <- hajek_share(subset_voterfile$afram,subset_voterfile$weight_indegree_sample)
        afram_hajek_realindegree <- hajek_share(subset_voterfile$afram,subset_voterfile$weight_real_indegree)
        
        #other etnicity ratio
        other_race_naive_mean <- length(subset_voterfile$race[subset_voterfile$race == "Other"]) / sample_size
        other_race_hajek_sampindegree <- hajek_share(subset_voterfile$other_race,subset_voterfile$weight_indegree_sample)
        other_race_hajek_realindegree <- hajek_share(subset_voterfile$other_race,subset_voterfile$weight_real_indegree)
        
        #unknown etnicity ratio
        unknown_race_naive_mean <- length(subset_voterfile$race[subset_voterfile$race == "Unknown"]) / sample_size
        unknown_race_hajek_sampindegree <- hajek_share(subset_voterfile$unknown_race,subset_voterfile$weight_indegree_sample)
        unknown_race_hajek_realindegree <- hajek_share(subset_voterfile$unknown_race,subset_voterfile$weight_real_indegree)
        
        if(length(runs) == 1) {
          sample_num <- sample
        } else {
          sample_num <- (run-1)*10+sample}
        
        df_sample <- c(sample_num,method,paste0("k",k),sample_size,
                       female_naive_mean,female_hajek_sampindegree,female_hajek_realindegree,
                       male_naive_mean,male_hajek_sampindegree,male_hajek_realindegree,
                       nogender_naive_mean,nogender_hajek_sampindegree,nogender_hajek_realindegree,
                       democrat_naive_mean,democrat_hajek_sampindegree,democrat_hajek_realindegree, 
                       republican_naive_mean,republican_hajek_sampindegree,republican_hajek_realindegree,
                       independent_naive_mean,independent_hajek_sampindegree,independent_hajek_realindegree,
                       caucasian_naive_mean,caucasian_hajek_sampindegree,caucasian_hajek_realindegree,
                       asian_naive_mean,asian_hajek_sampindegree,asian_hajek_realindegree,
                       hispanic_naive_mean,hispanic_hajek_sampindegree,hispanic_hajek_sampindegree,
                       afram_naive_mean,afram_hajek_sampindegree,afram_hajek_realindegree,
                       other_race_naive_mean,other_race_hajek_sampindegree,other_race_hajek_realindegree,
                       unknown_race_naive_mean,unknown_race_hajek_sampindegree,unknown_race_hajek_realindegree
        )
        
        
        
        attributes_BFS_df_weight_attributes <- rbind(attributes_BFS_df_weight_attributes, setNames(as.data.frame(t(df_sample)), names(attributes_BFS_k1_df_weight_attributes)))
        
      }
    }
    
  }}

write.csv(attributes_BFS_df_weight_attributes, file = "attributes_weights.csv")


## Step 12: Reweighting statistics: density, local clc and ethnic assortativity

# global clc treated separately because of computational reasons

add_edge_sampling_weights <- function(subgraph,
                                      last_layer, # 1 or 2
                                      zero_prob_weight = c("NA", "Inf")) {
  
  zero_prob_weight <- match.arg(zero_prob_weight)
  
  # Layers
  layer <- vertex_attr(subgraph, "layer") #for RW == k because seed and layer1 were false
  
  # Indegrees (if 0 then 1 as for nodal attributes)
  weight_indegree_sample <- subset_voterfile$weight_indegree_sample
  weight_real_indegree <- subset_voterfile$weight_real_indegree
  
  #edgelist
  el <- as_edgelist(subgraph, names = FALSE)
  p_e_weight_indegree_sample <- numeric(nrow(el))
  p_e_weight_real_indegree <- numeric(nrow(el))
  
  for (i in seq_len(nrow(el))) {
    u <- el[i, 1]
    v <- el[i, 2]
    
    lu <- layer[u]
    lv <- layer[v]
    
    if (lu == last_layer && lv == last_layer && !(method == "RW0.15")) {
      # both nodes in the last sampled layer
      p_e_weight_indegree_sample[i] <- 0
      p_e_weight_real_indegree[i] <- 0
    } else if (lu < lv && !(method == "RW0.15")) {
      # u sampled earlier
      p_e_weight_indegree_sample[i] <- weight_indegree_sample[u]
      p_e_weight_real_indegree[i] <- weight_real_indegree[u]
    } else if (lv < lu && !(method == "RW0.15")) {
      # v sampled earlier
      p_e_weight_indegree_sample[i] <- weight_indegree_sample[v]
      p_e_weight_real_indegree[i] <- weight_real_indegree[v]
    } else {
      # same layer, but not both in last_layer OR RW because layers do not exist
      p_e_weight_indegree_sample[i] <- weight_indegree_sample[u] + weight_indegree_sample[v]
      p_e_weight_real_indegree[i] <- weight_real_indegree[u] + weight_real_indegree[v]
    }
  }
  
  # Convert to inverse weights
  w_e_weight_indegree_sample <- ifelse(
    p_e_weight_indegree_sample == 0,
    if (zero_prob_weight == "Inf") Inf else NA_real_,
    1 / p_e_weight_indegree_sample
  )
  
  w_e_weight_real_indegree <- ifelse(
    p_e_weight_real_indegree == 0,
    if (zero_prob_weight == "Inf") Inf else NA_real_,
    1 / p_e_weight_real_indegree
  )
  
  #adds edge probabilities as edge attributes
  subgraph <- set_edge_attr(subgraph, "p_e_weight_indegree_sample", value = p_e_weight_indegree_sample)
  subgraph <- set_edge_attr(subgraph, "p_e_weight_real_indegree", value = p_e_weight_real_indegree)
  
  #adds inverse probabilities as weights as edge attributes
  subgraph <- set_edge_attr(subgraph, "w_e_weight_indegree_sample", value = w_e_weight_indegree_sample)
  subgraph <- set_edge_attr(subgraph, "w_e_weight_real_indegree", value = w_e_weight_real_indegree)
  
  return(subgraph)
}

#weighted density
weighted_edge_density <- function(g, weight_attr = "edge_w") {
  stopifnot(is_igraph(g))
  w <- edge_attr(g, weight_attr)
  if (is.null(w)) stop(sprintf("Edge attribute '%s' not found.", weight_attr))
  
  n <- vcount(g)
  denom <- n * (n - 1)
  
  sum(w) / denom
}

# weighted ethnic assortativity
weighted_assortativity_nominal <- function(g,
                                           group_attr = "group",
                                           weight_attr = "edge_w") {
  stopifnot(is_igraph(g))
  
  grp <- vertex_attr(g, group_attr)
  if (is.null(grp)) stop(sprintf("Vertex attribute '%s' not found.", group_attr))
  
  w <- edge_attr(g, weight_attr)
  if (is.null(w)) stop(sprintf("Edge attribute '%s' not found.", weight_attr))
  
  grp <- as.factor(grp)
  lev <- levels(grp)
  K <- length(lev)
  
  el <- as_edgelist(g, names = FALSE)
  
  M <- matrix(0, nrow = K, ncol = K, dimnames = list(lev, lev))
  
  for (i in seq_len(nrow(el))) {
    u <- el[i, 1]
    v <- el[i, 2]
    gu <- as.character(grp[u])
    gv <- as.character(grp[v])
    
    M[gu, gv] <- M[gu, gv] + w[i]
  }
  
  if (sum(M, na.rm = TRUE) == 0) return(NA_real_)
  
  e <- M / sum(M, na.rm = TRUE)
  a <- rowSums(e)
  b <- colSums(e)
  
  num <- sum(diag(e), na.rm = TRUE) - sum(a * b, na.rm = TRUE)
  den <- 1 - sum(a * b, na.rm = TRUE)
  
  if (den == 0) return(NA_real_)
  num / den
}

#weighted global clustering coef

weighted_global_clustering <- function(g, weight_attr = "edge_w") { # "noweight" if all weights are 1 for testing purposes
  stopifnot(is_igraph(g))
  
  w <- edge_attr(g, weight_attr)
  if (is.null(w)) stop(sprintf("Edge attribute '%s' not found.", weight_attr))
  
  # Collapse directed graph to undirected and sum reciprocal edge weights
  gu <- as_undirected(
    g,
    mode = "collapse",
    edge.attr.comb = list(edge_w = "sum", "ignore")
  )
  
  if (weight_attr == "noweight") {
    E(gu)$noweight <- 1
  }
  
  ew <- edge_attr(gu, weight_attr)
  ew[is.na(ew)] <- 0  # for stability
  
  # Precomp adjacency lists and edge id lookup
  adj <- adjacent_vertices(gu, V(gu), mode = "all")
  
  # Fast edge lookup
  edge_id <- function(u, v) {
    get_edge_ids(gu, c(u, v), directed = FALSE)
  }
  
  n <- vcount(gu)
  
  # Denominator: weighted triples
  s <- strength(gu, mode = "all", weights = ew)
  
  inc_edges <- incident_edges(gu, V(gu))
  sum_sq <- vapply(
    inc_edges,
    function(eids) {
      if (length(eids) == 0) return(0)
      sum(ew[eids]^2)
    },
    numeric(1)
  )
  
  triples_w <- sum((s^2 - sum_sq) / 2)
  if (triples_w == 0) return(NA_real_)
  
 # Numerator: streaming triangle enumeration
  tri_w <- 0
  
  for (i in seq_len(n)) {
    Ni <- adj[[i]]
    if (length(Ni) < 2) next
    
    Ni <- as.integer(Ni)
    
    # Only consider neighbors with higher index to avoid duplicates
    Ni <- Ni[Ni > i]
    if (length(Ni) < 2) next
    
    for (j_idx in seq_along(Ni)) {
      j <- Ni[j_idx]
      
      Nj <- adj[[j]]
      Nj <- as.integer(Nj)
      
      # Intersection of neighbors of i and j
      common <- intersect(Ni, Nj)
      
      # Only keep k > j to enforce ordering i < j < k
      common <- common[common > j]
      
      if (length(common) == 0) next
      
      e_ij <- edge_id(i, j)
      w_ij <- ew[e_ij]
      
      for (k in common) {
        e_ik <- edge_id(i, k)
        e_jk <- edge_id(j, k)
        
        tri_w <- tri_w + w_ij * ew[e_ik] * ew[e_jk]
      }
    }
  }
  
  3 * tri_w / triples_w
}


#get real indegrees for each vertex in complet graph
ind_full <- degree(graph_full_inner, mode = "in")
names(ind_full) <- V(graph_full_inner)$userid
degree_full_df <- data.frame(useridchr = names(ind_full),
                             real_indegree = ind_full)

user_id_graph_full_inner <- V(graph_full_inner)$userid
rm(graph_full_inner)



#combine degree from full graph with attribute df
voter_data_attributes_reduced <- voter_data_attributes_reduced %>%
  filter(useridchr %in% degree_full_df$useridchr) %>%
  left_join(., degree_full_df, by = "useridchr") %>%
  mutate(
    #gender
    female =  case_match(gender, "Female" ~ 1, .default = 0),
    male =  case_match(gender, "Male" ~ 1, .default = 0),
    nogender =  case_match(gender, "Unknown" ~ 1, .default = 0),
    #party
    democrat = case_match(party,"Democrat" ~1, .default = 0),
    republican = case_match(party,"Republican" ~1, .default = 0),
    independent = case_match(party,"Independent" ~1, .default = 0),
    #race
    caucasian = case_match(race,"Caucasian" ~1, .default = 0),
    asian = case_match(race,"Asian" ~1, .default = 0),
    hispanic = case_match(race,"Hispanic" ~1, .default = 0),
    afram = case_match(race,"African-American" ~1, .default = 0),
    #real indegree
    weight_real_indegree = if_else(real_indegree == 0,1,real_indegree))
head(voter_data_attributes_reduced)


# define hajek weight
hajek_share <- function(y,w){ #y is the estimate and w is the weight here indegree
  (sum(y/w, na.rm = TRUE)/sum(1/w, na.rm = TRUE))
}


#prepare df to save results
statistics_df_weight <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(statistics_df_weight) <- c("sample_num","sample_method","sample_depth","sample_size",
                                    "lcl_coef_naive_mean","lcl_coef_hajek_sampindegree","lcl_coef_hajek_realindegree",
                                    "density_naive","density_hajek_sampindegree","density_hajek_realindegree",
                                    "assortativity_naive","assortativity_hajek_sampindegree","assortativity_hajek_realindegree"
)

for (method in c("donated","BFS","BFS_biasdegree_linear","RW0.15")){ 
  
  for (k in c(1,2)) {
    
    if (k == 1){
      ntotal = 10373 
    } else{
      ntotal = 295426 
    }
    
    #load sampled data
    if (method == "BFS" & k == 1) {
      load("FILENAME.rdata") # unbiased BFS sample k = 1
      runs <- 1
    } else {
      if (method == "BFS_biasdegree_linear"  & k == 1){
        load("FILENAME.rdata") #biasdegree sample k = 1
        runs <- 1
      } else {
        if (method == "donated") {
          load(paste0("FILENAME",k,"restofFILENAME.rdata"))
          runs <- 1
        } else {
          runs <- c(1:10)
        }
      }
    }
    
    for (run in runs) {
      
      if (method == "BFS" & k == 2){
        load(paste0("FILENAME",run,".rdata"))
      } else {
        if (method == "BFS_biasdegree_linear"  & k == 2){
          load(paste0("FILENAME",run,".rdata"))
        } else {
            if (method == "RW0.15"){
              load(paste0("FILENAME",k,"...",ntotal,"...",run,".rdata"))
              # this is rw
              results_bfs <- results_rw #assign name for loop below
            }
          
        }
        
      }
      
      
      if (method == "donated"){
        num_results <- 1
      } else{
        num_results <- length(results_bfs)
      }
      
      for (sample in 1:num_results){
        
        subgraph <- edge_layer_detector(k,results_bfs, user_id_graph_full_inner,sample) #adds layer as atribute and returns subgraph
        
        #results_bfs[[sample]]$subgraph
        subset_voterfile <- voter_data_attributes_reduced %>%
          filter(useridchr %in% V(subgraph)$useridchr) %>%
          arrange(match(useridchr, V(subgraph)$useridchr)) %>%
          mutate(indegree_sample = degree(subgraph, mode = "in"),
                 weight_indegree_sample = if_else(indegree_sample == 0,1,indegree_sample),
                 is_seed = V(subgraph)$userid %in% user_id_graph_full_inner[results_bfs[[sample]]$seed_nodes],
                 lcl_coef = transitivity(subgraph, type = "local"),
                 outdegree_sample = degree(subgraph, mode = "out"),
                 last_layer = V(subgraph)$layer == k)
        
        # add different indegree estimates to edges
        subgraph <- add_edge_sampling_weights(subgraph,k,zero_prob_weight = "NA") #adds edge_incl_prob and edge_w, the inverse inclusion prob as edge attributes
        
        # add group attribute 
        V(subgraph)$group <- NA
        V(subgraph)$group[V(subgraph)$race == "Caucasian"] <- 1
        V(subgraph)$group[V(subgraph)$race == "Asian"]     <- 2
        V(subgraph)$group[V(subgraph)$race == "African-American"] <- 3
        V(subgraph)$group[V(subgraph)$race == "Hispanic"]     <- 4
        V(subgraph)$group[V(subgraph)$race == "Other"] <- 5
        V(subgraph)$group[V(subgraph)$race == "Unknown"]     <- 6
        
        #get statistics 
        sample_size <- vcount(subgraph)
        
        #local clustering coef
        lcl_coef_naive_mean <- mean(subset_voterfile$lcl_coef, na.rm = T)# same as transitivity local average
        lcl_coef_hajek_sampindegree <- hajek_share(subset_voterfile$lcl_coef[!is.na(subset_voterfile$lcl_coef)],subset_voterfile$weight_indegree_sample[!is.na(subset_voterfile$lcl_coef)])
        lcl_coef_hajek_realindegree <- hajek_share(subset_voterfile$lcl_coef[!is.na(subset_voterfile$lcl_coef)],subset_voterfile$weight_real_indegree[!is.na(subset_voterfile$lcl_coef)])
        
        #density
        density_naive <- edge_density(subgraph, loops = FALSE)
        density_hajek_sampindegree <- weighted_edge_density(subgraph, weight_attr = "w_e_weight_indegree_sample")
        density_hajek_realindegree <- weighted_edge_density(subgraph, weight_attr = "w_e_weight_real_indegree")
        
        
        #ethnic assortativity
        assortativity_naive <- assortativity_nominal(subgraph, V(subgraph)$group)
        assortativity_hajek_sampindegree <- weighted_assortativity_nominal(subgraph, group_attr = "group", weight_attr = "w_e_weight_indegree_sample")
        assortativity_hajek_realindegree <- weighted_assortativity_nominal(subgraph, group_attr = "group", weight_attr = "w_e_weight_real_indegree")
        
        
        
        # NOTE for global clc, density and assortativity: this can only be corrected by edge inclusion probability and not by node inclusion probability,
        #     in the logic of BFS edge inclusion probability can be determined by the nodal inclusion probability of the node that was sampled in the lower layer,
        #     edges between nodes sampled in the last layer have an inclusion probability of zero
        
        
        
        if(length(runs) == 1) {
          sample_num <- sample
        } else {
          sample_num <- (run-1)*10+sample}
        
        df_sample <- c(sample_num,method,paste0("k",k),sample_size,
                       lcl_coef_naive_mean,
                       lcl_coef_hajek_sampindegree,
                       lcl_coef_hajek_realindegree,
                       density_naive,
                       density_hajek_sampindegree,
                       density_hajek_realindegree,
                       assortativity_naive,
                       assortativity_hajek_sampindegree,
                       assortativity_hajek_realindegree
        )
        
        statistics_df_weight <- rbind(statistics_df_weight, setNames(as.data.frame(t(df_sample)), names(statistics_df_weight)))
        
      }
    }
    cat(paste0("run complete",method))
  }}

write.csv(statistics_df_weight, file = "statistics_df_weights_noedgeslastlayer_k1.csv")


## Step 13: Reweighting statistics: global clc (only for k1, k2 failed after 3.5 weeks)


#prepare df to save results
statistics_df_weight <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(statistics_df_weight) <- c("sample_num","sample_method","sample_depth","sample_size",
                                    "global_clc_hajek_naive","global_clc_hajek_sampindegree","global_clc_hajek_realindegree"
)

for (method in c("donated","BFS","BFS_biasdegree_linear","RW0.15")){ 
  
  for (k in c(1)) { #,2
    
    if (k == 1){
      ntotal = 10373 
    } else{
      ntotal = 295426 
    }
    
    #load sampled data
    if (method == "BFS" & k == 1) {
      load("FILENAME.rdata") # unbiased BFS sample k = 1
      runs <- 1
    } else {
      if (method == "BFS_biasdegree_linear"  & k == 1){
        load("FILENAME.rdata") #biasdegree sample k = 1
        runs <- 1
      } else {
        if (method == "donated") {
          load(paste0("FILENAME",k,".rdata"))
          runs <- 1
        } else {
          runs <- c(1:10)
        }
      }
    }
    
    for (run in runs) {
      
      if (method == "BFS" & k == 2){
        load(paste0("FILENAME",run,".rdata"))
      } else {
        if (method == "BFS_biasdegree_linear"  & k == 2){
          load(paste0("FILENAME",run,".rdata"))
        } else{
            if (method == "RW0.15"){
              load(paste0("FILENAME",k,"...",ntotal,"...",run,".rdata"))
              # this is rw
              results_bfs <- results_rw #assign name for loop below
            }
          }
        
      }
      
      
      if (method == "donated"){
        num_results <- 1 
      } else{
        num_results <- length(results_bfs)
      }
      
      for (sample in 1:num_results){
        
        subgraph <- edge_layer_detector(k,results_bfs, user_id_graph_full_inner,sample) #adds layer as attribute and returns subgraph
        if (method == "donated"){
          seed_list <- user_id_graph_full_inner[results_bfs$seed_nodes]
        } else{
          seed_list <- user_id_graph_full_inner[results_bfs[[sample]]$seed_nodes]
        }
        
        #results_bfs[[sample]]$subgraph
        subset_voterfile <- voter_data_attributes_reduced %>%
          filter(useridchr %in% V(subgraph)$useridchr) %>%
          arrange(match(useridchr, V(subgraph)$useridchr)) %>%
          mutate(indegree_sample = degree(subgraph, mode = "in"),
                 weight_indegree_sample = if_else(indegree_sample == 0,1,indegree_sample),
                 is_seed = V(subgraph)$userid %in% seed_list,
                 lcl_coef = transitivity(subgraph, type = "local"),
                 outdegree_sample = degree(subgraph, mode = "out"),
                 last_layer = V(subgraph)$layer == k)
        
        # add different indegree estimates to edges
        subgraph <- add_edge_sampling_weights(subgraph,k,zero_prob_weight = "NA") #adds edge_incl_prob and edge_w, the inverse inclusion prob as edge attributes
        
        #get statistics 
        sample_size <- vcount(subgraph)
        
        # #global clustering coef
        global_clc_hajek_naive <- transitivity(subgraph, type = "global")
        global_clc_hajek_sampindegree <-weighted_global_clustering(subgraph, k = k, method = method, weight_attr = "w_e_weight_indegree_sample", layer_attr = "layer")
        global_clc_hajek_realindegree <- weighted_global_clustering(subgraph, k = k, method = method, weight_attr = "w_e_weight_real_indegree", layer_attr = "layer")
        
        if(length(runs) == 1) {
          sample_num <- sample
        } else {
          sample_num <- (run-1)*10+sample}
        
        df_sample <- c(sample_num,method,paste0("k",k),sample_size,
                       global_clc_hajek_naive,
                       global_clc_hajek_sampindegree,
                       global_clc_hajek_realindegree
        )
        
        statistics_df_weight <- rbind(statistics_df_weight, setNames(as.data.frame(t(df_sample)), names(statistics_df_weight)))
        
      }
    }
    cat(paste0("run complete",method))
  }}

write.csv(statistics_df_weight, file = "statistics_df_weights_noedgeslastlayer_global_clc_k1.csv")


## Step 14: NMSE unweighted estimates

### attributes

#full network statistics
#gender
share_male_full <- full_graph_attributes$male /full_graph_attributes$sample_size
share_female_full <- full_graph_attributes$female /full_graph_attributes$sample_size
share_gender_unknown_full <- full_graph_attributes$gender_unknown /full_graph_attributes$sample_size

#ethnicity
share_caucasian_full <- full_graph_attributes$caucasian /full_graph_attributes$sample_size
share_asian_full <- full_graph_attributes$asian /full_graph_attributes$sample_size
share_afram_full <- full_graph_attributes$african_american /full_graph_attributes$sample_size
share_hispanic_full <- full_graph_attributes$hispanic /full_graph_attributes$sample_size
share_other_full <- full_graph_attributes$race_other /full_graph_attributes$sample_size
share_unknown_full <- full_graph_attributes$race_unknown /full_graph_attributes$sample_size

#party affiliation
share_democrat_full <- full_graph_attributes$democrat /full_graph_attributes$sample_size
share_republican_full <- full_graph_attributes$republican /full_graph_attributes$sample_size
share_independent_full <- full_graph_attributes$independent /full_graph_attributes$sample_size

#only works for no edges last layer
RMSE_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(RMSE_df) <- c("sample_method","sample_bias","sample_depth","gender_NMSE","ethnicity_NMSE","party_NMSE")

bias_types <- c("no bias","biasdegree_linear","partybias_democrat5percentmore","partybias_democrat10percentmore","partybias_democrat100percentmore",
                "partybias_republican5percentmore","partybias_republican10percentmore","partybias_republican100percentmore",
                "genderbias_women5percentmore","genderbias_women10percentmore","genderbias_women100percentmore",
                "genderbias_men5percentmore","genderbias_men10percentmore","genderbias_men100percentmore",
                "ethnicity_caucasian5percentmore","ethnicity_caucasian10percentmore","ethnicity_caucasian100percentmore",
                "ethnicity_asian5percentmore","ethnicity_asian10percentmore","ethnicity_asian100percentmore","RW0.15","donated")

for (bias in 1:length(bias_types)){
  
  for(k in 1:2){
    
    sample_depth <- k
    sample_method_bias <- bias_types[bias]
    
    
    if (!sample_method_bias %in% c("RW0.15", "donated", "no bias")){
      attributes_samples <- read.csv(paste0(input,"S100N564k",k,"network_sample_attributes_BFS_",bias_types[bias],".csv"), sep = ",")
    } else {
      if (k == 1){
        NRW <- 10373 
      } else {
        NRW <- 295426 
      }
      
      if (sample_method_bias == "RW0.15") {
        attributes_samples <- read.csv("FILENAME")
      } else {
          if (sample_method_bias == "no bias"){
            attributes_samples <- read.csv("FILENAME")
          } else{
            # "donated"
            attributes_samples <- read.csv("FILENAME") %>%
              mutate(sample_depth = paste0("k",k))
          }}
      }
    
    if (sample_method_bias == "") {
      sample_method_bias <-"no bias"
    }
    
    
    attributes_samples <- attributes_samples %>%
      mutate(share_female = female / sample_size,
             share_male = male / sample_size,
             share_nogender = gender_unknown / sample_size,
             share_caucasian = caucasian / sample_size,
             share_asian = asian / sample_size,
             share_afram = african_american / sample_size,
             share_hispanic = hispanic / sample_size,
             share_race_other = race_other / sample_size,
             share_race_unknown = race_unknown / sample_size,
             share_democrat = democrat / sample_size,
             share_republican = republican / sample_size,
             share_independent = independent / sample_size,
             gender_NMSE = ((share_male_full - share_male)^2 + (share_female_full - share_female)^2 + (share_gender_unknown_full - share_nogender)^2) / 
               (share_male_full^2 + share_female_full^2 + share_gender_unknown_full^2),
             ethnicity_NMSE = ((share_caucasian_full - share_caucasian)^2 + (share_asian_full - share_asian)^2 + (share_afram_full - share_afram)^2+
                                 (share_hispanic_full - share_hispanic)^2 + (share_other_full - share_race_other)^2 + (share_unknown_full - share_race_unknown)^2)/ 
               (share_caucasian_full^2 + share_asian_full^2 + share_afram_full^2 + share_hispanic_full^2 + share_other_full^2 + share_unknown_full^2),
             party_NMSE = ((share_democrat_full - share_democrat)^2 + (share_republican_full - share_republican)^2 + (share_independent_full - share_independent)^2) / 
               (share_democrat_full^2 + share_republican_full^2 + share_independent_full^2),
             sample_depth = sample_depth,
             sample_method = sample_method,
             sample_bias = sample_method_bias
      ) %>%
      select(sample_method,sample_bias,sample_depth,gender_NMSE,ethnicity_NMSE,party_NMSE)
    
    RMSE_df <- rbind(RMSE_df,attributes_samples)
    
    cat("run complete", sample_method_bias)
    
  }}

#when table retrieved then group by sampling method and calculate average

NMSE_df_summary <- RMSE_df %>% 
  group_by(sample_method,sample_bias,sample_depth) %>%
  summarise(gender_NMSE_av = mean(gender_NMSE),
            gender_NMSE_se = sd(gender_NMSE) / sqrt(n()),
            ethnicity_NMSE_av = mean(ethnicity_NMSE),
            ethnicity_NMSE_se = sd(ethnicity_NMSE) / sqrt(n()),
            party_NMSE_av = mean(party_NMSE),
            party_NMSE_se = sd(party_NMSE) / sqrt(n()),
            N_samples = n())


write.csv(NMSE_df_summary, "NMSE_attributes_biases_new.csv")



#### global network statistics ####

full_graph_statistics <- read.csv("whole_network_descriptives_in_voterfile.csv", sep = ",")

#full graph statistics
avg_degree_full <- full_graph_statistics$avg_indegree
density_full <- full_graph_statistics$density
avg_local_clustering_coeff_full <- full_graph_statistics$clustering_coeff
clustering_coeff_global_full <- full_graph_statistics$clustering_coeff_global_vector
assortativity_ethnic_full <- full_graph_statistics$assortativity_ethnic


# this is only for no edges last layer
RMSE_df <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(RMSE_df) <- c("sample_method_bias","sample_depth","avg_degree_NMSE","density_NMSE","avg_local_clustering_coeff_NMSE","clustering_coeff_global_NMSE","assortativity_ethnic_NMSE")


for (bias in 1:length(bias_types)){
  
  for(k in 1:2){
    
    sample_depth <- paste0("k",k)
    sample_method_bias <- bias_types[bias]
    
    
    if (!sample_method_bias %in% c("RW0.15", "donated", "no bias")){#
      
      samples_statistics <- read.csv("FILENAME")
      
    } else {
      if (k == 1){
        NRW <- 10373 
      } else {
        NRW <- 295426 
      }
      
      if (sample_method_bias == "RW0.15") {
        samples_statistics <- read.csv("FILENAME")
      } else { if (sample_method_bias == "no bias"){
          samples_statistics <- read.csv("FILENAME")
          
        } else{# "donated"
          samples_statistics <- read.csv("FILENAME") %>%
            mutate(sample_depth = paste0("k",k),
                   avg_indegree = avg_degree) 
        }}
      }
    
    
    samples_statistics <- samples_statistics %>%
      mutate(avg_degree_NMSE = (avg_degree_full - avg_indegree)^2 / 
               avg_degree_full^2,
             density_NMSE = (density_full - density)^2 / 
               density_full^2,
             avg_local_clustering_coeff_NMSE = (avg_local_clustering_coeff_full - clustering_coeff)^2 / 
               avg_local_clustering_coeff_full^2,
             clustering_coeff_global_full_NMSE = (clustering_coeff_global_full - clustering_coeff_global_vector)^2 / 
               clustering_coeff_global_full^2,
             assortativity_ethnic_NMSE = (assortativity_ethnic_full - assortativity_ethnic)^2 / 
               assortativity_ethnic_full^2,
             sample_depth = sample_depth,
             sample_method_bias = sample_method_bias
      ) %>%
      select(sample_method_bias,sample_depth,avg_degree_NMSE,density_NMSE,avg_local_clustering_coeff_NMSE,clustering_coeff_global_full_NMSE,assortativity_ethnic_NMSE)
    
    RMSE_df <- rbind(RMSE_df,samples_statistics)
    
    cat("run complete", sample_method_bias)
    
  }}

#when table retrieved then group by sampling method and calculate average

NMSE_df_summary <- RMSE_df %>% 
  group_by(sample_method_bias,sample_depth) %>%
  summarise(avg_degree_NMSE_av = mean(avg_degree_NMSE),
            avg_degree_NMSE_se = sd(avg_degree_NMSE) / sqrt(n()),
            density_NMSE_av = mean(density_NMSE),
            density_NMSE_se = sd(density_NMSE) / sqrt(n()),
            avg_local_clustering_coeff_NMSE_av = mean(avg_local_clustering_coeff_NMSE),
            avg_local_clustering_coeff_NMSE_se = sd(avg_local_clustering_coeff_NMSE) / sqrt(n()),
            clustering_coeff_global_full_NMSE_av = mean(clustering_coeff_global_full_NMSE),
            clustering_coeff_global_full_NMSE_se = sd(clustering_coeff_global_full_NMSE) / sqrt(n()),
            assortativity_ethnic_NMSE_av = mean(assortativity_ethnic_NMSE),
            assortativity_ethnic_NMSE_se = sd(assortativity_ethnic_NMSE) / sqrt(n()),
            N_samples = n())


write.csv(NMSE_df_summary, "NMSE_statistics_biases.csv")


## Step 15: NMSE reweighted estimates

### attributes ###

attributes_samples <- read.csv("FILENAME")

NMSE_df_summary <- as.data.frame(attributes_samples) %>%
  mutate(
    #gender
    gender_NMSE_naive_mean = ((share_male_full - male_naive_mean)^2 + (share_female_full - female_naive_mean)^2 + (share_gender_unknown_full - nogender_naive_mean)^2) /
      (share_male_full^2 + share_female_full^2 + share_gender_unknown_full^2),
    gender_NMSE_hajek_sampindegree = ((share_male_full - male_hajek_sampindegree)^2 + (share_female_full - female_hajek_sampindegree)^2 + (share_gender_unknown_full - nogender_hajek_sampindegree)^2) /
      (share_male_full^2 + share_female_full^2 + share_gender_unknown_full^2),
    gender_NMSE_hajek_realindegree = ((share_male_full - male_hajek_realindegree)^2 + (share_female_full - female_hajek_realindegree)^2 + (share_gender_unknown_full - nogender_hajek_realindegree)^2) /
      (share_male_full^2 + share_female_full^2 + share_gender_unknown_full^2),
    #ethnicity
    ethnicity_NMSE_naive_mean = ((share_caucasian_full - caucasian_naive_mean)^2 + (share_asian_full - asian_naive_mean)^2 + (share_afram_full - afram_naive_mean)^2+
                                   (share_hispanic_full - hispanic_naive_mean)^2 + (share_other_full - other_race_naive_mean)^2 + (share_unknown_full - unknown_race_naive_mean)^2)/
      (share_caucasian_full^2 + share_asian_full^2 + share_afram_full^2 + share_hispanic_full^2 + share_other_full^2 + share_unknown_full^2),
    ethnicity_NMSE_hajek_sampindegree = ((share_caucasian_full - caucasian_hajek_sampindegree)^2 + (share_asian_full - asian_hajek_sampindegree)^2 + (share_afram_full - afram_hajek_sampindegree)^2+
                                           (share_hispanic_full - hispanic_hajek_sampindegree)^2 + (share_other_full - other_race_hajek_sampindegree)^2 + (share_unknown_full - unknown_race_hajek_sampindegree)^2)/
      (share_caucasian_full^2 + share_asian_full^2 + share_afram_full^2 + share_hispanic_full^2 + share_other_full^2 + share_unknown_full^2),
    ethnicity_NMSE_hajek_realindegree = ((share_caucasian_full - caucasian_hajek_realindegree)^2 + (share_asian_full - asian_hajek_realindegree)^2 + (share_afram_full - afram_hajek_realindegree)^2+
                                           (share_hispanic_full - hispanic_hajek_sampindegree)^2 + (share_other_full - other_race_hajek_realindegree)^2 + (share_unknown_full - unknown_race_hajek_realindegree)^2)/
      (share_caucasian_full^2 + share_asian_full^2 + share_afram_full^2 + share_hispanic_full^2 + share_other_full^2 + share_unknown_full^2),
    #party
    party_NMSE_naive_mean = ((share_democrat_full - democrat_naive_mean)^2 + (share_republican_full - republican_naive_mean)^2 + (share_independent_full - independent_naive_mean)^2) /
      (share_democrat_full^2 + share_republican_full^2 + share_independent_full^2),
    party_NMSE_hajek_sampindegree = ((share_democrat_full - democrat_hajek_sampindegree)^2 + (share_republican_full - republican_hajek_sampindegree)^2 + (share_independent_full - independent_hajek_sampindegree)^2) /
      (share_democrat_full^2 + share_republican_full^2 + share_independent_full^2),
    party_NMSE_hajek_realindegree = ((share_democrat_full - democrat_hajek_realindegree)^2 + (share_republican_full - republican_hajek_realindegree)^2 + (share_independent_full - independent_hajek_realindegree)^2) /
      (share_democrat_full^2 + share_republican_full^2 + share_independent_full^2)
  ) %>%
  select(sample_method,sample_depth,
         gender_NMSE_naive_mean,gender_NMSE_hajek_sampindegree,gender_NMSE_hajek_realindegree,
         ethnicity_NMSE_naive_mean,ethnicity_NMSE_hajek_sampindegree,ethnicity_NMSE_hajek_realindegree,
         party_NMSE_naive_mean,party_NMSE_hajek_sampindegree,party_NMSE_hajek_realindegree) %>%
  group_by(sample_method,sample_depth) %>%
  summarise(
    #gender
    gender_NMSE_naive_mean_av = mean(gender_NMSE_naive_mean),
    gender_NMSE_naive_mean_se = sd(gender_NMSE_naive_mean) / sqrt(n()),
    gender_NMSE_hajek_sampindegree_av = mean(gender_NMSE_hajek_sampindegree),
    gender_NMSE_hajek_sampindegree_se = sd(gender_NMSE_hajek_sampindegree) / sqrt(n()),
    gender_NMSE_hajek_realindegree_av = mean(gender_NMSE_hajek_realindegree),
    gender_NMSE_hajek_realindegree_se = sd(gender_NMSE_hajek_realindegree) / sqrt(n()),
    #ethnicity
    ethnicity_NMSE_naive_mean_av = mean(ethnicity_NMSE_naive_mean),
    ethnicity_NMSE_naive_mean_se = sd(ethnicity_NMSE_naive_mean) / sqrt(n()),
    ethnicity_NMSE_hajek_sampindegree_av = mean(ethnicity_NMSE_hajek_sampindegree),
    ethnicity_NMSE_hajek_sampindegree_se = sd(ethnicity_NMSE_hajek_sampindegree) / sqrt(n()),
    ethnicity_NMSE_hajek_realindegree_av = mean(ethnicity_NMSE_hajek_realindegree),
    ethnicity_NMSE_hajek_realindegree_se = sd(ethnicity_NMSE_hajek_realindegree) / sqrt(n()),
    #party
    party_NMSE_naive_mean_av = mean(party_NMSE_naive_mean),
    party_NMSE_naive_mean_se = sd(party_NMSE_naive_mean) / sqrt(n()),
    party_NMSE_hajek_sampindegree_av = mean(party_NMSE_hajek_sampindegree),
    party_NMSE_hajek_sampindegree_se = sd(party_NMSE_hajek_sampindegree) / sqrt(n()),
    party_NMSE_hajek_realindegree_av = mean(party_NMSE_hajek_realindegree),
    party_NMSE_hajek_realindegree_se = sd(party_NMSE_hajek_realindegree) / sqrt(n()),
    #samples
    N_samples = n()) %>%
  select(sample_method,sample_depth,
         gender_NMSE_naive_mean_av,gender_NMSE_naive_mean_se,gender_NMSE_hajek_sampindegree_av,gender_NMSE_hajek_sampindegree_se,
         gender_NMSE_hajek_realindegree_av,gender_NMSE_hajek_realindegree_se,
         ethnicity_NMSE_naive_mean_av,ethnicity_NMSE_naive_mean_se,ethnicity_NMSE_hajek_sampindegree_av,ethnicity_NMSE_hajek_sampindegree_se,
         ethnicity_NMSE_hajek_realindegree_av,ethnicity_NMSE_hajek_realindegree_se,
         party_NMSE_naive_mean_av,party_NMSE_naive_mean_se,party_NMSE_hajek_sampindegree_av,party_NMSE_hajek_sampindegree_se,
         party_NMSE_hajek_realindegree_av,party_NMSE_hajek_realindegree_se,
         N_samples)


write.csv(NMSE_df_summary, "NMSE_attributes_biases_reweighted.csv")


#### global network statistics ####

full_graph_statistics <- read.csv("whole_network_descriptives_in_voterfile.csv", sep = ",")

#full graph statistics
avg_degree_full <- full_graph_statistics$avg_indegree
density_full <- full_graph_statistics$density
avg_local_clustering_coeff_full <- full_graph_statistics$clustering_coeff
clustering_coeff_global_full <- full_graph_statistics$clustering_coeff_global_vector
assortativity_ethnic_full <- full_graph_statistics$assortativity_ethnic

statistics_samples <- read.csv("FILENAME") 

NMSE_df_summary <- as.data.frame(statistics_samples) %>%
  mutate(
    #local clustering coef
    lcl_coef_NMSE_naive = (avg_local_clustering_coeff_full - lcl_coef_naive_mean)^2 / 
      avg_local_clustering_coeff_full^2, ,
    lcl_coef_NMSE_hajek_sampindegree = (avg_local_clustering_coeff_full - lcl_coef_hajek_sampindegree)^2 / 
      avg_local_clustering_coeff_full^2, 
    lcl_coef_NMSE_hajek_realindegree = (avg_local_clustering_coeff_full - lcl_coef_hajek_realindegree)^2 / 
      avg_local_clustering_coeff_full^2,
    #density
    density_NMSE_naive = (density_full - density_naive)^2 / 
      density_full^2,
    density_NMSE_hajek_sampindegree = (density_full - density_hajek_sampindegree)^2 / 
      density_full^2,
    density_NMSE_hajek_realindegree = (density_full - density_hajek_realindegree)^2 / 
      density_full^2,
    #ethnic assortativity
    assortativity_ethnic_NMSE_naive = (assortativity_ethnic_full - assortativity_naive)^2 / 
      assortativity_ethnic_full^2,
    assortativity_ethnic_NMSE_hajek_sampindegree = (assortativity_ethnic_full - assortativity_hajek_sampindegree)^2 / 
      assortativity_ethnic_full^2,
    assortativity_ethnic_NMSE_hajek_realindegree = (assortativity_ethnic_full - assortativity_hajek_realindegree)^2 / 
      assortativity_ethnic_full^2,
    #global clustering coefficient
    clustering_coeff_global_NMSE_naive = (clustering_coeff_global_full - global_clc_hajek_naive)^2 / 
      clustering_coeff_global_full^2,
    clustering_coeff_global_hajek_sampindegree = (clustering_coeff_global_full - global_clc_hajek_sampindegree)^2 / 
      clustering_coeff_global_full^2,
    clustering_coeff_global_NMSE_hajek_realindegree = (clustering_coeff_global_full - global_clc_hajek_realindegree)^2 / 
      clustering_coeff_global_full^2
  ) %>%
  select(sample_method,sample_depth,
         lcl_coef_NMSE_naive,lcl_coef_NMSE_hajek_sampindegree,lcl_coef_NMSE_hajek_realindegree,
         density_NMSE_naive,density_NMSE_hajek_sampindegree,density_NMSE_hajek_realindegree,
         assortativity_ethnic_NMSE_naive,assortativity_ethnic_NMSE_hajek_sampindegree,
         assortativity_ethnic_NMSE_hajek_realindegree,
         clustering_coeff_global_NMSE_naive,clustering_coeff_global_hajek_sampindegree,
         clustering_coeff_global_NMSE_hajek_realindegree) %>%
  group_by(sample_method,sample_depth) %>%
  summarise(# local clc
    lcl_coef_NMSE_naive_mean_av = mean(lcl_coef_NMSE_naive),
    lcl_coef_NMSE_naive_mean_se = sd(lcl_coef_NMSE_naive) / sqrt(n()),
    lcl_coef_NMSE_hajek_sampindegree_av = mean(lcl_coef_NMSE_hajek_sampindegree),
    lcl_coef_NMSE_hajek_sampindegree_se = sd(lcl_coef_NMSE_hajek_sampindegree) / sqrt(n()),
    lcl_coef_NMSE_hajek_realindegree_av = mean(lcl_coef_NMSE_hajek_realindegree),
    lcl_coef_NMSE_hajek_realindegree_se = sd(lcl_coef_NMSE_hajek_realindegree) / sqrt(n()),
    # density
    density_NMSE_naive_mean_av = mean(density_NMSE_naive),
    density_NMSE_naive_mean_se = sd(density_NMSE_naive) / sqrt(n()),
    density_NMSE_hajek_sampindegree_av = mean(density_NMSE_hajek_sampindegree),
    density_NMSE_hajek_sampindegree_se = sd(density_NMSE_hajek_sampindegree) / sqrt(n()),
    density_NMSE_hajek_realindegree_av = mean(density_NMSE_hajek_realindegree),
    density_NMSE_hajek_realindegree_se = sd(density_NMSE_hajek_realindegree) / sqrt(n()),
    # ethnic assortativity
    assortativity_ethnic_NMSE_naive_mean_av = mean(assortativity_ethnic_NMSE_naive),
    assortativity_ethnic_NMSE_naive_mean_se = sd(assortativity_ethnic_NMSE_naive) / sqrt(n()),
    assortativity_ethnic_NMSE_hajek_sampindegree_av = mean(assortativity_ethnic_NMSE_hajek_sampindegree),
    assortativity_ethnic_NMSE_hajek_sampindegree_se = sd(assortativity_ethnic_NMSE_hajek_sampindegree) / sqrt(n()),
    assortativity_ethnic_NMSE_hajek_realindegree_av = mean(assortativity_ethnic_NMSE_hajek_realindegree),
    assortativity_ethnic_NMSE_hajek_realindegree_se = sd(assortativity_ethnic_NMSE_hajek_realindegree) / sqrt(n()),
    # global clc
    clustering_coeff_global_NMSE_naive_mean_av = mean(clustering_coeff_global_NMSE_naive),
    clustering_coeff_global_NMSE_naive_mean_se = sd(clustering_coeff_global_NMSE_naive) / sqrt(n()),
    clustering_coeff_global_NMSE_hajek_sampindegree_av = mean(clustering_coeff_global_hajek_sampindegree),
    clustering_coeff_global_NMSE_hajek_sampindegree_se = sd(clustering_coeff_global_hajek_sampindegree) / sqrt(n()),
    clustering_coeff_global_NMSE_hajek_realindegree_av = mean(clustering_coeff_global_NMSE_hajek_realindegree),
    clustering_coeff_global_NMSE_hajek_realindegree_se = sd(clustering_coeff_global_NMSE_hajek_realindegree) / sqrt(n()),
    #samples
    N_samples = n()) %>%
  select(sample_method,sample_depth,
         lcl_coef_NMSE_naive_mean_av,lcl_coef_NMSE_naive_mean_se,lcl_coef_NMSE_hajek_sampindegree_av,lcl_coef_NMSE_hajek_sampindegree_se,
         lcl_coef_NMSE_hajek_realindegree_av,lcl_coef_NMSE_hajek_realindegree_se,
         density_NMSE_naive_mean_av,density_NMSE_naive_mean_se,density_NMSE_hajek_sampindegree_av,density_NMSE_hajek_sampindegree_se,
         density_NMSE_hajek_realindegree_av,density_NMSE_hajek_realindegree_se,
         assortativity_ethnic_NMSE_naive_mean_av,assortativity_ethnic_NMSE_naive_mean_se,assortativity_ethnic_NMSE_hajek_sampindegree_av,assortativity_ethnic_NMSE_hajek_sampindegree_se,
         assortativity_ethnic_NMSE_hajek_realindegree_av,assortativity_ethnic_NMSE_hajek_realindegree_se,
         clustering_coeff_global_NMSE_naive_mean_av, clustering_coeff_global_NMSE_naive_mean_se,
         clustering_coeff_global_NMSE_hajek_sampindegree_av,clustering_coeff_global_NMSE_hajek_sampindegree_se,
         clustering_coeff_global_NMSE_hajek_realindegree_av,clustering_coeff_global_NMSE_hajek_realindegree_se,
         N_samples)

write.csv(NMSE_df_summary, "NMSE_statistics_biases_reweighted.csv")


## Step 16: Indegree inclusion share correlation for Appendix


#get real indegrees for each vertex in complete graph
full_df <- data.frame(name_vertex = as.character(V(graph_full_inner)$useridchr),
                      indegree_full = degree(graph_full_inner, mode = "in"))
rm(graph_full_inner)

# exclusion of indegree bias for this purpose
bias_types <- c("unbiased","partybias/democrat5percentmore","partybias/democrat10percentmore","partybias/democrat100percentmore",
                "partybias/republican5percentmore","partybias/republican10percentmore","partybias/republican100percentmore",
                "genderbias/women5percentmore","genderbias/women10percentmore","genderbias/women100percentmore",
                "genderbias/men5percentmore","genderbias/men10percentmore","genderbias/men100percentmore",
                "ethnicity/caucasian5percentmore","ethnicity/caucasian10percentmore","ethnicity/caucasian100percentmore",
                "ethnicity/asian5percentmore","ethnicity/asian10percentmore","ethnicity/asian100percentmore") # all the biased and unbiased samples except indegree bias combined (for power)


# k1 

for (method in bias_types) {
  
  if (method == "unbiased") { #load unbiased
    load("S100N564k1network_sample_results_BFS.rdata")
    
    num_results <- length(results_bfs)
    
    for (sample in 1:num_results){
      subgraph <- results_bfs[[sample]]$subgraph
      df_sample <- data.frame(name_vertex = as.character(V(subgraph)$useridchr),
                              method = 1)
      colnames(df_sample) <- c("name_vertex",paste0(sub(".*/", "", method),"_",sample))                        
      
      full_df <- left_join(full_df, df_sample, by = "name_vertex")
      
    }
  } else{
    for (run in 1:2){
      load("FILENAME.rdata") # load all bias types
      
      num_results <- length(results_bfs)
      
      
      for (sample in 1:num_results){
        subgraph <- results_bfs[[sample]]$subgraph
        df_sample <- data.frame(name_vertex = as.character(V(subgraph)$useridchr),
                                method = 1)
        
        if (run == 2){sample <- sample + 50}                        
        colnames(df_sample) <- c("name_vertex",paste0(sub(".*/", "", method),"_",sample))                        
        
        full_df <- left_join(full_df, df_sample, by = "name_vertex")
      }}
    
  }}

write.csv(full_df, "k1_indegree_inclusion.csv")

full_df <- full_df %>%
  mutate(across(where(is.numeric), ~coalesce(.x, 0)),
         share_samples = rowSums(pick(-c(name_vertex, indegree_full)), na.rm = TRUE) /(ncol(.)-2)) %>%
  select(name_vertex,share_samples,indegree_full)

cor(full_df$indegree_full, full_df$share_samples) # 0.7045

png(filename = "inclusion_indegree_k1.png",
    width = 30, height = 30, res = 300, units = "cm")

ggplot(full_df, aes(x = share_samples, y = indegree_full)) +
  geom_point()+
  geom_smooth(
    se = TRUE,
    color = "#D20064",
    fill = "#F4A3C4",
    linewidth = 1.2
  ) +
  labs(
    x = "Inclusion probability",
    y = "Indegree in complete network"
  ) +
  ylim(0,15000)+
  theme_minimal(base_size = 14)

dev.off()


## Step 17: Plots main paper: Attributes

# complete graph data

full_202109_attributes <- read.csv("whole_network_attributes.csv", sep = ",")

full_202109_attributes_long <- full_202109_attributes %>% 
  select(-X, -sample_size,-sample_depth, -sample_method,-sample_num) %>% 
  pivot_longer(cols = everything(), names_to = "attribute", values_to = "full_network") %>% 
  filter(!(attribute %in% c("gender_unknown","race_unknown","race_other")))
# read and combine 

input <- "REPOSITORY_PATH"

bias_types <- c("biasdegree_linear","partybias_democrat5percentmore","partybias_democrat10percentmore","partybias_democrat100percentmore",
                "partybias_republican5percentmore","partybias_republican10percentmore","partybias_republican100percentmore",
                "genderbias_women5percentmore","genderbias_women10percentmore","genderbias_women100percentmore",
                "genderbias_men5percentmore","genderbias_men10percentmore","genderbias_men100percentmore",
                "ethnicity_caucasian5percentmore","ethnicity_caucasian10percentmore","ethnicity_caucasian100percentmore",
                "ethnicity_asian5percentmore","ethnicity_asian10percentmore","ethnicity_asian100percentmore") #length: 19

#real donation
donation_k1 <- read.csv(paste0(input,"FILENAME_attributes_BFS_donated.csv")) %>% 
  mutate(sample_method = "real donation",
         size = "k1")

# simulated donations
samples_k1 <- read.csv(paste0(input, "FILENAME_attributes_BFS.csv")) %>% 
  mutate(seed_selection = "no bias")

for (bias in 1:length(bias_types)) {
  
  bias_type <- bias_types[bias]
  
  sample <- read.csv(paste0(input, "FILENAME_attributes_BFS_",bias_type,".csv")) %>% 
    mutate(seed_selection = bias_type)
  
  samples_k1 <- rbind(samples_k1,sample)
  
  rm(sample)
}

samples_k1 <- samples_k1 %>% 
  select(-X)

## random walk

rwj015_k1 <- read.csv(paste0(input, "FILENAME_attributes_RW.csv")) %>% 
  mutate(seed_selection = "RW j = 0.15",
         bias_category = "RW",
         bias_explicit = "RW j = 0.15")


rw_k1 <- rwj015_k1 %>% 
  rename(sample_number = sample_num) %>% 
  select(-X, -sample_depth, -sample_method) %>% 
  pivot_longer(cols = !c(seed_selection, bias_category, bias_explicit, sample_number,sample_size), names_to = "attribute", values_to = "value") %>% 
  left_join(full_202109_attributes_long, by = "attribute")


samples_k1_long <- samples_k1 %>% 
  rename(sample_number = sample_num) %>% 
  mutate(bias_category = sub("_.*", "", seed_selection),
         bias_explicit = sub(".*_", "", seed_selection)) %>% 
  select(-sample_depth, -sample_method) %>% 
  pivot_longer(cols = !c(seed_selection, bias_category, bias_explicit, sample_number,sample_size), names_to = "attribute", values_to = "value") %>% 
  left_join(full_202109_attributes_long, by = "attribute") %>% 
  rbind(.,rw_k1) %>%
  filter(!(attribute %in% c("gender_unknown","race_unknown","race_other"))) %>% 
  mutate(share = value / sample_size,
         share_full = full_network / full_202109_attributes$sample_size)

samples_k1_long$bias_explicit[samples_k1_long$bias_explicit == "linear"] <- "popularity bias (indegree linear)"


full_202109_attributes_percent <- full_202109_attributes_long %>% 
  mutate(full_network = full_network /full_202109_attributes$sample_size) %>% 
  pivot_wider(names_from = attribute, values_from = full_network)

options(
  scipen = 999,
  digits = 4)

samples_k1_long <- samples_k1_long %>% 
  mutate(bias_explicit = factor(bias_explicit, levels = unique(bias_explicit)),
         bias_category = factor(bias_category, levels = unique(bias_category)),
         attribute = factor(attribute, levels = unique(attribute)))

scale_color <- rev(c("#0072B2","#E69F00","#56B4E9","#F0E442","#009E73","#CC79A7")) #colorblind palette
color_donated <- "black" 
color_line <- "#D55E00"

samples_k1_long_all <- samples_k1_long %>% 
  filter(attribute %in% c("female","male","caucasian","asian","hispanic","african_american","democrat","republican")) %>% 
  left_join(., donation_k1 %>% select(sample_size,female,male, caucasian, asian,african_american, hispanic, democrat, republican) %>% 
              pivot_longer(!sample_size, names_to = "attribute", values_to = "num") %>%
              mutate(share_donated = num/ sample_size) %>% select(attribute, share_donated), by = "attribute") %>% 
  mutate(share_donated = if_else(bias_category != "RW",share_donated,NA),
         attribute = case_match(attribute,
                                "female" ~ "Female",
                                "male" ~ "Male",
                                "caucasian" ~ "Caucasian",
                                "asian" ~ "Asian",
                                "hispanic" ~ "Hispanic",
                                "african_american" ~ "African American",
                                "democrat" ~ "Democrat",
                                "republican" ~ "Republican"))

# all attributes
plot_k1_all <- ggplot(samples_k1_long_all, aes(y = bias_explicit)) +
  geom_violin(aes(x = share, fill = bias_category)) +
  scale_y_discrete(limits = rev) +
  geom_hline(yintercept = 1.5, color = "black") +
  geom_vline(
    aes(xintercept = share_full, color = "hline"),
    linewidth = 1.2
  ) + 
  geom_point(
    aes(x = share_donated, color = "point"),
    fill = color_donated,
    size = 4.5,
    shape = 18,
    inherit.aes = T
  ) +
  scale_color_manual(
    name = "Reference",
    values = c(
      "hline" = color_line,
      "point" = color_donated
    ),
    labels = c(
      "hline" = "Original Network",
      "point" = "BFS, Donated Seeds"
    )
  ) +
  scale_x_continuous(labels = percent_format()) +
  scale_fill_manual(
    "Sampling Method\nand Seed\nSelection Bias",
    labels = c(
      "BFS, No Bias",
      "BFS, Popularity Bias\n(Indegree)",
      "BFS, Party Bias",
      "BFS, Gender Bias",
      "BFS, Racial Bias",
      "Random Walk"
    ),
    values = scale_color
  ) +
  xlab("Share Attribute") +
  ylab("Sampling Bias") +
  facet_wrap2(
    . ~ attribute,
    scales = "free_x",
    strip = strip_themed(
      background_x = elem_list_rect(
        fill = c(
          "#CDB2D6", "#CDB2D6",
          "#C1DFDF", "#C1DFDF", "#C1DFDF", "#C1DFDF",
          "#B2D8FF", "#B2D8FF"
        ) #colors per facet level
      )
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    strip.text = element_text(size = 15),
    legend.title = element_text(size = 20)
  )

atts <- unique(samples_k1_long_all$attribute)

plots <- lapply(seq_along(atts), function(i){
  
  att <- atts[i]
  
  p <-
    plot_k1_all %+%
    filter(samples_k1_long_all, attribute == att) +
    facet_null() +
    labs(x = att) +
    theme(
      legend.position = "none",
    )
  
  ## remove y axis from middle/right columns
  if(!(i %in% c(1,4,7))){
    p <- p +
      theme(
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  p
})

# Extract legend
legend <- wrap_elements(
  cowplot::get_legend(
    plot_k1_all +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)
      )
  )
)

# Assemble
plot_attributes_biases_all <- wrap_plots(
  c(plots, list(legend)),
  ncol = 3,
  byrow = TRUE
)

setEPS()
postscript(
  "k1_compare_biases_attributes_all.eps",
  width = 33 / 2.54,
  height = 45 / 2.54,
  #units = "cm",
  paper = "special"
)

plot_attributes_biases_all

dev.off()


## reweighting plots

attributes_reweighted <- read.csv("attributes_weights.csv", sep = ";")


attributes_reweighted_long <- attributes_reweighted %>%
  select(-starts_with("X")) %>%
  filter(sample_depth == "k1") %>% 
  pivot_longer(
    cols = -c(sample_num,sample_method,sample_depth,sample_size),
    names_to = c("attribute", "statistic"),
    names_pattern = "([^_]+)_(.*)"
  ) %>%
  pivot_wider(names_from = statistic, values_from = value) %>% 
  mutate(attribute = case_match(attribute,
                                "afram" ~ "african_american", .default = attribute)) %>% 
  left_join(., full_202109_attributes_long %>%  mutate(full_network = full_network / 911331), by = "attribute") %>% 
  select(!starts_with("race_")) %>% 
  pivot_longer(c("naive_mean", starts_with("hajek")),
               names_to = "weighting",
               values_to = "value") 

attributes_reweighted_long_donated <- attributes_reweighted_long %>% 
  filter(sample_method == "donated") 

attributes_reweighted_long <- attributes_reweighted_long %>%
  filter(sample_method != "donated") %>%
  left_join(.,attributes_reweighted_long_donated %>% select(donated = value,attribute,weighting), by = c("attribute","weighting")) %>% 
  filter(attribute %in% c("female","male","caucasian","asian","hispanic","african_american","democrat","republican")) %>%
  mutate(donated = as.numeric(gsub(",", ".", gsub("\\.", "", donated))),
         value = as.numeric(gsub(",", ".", gsub("\\.", "", value))),
         sample_method = case_match(sample_method,
                                    "RW0.15" ~ "Random Walk",
                                    "BFS_biasdegree_linear" ~ "BFS, Popularity Bias",
                                    "BFS" ~ "BFS, No Bias"),
         weighting = case_match(weighting,
                                "naive_mean" ~ "unweighted",
                                "hajek_sampindegree" ~ "weighted by sample indegree",
                                "hajek_realindegree" ~ "weighted by real indegree"),
         weighting = factor(weighting, levels = rev(c("unweighted","weighted by sample indegree","weighted by real indegree"))),
         sample_method = factor(sample_method, levels = unique(sample_method)),
         donated = if_else(sample_method != "Random Walk",donated,NA),
         attribute = case_match(attribute,
                                "female" ~ "Female",
                                "male" ~ "Male",
                                "caucasian" ~ "Caucasian",
                                "asian" ~ "Asian",
                                "hispanic" ~ "Hispanic",
                                "african_american" ~ "African American",
                                "democrat" ~ "Democrat",
                                "republican" ~ "Republican"),
         attribute = factor(attribute, levels = c("Female","Male","Democrat","Republican","Caucasian","Asian","Hispanic","African American")))


plot_k1_attributes_reweighting_part1 <- ggplot(attributes_reweighted_long %>% filter(attribute %in% c("Female","Male","Democrat","Republican")),
                                               aes(y = interaction(weighting)))+
  geom_violin(aes(x = value, fill = sample_method))+
  geom_vline(aes(xintercept = full_network, color = "hline"), linewidth = 1.2)+
  geom_point(aes(x = donated, color = "point"),fill = color_donated, size = 4.5, shape = 18,
             inherit.aes = T) +
  scale_color_manual(
    name = "",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Original Network", "point" = "BFS, Donated Seeds")
  ) +
  scale_x_continuous(labels = percent_format()) +
  scale_fill_manual("", values = c("#CC79A7","#009E73","#0072B2"))+
  xlab("")+
  ylab("Reweighting")+
  facet_grid2(sample_method~attribute, scales = "free_x", switch = "x",
              strip = strip_themed(background_x = elem_list_rect(fill = c("white","white","white"), color = c("white","white","white")),
                                   background_y = elem_list_rect(fill = c("white","white","white"), color = c("#CC79A7","#009E73","#0072B2"))))+
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size=20),strip.placement = "outside", strip.text = element_text(size = 15), legend.title = element_text(size=24))

plot_k1_attributes_reweighting_part2 <- ggplot(attributes_reweighted_long %>% filter(attribute %in% c("Caucasian","Asian","Hispanic","African American")),
                                               aes(y = interaction(weighting)))+ 
  geom_violin(aes(x = value, fill = sample_method))+
  geom_vline(aes(xintercept = full_network, color = "hline"), linewidth = 1.2)+
  geom_point(aes(x = donated, color = "point"),fill = color_donated, size = 4.5, shape = 18,
             inherit.aes = T) +
  scale_color_manual(
    name = "",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Original Network", "point" = "BFS, Donated Seeds")
  ) +
  scale_x_continuous(labels = percent_format()) +
  scale_fill_manual("", values = c("#CC79A7","#009E73","#0072B2"))+
  xlab("")+
  ylab("Reweighting")+
  facet_grid2(sample_method~attribute, scales = "free_x", switch = "x",
              strip = strip_themed(background_x = elem_list_rect(fill = c("white","white","white"), color = c("white","white","white")),
                                   background_y = elem_list_rect(fill = c("white","white","white"), color = c("#CC79A7","#009E73","#0072B2"))))+
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size=20),strip.placement = "outside", strip.text = element_text(size = 15), legend.title = element_text(size=24))


png(filename = "k1_attributes_reweighting.png", width = 39, height = 43, units = "cm", res = 300)

(plot_k1_attributes_reweighting_part1 / plot_k1_attributes_reweighting_part2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical")

dev.off()

setEPS()
postscript(
  "k1_attributes_reweighting.eps",
  width = 40 / 2.54,
  height = 43 / 2.54,
  #units = "cm",
  paper = "special"
)

(plot_k1_attributes_reweighting_part1 / plot_k1_attributes_reweighting_part2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "vertical")

dev.off()


## Step 18: Plots main paper: Statistics


# complete graph data

full_202109_stats <- read.csv("whole_network_descriptives_in_voterfile.csv", sep = ",")

full_202109_stats_long <- full_202109_stats%>% 
  select(avg_indegree, density, clustering_coeff, assortativity_ethnic, clustering_coeff_global_vector) %>% 
  pivot_longer(cols = c(avg_indegree, density, clustering_coeff, assortativity_ethnic, clustering_coeff_global_vector), names_to = "statistic", values_to = "full_network")



# read and combine 

input <- "REPOSITORY"

bias_types <- c("biasdegree_linear","partybias_democrat5percentmore","partybias_democrat10percentmore","partybias_democrat100percentmore",
                "partybias_republican5percentmore","partybias_republican10percentmore","partybias_republican100percentmore",
                "genderbias_women5percentmore","genderbias_women10percentmore","genderbias_women100percentmore",
                "genderbias_men5percentmore","genderbias_men10percentmore","genderbias_men100percentmore",
                "ethnicity_caucasian5percentmore","ethnicity_caucasian10percentmore","ethnicity_caucasian100percentmore",
                "ethnicity_asian5percentmore","ethnicity_asian10percentmore","ethnicity_asian100percentmore") #length: 19

# k1
#real donation
donation_k1 <- read.csv("S3N564k1network_sample_descriptives_BFS_donated.csv") %>% 
  mutate(sample_method = "real donation",
         size = "k1")

# simulated donations
samples_k1 <- read.csv(paste0(input, "S100N564k1network_sample_descriptives_BFS.csv")) %>% 
  mutate(seed_selection = "no bias")

for (bias in 1:length(bias_types)) {
  
  bias_type <- bias_types[bias]
  
  sample <- read.csv(paste0(input, "S100N564k1network_sample_descriptives_BFS_",bias_type,".csv")) %>% 
    mutate(seed_selection = bias_type)
  
  samples_k1 <- rbind(samples_k1,sample)
  
  rm(sample)
}
samples_k1 <- samples_k1 %>% 
  select(-X)

## random walk

rwj015_k1 <- read.csv(paste0(input, "S100N10373k1j0.15network_sample_descriptives_RW.csv")) %>% 
  mutate(seed_selection = "RW j = 0.15",
         bias_category = "RW",
         bias_explicit = "RW j = 0.15")

rw_k1 <-rbind(rwj015_k1) %>% 
  rename(sample_number = sample_num) %>% 
  select(seed_selection,bias_category,bias_explicit,sample_number,
         avg_indegree,density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic) %>% 
  pivot_longer(cols = c(avg_indegree, density, clustering_coeff, assortativity_ethnic, clustering_coeff_global_vector), names_to = "statistic", values_to = "value") %>% 
  left_join(full_202109_stats_long, by = "statistic")


samples_k1_long <- samples_k1 %>% 
  mutate(sample_number = X.1,
         bias_category = sub("_.*", "", seed_selection),
         bias_explicit = sub(".*_", "", seed_selection)) %>% 
  select(seed_selection,bias_category,bias_explicit,sample_number,
         avg_indegree,density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic) %>% 
  pivot_longer(cols = c(avg_indegree, density, clustering_coeff, assortativity_ethnic, clustering_coeff_global_vector), names_to = "statistic", values_to = "value") %>% 
  left_join(full_202109_stats_long, by = "statistic") %>% 
  rbind(.,rw_k1) %>% 
  mutate(bias_explicit = case_match(bias_explicit,
                                    "linear" ~ "popularity bias (indegree linear)", .default = bias_explicit),
         bias_explicit = factor(bias_explicit, levels = rev(unique(bias_explicit))),
         bias_category = factor(bias_category, levels = rev(unique(bias_category))),
         statistic = factor(statistic, levels = unique(statistic))) 

levels_bias <- unique(samples_k1_long$bias_explicit)
levels_bias_category <- unique(samples_k1_long$bias_category)  

write.csv(samples_k1_long, "statistics_samples_k1.csv")


scale_color <- rev(c("#0072B2","#E69F00","#56B4E9","#F0E442","#009E73","#CC79A7")) #colorblind palette
color_donated <- "black" 
color_line <- "#D55E00" 

options(
  scipen = 999,
  digits = 4)

# subsets of data for each statistic
samples_k1_long_density <- samples_k1_long %>% 
  filter(statistic == "density") %>% 
  mutate(bias_explicit = factor(bias_explicit, levels = levels_bias),
         bias_category = factor(bias_category, levels = levels_bias_category)) %>% 
  mutate(value_full = full_202109_stats$density,
         value_donated = if_else(bias_category != "RW", donation_k1$density, NA))

samples_k1_long_degree <- samples_k1_long %>% 
  filter(statistic == "avg_indegree") %>% 
  mutate(bias_explicit = factor(bias_explicit, levels = levels_bias),
         bias_category = factor(bias_category, levels = levels_bias_category)) %>% 
  mutate(value_full = NA,
         value_donated = if_else(bias_category != "RW", donation_k1$avg_degree, NA))

samples_k1_long_assortativity_ethnic <- samples_k1_long %>% 
  filter(statistic == "assortativity_ethnic") %>% 
  mutate(bias_explicit = factor(bias_explicit, levels = levels_bias),
         bias_category = factor(bias_category, levels = levels_bias_category)) %>% 
  mutate(value_full = full_202109_stats$assortativity_ethnic,
         value_donated = if_else(bias_category != "RW", donation_k1$assortativity_ethnic, NA))

samples_k1_long_clusteringcoeff_global <- samples_k1_long %>% 
  filter(statistic == "clustering_coeff_global_vector") %>% 
  mutate(bias_explicit = factor(bias_explicit, levels = levels_bias),
         bias_category = factor(bias_category, levels = levels_bias_category)) %>% 
  mutate(value_full = full_202109_stats$clustering_coeff_global_vector,
         value_donated = if_else(bias_category != "RW", donation_k1$clustering_coeff_global_vector, NA))

samples_k1_long_clusteringcoeff_local <- samples_k1_long %>% 
  filter(statistic == "clustering_coeff") %>% 
  mutate(bias_explicit = factor(bias_explicit, levels = levels_bias),
         bias_category = factor(bias_category, levels = levels_bias_category))  %>% 
  mutate(value_full = full_202109_stats$clustering_coeff,
         value_donated = if_else(bias_category != "RW", donation_k1$clustering_coeff, NA))


density_plot_k1 <- 
  ggplot(samples_k1_long_density, aes(x = interaction(bias_explicit)))+
  geom_violin(aes(y = value, fill = bias_category))+
  geom_hline(aes(yintercept = full_202109_stats$density[1], color = "hline"), linewidth = 1.2)+
  geom_vline(xintercept = 1.5, color = "black") +
  scale_y_continuous(labels = function(x) sprintf("%7.4f", x)) + #
  scale_x_discrete(limits = rev) +
  coord_flip()+
  geom_point(aes(y = value_donated, color = "point"), size = 3.5, shape = 23,
             fill = color_donated,
             inherit.aes = T) +
  scale_color_manual(
    name = "Reference",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Value full network", "point" = "Donated seeds"),
    breaks = c("hline", "point"),
  ) +
  scale_fill_manual(values = scale_color, name = "Bias category", labels = c("No bias","Popularity Bias (Indegree)","Party Bias", "Gender Bias","Racial Bias","Random Walk"))+
  ylab("Density")+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2)) +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), strip.text = element_text(size = 15), legend.title = element_text(size=24))

degree_plot_k1 <- 
  ggplot(samples_k1_long_degree, aes(x = interaction(bias_explicit)))+
  geom_violin(aes(y = value, fill = bias_category))+
  scale_x_discrete(limits = rev) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_hline(aes(yintercept = full_202109_stats$avg_indegree[1], color = "hline"), linewidth = 1.2)+
  coord_flip()+
  scale_y_continuous(limits = c(0, 21.44), breaks = c(0.00,2.00,4.00,6.00,20.67),labels = function(x) sprintf("%7.4f", x)) + #
  scale_y_break(c(6,20)) +
  ylab("Average Degree")+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2)) +
  geom_point(aes(y = value_donated, color = "point"), size = 3.5, shape = 23,
             fill = color_donated,
             inherit.aes = T) +
  scale_color_manual(
    name = "Reference",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Value full network", "point" = "Donated seeds"),
    breaks = c("hline", "point"),
  ) +
  scale_fill_manual(values = scale_color, name = "Bias category", labels = c("no bias","popularity bias (indegree)","party bias", "gender bias","racial bias","Random Walk"))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), strip.text = element_text(size = 15), legend.title = element_text(size=24),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.title.x.top = element_blank())


assortativity_ethnic_plot_k1 <- 
  ggplot(samples_k1_long_assortativity_ethnic, aes(x = interaction(bias_explicit)))+
  geom_violin(aes(y = value, fill = bias_category))+
  scale_x_discrete(limits = rev) +
  scale_y_continuous(breaks = c(0.0,0.1,0.2,0.3),labels = function(x) sprintf("%6.1f", x))+
  geom_vline(xintercept = 1.5, color = "black") +
  geom_hline(aes(yintercept = full_202109_stats$assortativity_ethnic[1], color = "hline"), linewidth = 1.2) +
  coord_flip()+
  geom_point(aes(y = value_donated, color = "point"),fill = color_donated, size = 4.5, shape = 18,
             inherit.aes = T) +
  scale_color_manual(
    name = "Reference",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Value full network", "point" = "Donated seeds")
  ) +
  scale_fill_manual(values = scale_color, name = "Bias category", labels = c("no bias","popularity bias (indegree)","party bias", "gender bias","racial bias","Random Walk"))+
  ylab("Racial Assortativity")+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2)) +
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), strip.text = element_text(size = 15), legend.title = element_text(size=24))


clusteringcoeff_global_plot_k1 <- 
  ggplot(samples_k1_long_clusteringcoeff_global,
         aes(x = interaction(bias_explicit)))+
  geom_violin(aes(y = value, fill = bias_category))+
  scale_fill_manual(values = scale_color, name = "Bias category",
                    labels = c("no bias","popularity bias (indegree)","party bias", "gender bias","racial bias","Random Walk"))+
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0, 0.09), breaks = c(0.000,0.001, 0.002,0.012,0.06,0.08),labels = function(x) sprintf("%7.3f", x))+
  scale_y_break(c(0.0025,0.011), scales = .7) +
  scale_y_break(c(0.013,0.063), scales = .3) +
  geom_vline(xintercept = 1.5, color = "black") +
  geom_hline(aes(yintercept = full_202109_stats$clustering_coeff_global_vector[1], color = "hline"), linewidth = 1.2)+
  geom_point(aes(y = value_donated, color = "point"), size = 3.5, shape = 23,
             fill = color_donated,
             inherit.aes = T) +
  scale_color_manual(
    name = "Reference",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Value full network", "point" = "Donated seeds")
  ) +
  ylab("Global Clustering Coef.")+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2)) +
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), strip.text = element_text(size = 15), legend.title = element_text(size=24),
        axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.title.x.top = element_blank())

clusteringcoeff_local_plot_k1 <- 
  ggplot(samples_k1_long_clusteringcoeff_local, aes(x = interaction(bias_explicit)))+
  geom_violin(aes(y = value, fill = bias_category))+
  scale_fill_manual(values = scale_color, labels = c("No Bias","Popularity Bias (Indegree)","Party Bias", "Gender Bias","Racial Bias","Random Walk"),
                    name = "Bias Category")+
  scale_x_discrete(limits = rev) +
  scale_y_continuous(limits = c(0, 0.3), breaks = c(0.1,0.2,0.3),labels = function(x) sprintf("%7.1f", x))+
  geom_vline(xintercept = 1.5, color = "black") +
  geom_hline(aes(yintercept = full_202109_stats$clustering_coeff[1], color = "hline"), linewidth = 1.2)+
  geom_point(aes(y = donation_k1$clustering_coeff_global_vector, color = "point"), size = 3.5, shape = 23,
             fill = color_donated,
             inherit.aes = T) +
  scale_color_manual(
    name = "Reference",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Value\rOriginal Network", "point" = "Donated Seeds")
  ) +
  ylab("Local Clustering Coef.")+
  guides(fill = guide_legend(order = 1),
         color = guide_legend(order = 2)) +
  coord_flip()+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), strip.text = element_text(size = 15), legend.title = element_text(size=24))


legend_plot <- ggplot(
  samples_k1_long_clusteringcoeff_local,
  aes(x = interaction(bias_explicit))
) +
  geom_violin(aes(y = value, fill = bias_category)) +
  geom_vline(
    aes(xintercept = 1.5, color = "hline"),
    linewidth = 1.2
  ) +
  geom_point(
    aes(y = donation_k1$clustering_coeff_global_vector, color = "point"),
    size = 3.5,
    shape = 23,
    fill = color_donated
  ) +
  scale_fill_manual(
    values = scale_color,
    labels = c(
      "BFS, No Bias", "BFS, Popularity Bias (Indegree)", "BFS, Party Bias",
      "BFS, Gender Bias", "BFS, Racial Bias", "Random Walk"
    ),
    name = "Sampling Method and\nSeed Selection Bias"
  ) +
  scale_color_manual(
    name = "Reference",
    values = c(
      hline = color_line,
      point = color_donated
    ),
    labels = c(
      hline = "Original Network",
      point = "BFS, Donated Seeds"
    )
  ) +
  theme_void() +
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), strip.text = element_text(size = 15), legend.title = element_text(size=24))

theme_right <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank()
)


degree_plot_k1 <- degree_plot_k1 + theme_right
clusteringcoeff_global_plot_k1 <- clusteringcoeff_global_plot_k1 + theme_right

legend <- get_legend(
  legend_plot +
    guides(fill = guide_legend(nrow = 6),
           color = guide_legend(nrow = 2)) +
    theme(
      legend.position = "right"
    )
)
legend_plot <- ggdraw(legend)

p1 <- density_plot_k1 + theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.title.y = element_blank())
p1 <- ggdraw() +
  draw_plot(p1, x = 0.02, y = 0.02, width = 0.969, height = 0.969)

p2 <- degree_plot_k1 + theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.title.y = element_blank())

p3 <- clusteringcoeff_local_plot_k1 + theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.title.y = element_blank())
p3 <- ggdraw() +
  draw_plot(p3, x = 0.02, y = 0.02, width = 0.969, height = 0.969)

p4 <- clusteringcoeff_global_plot_k1 + theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.title.y = element_blank())

p5 <- assortativity_ethnic_plot_k1 + theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.title.y = element_blank())
p5 <- ggdraw() +
  draw_plot(p5, x = 0.02, y = 0.02, width = 0.969, height = 0.969)


png(filename = "combined_statistics_k1.png", width = 33, height = 45, res = 900, units = "cm")

aplot::plot_list(p1,p2,p3,p4,p5,legend_plot, nrow = 3, widths = c(1.7,1), heights = c(1,1,1))

dev.off()


setEPS()
postscript(
  "combined_statistics_k1.eps",
  width = 33 / 2.54,
  height = 45 / 2.54,
  paper = "special"
)

aplot::plot_list(p1,p2,p3,p4,p5,legend_plot, nrow = 3, widths = c(1.7,1), heights = c(1,1,1))

dev.off()


## plot reweight

statistics_reweighted <- read.csv("statistics_df_weights_noedgeslastlayer_k1.csv", sep = ",") %>% 
  rbind(read.csv("statistics_df_weights_noedgeslastlayer_global_clc_k1.csv", sep = ","))


statistics_reweighted_long <- statistics_reweighted %>%
  select(-starts_with("X"), -sample_size.y) %>%
  filter(sample_depth == "k1") %>% 
  pivot_longer(
    cols = -c(sample_num,sample_method,sample_depth,sample_size.x),
    names_to = c("statistic", "weighting"),
    names_pattern = "^((?:global|lcl)_[^_]+|[^_]+)_(.*)$"
  ) %>%
  left_join(., full_202109_stats_long %>%  mutate(full_network = full_network,
                                                  statistic = case_match(statistic,
                                                                         "assortativity_ethnic" ~"assortativity",
                                                                         "density" ~ "density",
                                                                         "clustering_coeff" ~"lcl_coef",
                                                                         "clustering_coeff_global_vector" ~ "global_clc",
                                                                         .default = statistic)), by = "statistic") 

statistics_reweighted_long_donated <- statistics_reweighted_long %>% 
  filter(sample_method == "donated") 

statistics_reweighted_long <- statistics_reweighted_long %>%
  filter(sample_method != "donated") %>%
  left_join(.,statistics_reweighted_long_donated %>% select(donated = value,statistic,weighting), by = c("statistic","weighting")) %>% 
  mutate(
    sample_method = case_match(sample_method,
                               "RW0.15" ~ "Random Walk",
                               "BFS_biasdegree_linear" ~ "BFS, Popularity Bias",
                               "BFS" ~ "BFS, No Bias"),
    weighting = case_match(weighting,
                           "naive" ~ "unweighted",
                           "naive_mean" ~ "unweighted",
                           "hajek_naive" ~ "unweighted",
                           "hajek_sampindegree" ~ "weighted by sample indegree",
                           "hajek_realindegree" ~ "weighted by real indegree"),
    weighting = factor(weighting, levels = rev(c("unweighted","weighted by sample indegree","weighted by real indegree"))),
    sample_method = factor(sample_method, levels = unique(sample_method)),
    statistic = case_match(statistic, "density" ~"Density",
                           "assortativity" ~ "Racial Assortativity",
                           "lcl_coef" ~ "Local Clustering Coef.",
                           "global_clc" ~ "Global Clustering Coef.", .default = statistic),
    statistic = factor(statistic, levels = c("Density","Local Clustering Coef.","Racial Assortativity","Global Clustering Coef.")),
    donated = if_else(sample_method != "Random Walk",donated,NA)) 


plot_k1_reweighting <- ggplot(statistics_reweighted_long, aes(y = interaction(weighting)))+ 
  geom_violin(aes(x = value, fill = sample_method))+
  #scale_x_discrete(limits = rev) +
  geom_vline(aes(xintercept = full_network, color = "hline"), linewidth = 1.2)+
  geom_point(aes(x = donated, color = "point"),fill = color_donated, size = 4.5, shape = 18,
             inherit.aes = T) +
  scale_color_manual(
    name = "",
    values = c("hline" = color_line, "point" = color_donated),
    labels = c("hline" = "Original Network", "point" = "BFS, Donated Seeds")
  ) +
  scale_fill_manual("", values = c("#CC79A7","#009E73","#0072B2"))+
  xlab("")+
  ylab("Reweighting")+
  facet_grid2(sample_method~statistic, scales = "free_x",
              switch = "x",
              strip = strip_themed(background_x = elem_list_rect(fill = c("white","white","white"), color = c("white","white","white")),
                                   background_y = elem_list_rect(fill = c("white","white","white"), color = c("#CC79A7","#009E73","#0072B2"))))+
  guides(fill  = guide_legend(order = 1),
         color = guide_legend(order = 2))+
  theme_bw()+
  theme(axis.text=element_text(size=15), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20), legend.position = "bottom", legend.box = "vertical",
        legend.text = element_text(size=20),strip.placement = "outside", strip.text = element_text(size = 15), legend.title = element_text(size=24)) +
  facetted_pos_scales(
    x = list(
      statistic == "Racial Assortativity" ~ scale_x_continuous(
        breaks = c(0.1, 0.2, 0.3)
      )
    )
  )

setEPS()
postscript(
  "k1_compare_statistics_reweighting.eps",
  width = 40 / 2.54,
  height = 22.5 / 2.54,
  #units = "cm",
  paper = "special"
)

plot_k1_reweighting

dev.off()

## Step 19: Tables for Appendix: Descriptives

options(#arrow.skip_nul = TRUE,
  scipen = 999,
  digits = 4)

input <- "REPOSITORY_PATH"

bias_types <- c("biasdegree_linear","partybias_democrat5percentmore","partybias_democrat10percentmore","partybias_democrat100percentmore",
                "partybias_republican5percentmore","partybias_republican10percentmore","partybias_republican100percentmore",
                "genderbias_women5percentmore","genderbias_women10percentmore","genderbias_women100percentmore",
                "genderbias_men5percentmore","genderbias_men10percentmore","genderbias_men100percentmore",
                "ethnicity_caucasian5percentmore","ethnicity_caucasian10percentmore","ethnicity_caucasian100percentmore",
                "ethnicity_asian5percentmore","ethnicity_asian10percentmore","ethnicity_asian100percentmore") #length: 19

# k1
#real donation
donation_k1 <- read.csv("S3N564k1network_sample_descriptives_BFS_donated.csv") %>% 
  mutate(seed_selection = "real donation",
         size = "k1") %>% 
  select(seed_selection, sample_number, num_nodes, num_edges, num_components, largest_component_size)

# simulated donations
samples_k1 <- read.csv(paste0(input, "S100N564k1network_sample_descriptives_BFS.csv")) %>% 
  mutate(seed_selection = "no bias")

for (bias in 1:length(bias_types)) {
  
  bias_type <- bias_types[bias]
  
  sample <- read.csv(paste0(input, "S100N564k1network_sample_descriptives_BFS_",bias_type,".csv")) %>% 
    mutate(seed_selection = bias_type)
  
  samples_k1 <- rbind(samples_k1,sample)
  
  rm(sample)
}
samples_k1 <- samples_k1 %>% 
  select(seed_selection, sample_number, num_nodes, num_edges, num_components, largest_component_size)

# random walk
rwj015_k1 <- read.csv(paste0(input, "S100N10373k1j0.15network_sample_descriptives_RW.csv")) %>% 
  mutate(seed_selection = "RW j = 0.15",
         bias_category = "RW",
         bias_explicit = "RW j = 0.15")

rw_k1 <- rwj015_k1 %>% 
  rename(sample_number = sample_num) %>% 
  select(seed_selection,sample_number, num_nodes, num_edges, num_components, largest_component_size)

samples_k1_long <- rbind(donation_k1,samples_k1,rw_k1) %>% 
  mutate(seed_selection= if_else(grepl("_", seed_selection), sub("^[^_]*_", "", seed_selection),seed_selection)) %>%
  mutate(largest_component_share = largest_component_size / num_nodes) %>% 
  group_by(seed_selection) %>% 
  summarise(num_nodes_av = mean(num_nodes),
            num_nodes_se = sd(num_nodes) / sqrt(n()),
            num_edges_av = mean(num_edges),
            num_edges_se = sd(num_edges) / sqrt(n()),
            num_components_av = mean(num_components),
            num_components_se = sd(num_components) / sqrt(n()),
            largest_component_share_av = as.character(format(round(mean(largest_component_share),3),nsmall = 3)),
            largest_component_share_se = as.character(format(round(sd(largest_component_share) / sqrt(n()),3),nsmall = 3)), #
            largest_component_share_se = if_else(largest_component_share_se == "NA","",
                                                 if_else(largest_component_share_se == "0.000", "(>.001)", paste0("(",str_sub(largest_component_share_se, 2),")")))
  ) %>% 
  mutate(num_nodes = if_else(grepl("RW|dona", seed_selection),paste0(format(round(num_nodes_av,0),nsmall = 0)),
                             paste0(format(round(num_nodes_av,0),nsmall = 0)," (",format(round(num_nodes_se,0),nsmall = 0),")")),
         num_edges = if_else(grepl("dona", seed_selection), paste0(format(round(num_edges_av,0),nsmall = 0)),paste0(format(round(num_edges_av,0),nsmall = 0)," (",format(round(num_edges_se,0),nsmall = 0),")")),
         num_components = paste0(format(round(num_components_av,0),nsmall = 0)," (",format(round(num_components_se,0),nsmall = 0),")"),
         share_in_largest_component = paste0(largest_component_share_av," ",largest_component_share_se)
  ) %>% 
  select(seed_selection,num_nodes,num_edges,
         share_in_largest_component) %>% 
  arrange(match(seed_selection, c("no bias","real donation","linear",
                                  "women5percentmore","women10percentmore","women100percentmore",
                                  "men5percentmore","men10percentmore","men100percentmore",
                                  "caucasian5percentmore","caucasian10percentmore","caucasian100percentmore",
                                  "asian5percentmore","asian10percentmore", "asian100percentmore",
                                  "democrat5percentmore","democrat10percentmore","democrat100percentmore",
                                  "republican5percentmore","republican10percentmore","republican100percentmore",
                                  "RW j = 0.15"
  ))) %>% 
  mutate(seed_selection = case_match(seed_selection,
                                     "no bias" ~ "no bias",
                                     "linear" ~ "indegree bias linear",
                                     "real donation" ~"seeds from data donation",
                                     "asian100percentmore" ~ "asians 100 percent more",
                                     "asian10percentmore" ~ "asians 10 percent more",
                                     "asian5percentmore" ~ "asians 5 percent more",
                                     "caucasian100percentmore" ~ "caucasians 100 percent more",
                                     "caucasian10percentmore" ~ "caucasians 10 percent more",
                                     "caucasian5percentmore" ~ "caucasians 5 percent more",
                                     "women100percentmore" ~ "women 100 percent more",
                                     "women10percentmore" ~ "women 10 percent more",
                                     "women5percentmore" ~ "women 5 percent more",
                                     "men100percentmore" ~ "men 100 percent more",
                                     "men10percentmore" ~ "men 10 percent more",
                                     "men5percentmore" ~ "men 5 percent more",
                                     "democrat100percentmore" ~ "democrats 100 percent more",
                                     "democrat10percentmore" ~ "democrats 10 percent more",
                                     "democrat5percentmore" ~ "democrats 5 percent more",
                                     "republican100percentmore" ~ "republicans 100 percent more",
                                     "republican10percentmore" ~ "republicans 10 percent more",
                                     "republican5percentmore" ~ "republicans 5 percent more",
                                     "RW j = 0.15" ~ "RW, j = 0.15",
                                     .default = NA)) %>% 
  rename_with(~ str_replace_all(.x, "_", " ")) 




int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate than democrats",
                   "Random Walk as silver standard, N = 10,373")


kbl <- samples_k1_long %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:descriptives1layer}Sample descriptives over 100 samples for sampling 1 layer via BFS and for samples collected via RW. Means with Standard Errors in braces.")

pos <- seq(3, nrow(samples_k1_long), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(samples_k1_long), "}{l}{\\textit{", int_headlines[i], "}} \\\\ \n"
      )
    )
}



kbl %>%
  kable_styling(latex_options = "HOLD_position")

#K2


# k2
#real donation
donation_k2 <- read.csv("S3N564k2network_sample_descriptives_BFS_donated.csv") %>% 
  mutate(seed_selection = "real donation",
         size = "k2") %>% 
  select(seed_selection, sample_number, num_nodes, num_edges, num_components, largest_component_size)

# simulated donations
samples_k2 <- read.csv(paste0(input, "S100N564k2network_sample_descriptives_BFS.csv")) %>% 
  mutate(seed_selection = "no bias")

for (bias in 1:length(bias_types)) {
  
  bias_type <- bias_types[bias]
  
  sample <- read.csv(paste0(input, "S100N564k2network_sample_descriptives_BFS_",bias_type,".csv")) %>% 
    mutate(seed_selection = bias_type)
  
  samples_k2 <- rbind(samples_k2,sample)
  
  rm(sample)
}
samples_k2 <- samples_k2 %>% 
  select(seed_selection, sample_number, num_nodes, num_edges, num_components, largest_component_size)

## random walk
rwj015_k2 <- read.csv(paste0(input, "S100N295426k2j0.15network_sample_descriptives_RW.csv")) %>% 
  mutate(seed_selection = "RW j = 0.15",
         bias_category = "RW",
         bias_explicit = "RW j = 0.15")

rw_k2 <-rwj015_k2 %>% 
  rename(sample_number = sample_num) %>% 
  select(seed_selection,sample_number, num_nodes, num_edges, num_components, largest_component_size)

samples_k2_long <- rbind(donation_k2,samples_k2,rw_k2) %>% 
  mutate(seed_selection= if_else(grepl("_", seed_selection), sub("^[^_]*_", "", seed_selection),seed_selection)) %>%
  mutate(largest_component_share = largest_component_size / num_nodes) %>% 
  group_by(seed_selection) %>% 
  summarise(num_nodes_av = mean(num_nodes),
            num_nodes_se = sd(num_nodes) / sqrt(n()),
            num_edges_av = mean(num_edges),
            num_edges_se = sd(num_edges) / sqrt(n()),
            num_components_av = mean(num_components),
            num_components_se = sd(num_components) / sqrt(n()),
            largest_component_share_av = as.character(format(round(mean(largest_component_share),3),nsmall = 3)),
            largest_component_share_se = as.character(format(round(sd(largest_component_share) / sqrt(n()),3),nsmall = 3)),
            largest_component_share_se = if_else(largest_component_share_se == "NA","",
                                                 if_else(largest_component_share_se == "0.000", "(>.001)", paste0("(",str_sub(largest_component_share_se, 2),")")))
  ) %>% 
  mutate(num_nodes = if_else(grepl("RW|dona", seed_selection),paste0(format(round(num_nodes_av,0),nsmall = 0)),
                             paste0(format(round(num_nodes_av,0),nsmall = 0)," (",format(round(num_nodes_se,0),nsmall = 0),")")),
         num_edges = if_else(grepl("dona", seed_selection), paste0(format(round(num_edges_av,0),nsmall = 0)),paste0(format(round(num_edges_av,0),nsmall = 0)," (",format(round(num_edges_se,0),nsmall = 0),")")),
         num_components = paste0(format(round(num_components_av,0),nsmall = 0)," (",format(round(num_components_se,0),nsmall = 0),")"),
         share_in_largest_component = paste0(largest_component_share_av," ",largest_component_share_se)
  ) %>%
  select(seed_selection,num_nodes,num_edges,share_in_largest_component) %>% 
  arrange(match(seed_selection, c("no bias","real donation","linear",
                                  "women5percentmore","women10percentmore","women100percentmore",
                                  "men5percentmore","men10percentmore","men100percentmore",
                                  "caucasian5percentmore","caucasian10percentmore","caucasian100percentmore",
                                  "asian5percentmore","asian10percentmore", "asian100percentmore",
                                  "democrat5percentmore","democrat10percentmore","democrat100percentmore",
                                  "republican5percentmore","republican10percentmore","republican100percentmore",
                                  "RW j = 0.15"
  ))) %>% 
  mutate(seed_selection = case_match(seed_selection,
                                     "no bias" ~ "no bias",
                                     "linear" ~ "indegree bias linear",
                                     "real donation" ~"seeds from data donation",
                                     "asian100percentmore" ~ "asians 100 percent more",
                                     "asian10percentmore" ~ "asians 10 percent more",
                                     "asian5percentmore" ~ "asians 5 percent more",
                                     "caucasian100percentmore" ~ "caucasians 100 percent more",
                                     "caucasian10percentmore" ~ "caucasians 10 percent more",
                                     "caucasian5percentmore" ~ "caucasians 5 percent more",
                                     "women100percentmore" ~ "women 100 percent more",
                                     "women10percentmore" ~ "women 10 percent more",
                                     "women5percentmore" ~ "women 5 percent more",
                                     "men100percentmore" ~ "men 100 percent more",
                                     "men10percentmore" ~ "men 10 percent more",
                                     "men5percentmore" ~ "men 5 percent more",
                                     "democrat100percentmore" ~ "democrats 100 percent more",
                                     "democrat10percentmore" ~ "democrats 10 percent more",
                                     "democrat5percentmore" ~ "democrats 5 percent more",
                                     "republican100percentmore" ~ "republicans 100 percent more",
                                     "republican10percentmore" ~ "republicans 10 percent more",
                                     "republican5percentmore" ~ "republicans 5 percent more",
                                     "RW j = 0.15" ~ "RW, j = 0.15",
                                     .default = NA)) %>% 
  rename_with(~ str_replace_all(.x, "_", " "))


int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate than democrats",
                   "Random Walk as silver standard, N = 295,426")


kbl <- samples_k2_long %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:descriptives2layers}Sample descriptives over 100 samples for sampling 2 layers via BFS and for samples collected via RW. Means with Standard Errors in braces")

pos <- seq(3, nrow(samples_k2_long), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(samples_k2_long), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}



kbl %>%
  kable_styling(latex_options = "HOLD_position")

## Step 20: Tables for Appendix: Attributes 

format_mean_se <- function(mean, se, digits = 3) {
  ifelse(
    is.na(se),
    sprintf(paste0("%.", digits, "f"), mean),
    sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, se)
  )
}

########## k1 #############

sample_summary <- samples_k1_long %>% 
  select(sample_number, attribute, share, seed_selection, share_full) %>% 
  group_by(seed_selection, attribute) %>%
  summarise(
    mean = mean(share, na.rm = TRUE),
    se = sd(share, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

full_summary <- samples_k1_long %>%
  distinct(attribute, share_full) %>%
  transmute(
    seed_selection = "full network",
    attribute,
    mean = share_full,
    se = NA_real_
  )

donation_summary <- donation_k1 %>%
  pivot_longer(
    cols = c(
      female, male,
      caucasian, african_american, asian, hispanic,
      democrat, republican, independent
    ),
    names_to = "attribute",
    values_to = "count"
  ) %>%
  transmute(
    seed_selection = "seeds from data donation",
    attribute,
    mean = count / sample_size,
    se = NA_real_
  )


attribute_summary <- bind_rows(
  full_summary,
  donation_summary,
  sample_summary
) %>%
  mutate(
    seed_selection = case_match(
      seed_selection,
      "no bias" ~ "no bias",
      "biasdegree_linear" ~ "indegree bias linear",
      "genderbias_women5percentmore" ~ "women 5 percent more",
      "genderbias_women10percentmore" ~ "women 10 percent more",
      "genderbias_women100percentmore" ~ "women 100 percent more",
      "genderbias_men5percentmore" ~ "men 5 percent more",
      "genderbias_men10percentmore" ~ "men 10 percent more",
      "genderbias_men100percentmore" ~ "men 100 percent more",
      "ethnicity_caucasian5percentmore" ~ "caucasians 5 percent more",
      "ethnicity_caucasian10percentmore" ~ "caucasians 10 percent more",
      "ethnicity_caucasian100percentmore" ~ "caucasians 100 percent more",
      "ethnicity_asian5percentmore" ~ "asians 5 percent more",
      "ethnicity_asian10percentmore" ~ "asians 10 percent more",
      "ethnicity_asian100percentmore" ~ "asians 100 percent more",
      "partybias_democrat5percentmore" ~ "democrats 5 percent more",
      "partybias_democrat10percentmore" ~ "democrats 10 percent more",
      "partybias_democrat100percentmore" ~ "democrats 100 percent more",
      "partybias_republican5percentmore" ~ "republicans 5 percent more",
      "partybias_republican10percentmore" ~ "republicans 10 percent more",
      "partybias_republican100percentmore" ~ "republicans 100 percent more",
      "RW j = 0.15" ~ "RW, j = 0.15",
      .default = seed_selection
    ),
    seed_selection = factor(
      seed_selection,
      levels = c(
        "full network",
        "no bias",
        "seeds from data donation",
        "indegree bias linear",
        "women 5 percent more",
        "women 10 percent more",
        "women 100 percent more",
        "men 5 percent more",
        "men 10 percent more",
        "men 100 percent more",
        "caucasians 5 percent more",
        "caucasians 10 percent more",
        "caucasians 100 percent more",
        "asians 5 percent more",
        "asians 10 percent more",
        "asians 100 percent more",
        "democrats 5 percent more",
        "democrats 10 percent more",
        "democrats 100 percent more",
        "republicans 5 percent more",
        "republicans 10 percent more",
        "republicans 100 percent more",
        "RW, j = 0.15"
      )
    )
  ) %>% 
  arrange(seed_selection)


attribute_table <- attribute_summary %>%
  mutate(value = format_mean_se(mean, se, 3)) %>%
  select(seed_selection, attribute, value) %>%
  pivot_wider(
    names_from = attribute,
    values_from = value
  )


####
int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate  than democrats",
                   "Random Walk as silver standard, N = 10,373")

# part 1: gender and party 
kbl <- attribute_table %>%
  select(`Seed selection` = seed_selection, Female = female, Male = male, Democrat = democrat, Republican = republican) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_biases_attributes_part1}Means and standard errors for attribute shares (in percent) for gender and party affiliation for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 4, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


# part 2: ethnicity
kbl <- attribute_table %>%
  select(`Seed selection` = seed_selection, Caucasian = caucasian, Asian = asian, Hispanic = hispanic, `African American` = african_american) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_biases_attributes_part2}Means and standard errors for attribute shares (in percent) for ethnicity for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 4, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


########## k2 #############

sample_summary <- samples_k2_long %>% 
  select(sample_number, attribute, share, seed_selection, share_full) %>% 
  group_by(seed_selection, attribute) %>%
  summarise(
    mean = mean(share, na.rm = TRUE),
    se = sd(share, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

full_summary <- samples_k2_long %>%
  distinct(attribute, share_full) %>%
  transmute(
    seed_selection = "full network",
    attribute,
    mean = share_full,
    se = NA_real_
  )

donation_summary <- donation_k2 %>%
  pivot_longer(
    cols = c(
      female, male,
      caucasian, african_american, asian, hispanic,
      democrat, republican, independent
    ),
    names_to = "attribute",
    values_to = "count"
  ) %>%
  transmute(
    seed_selection = "seeds from data donation",
    attribute,
    mean = count / sample_size,
    se = NA_real_
  )


attribute_summary <- bind_rows(
  full_summary,
  donation_summary,
  sample_summary
) %>%
  mutate(
    seed_selection = case_match(
      seed_selection,
      "no bias" ~ "no bias",
      "biasdegree_linear" ~ "indegree bias linear",
      "genderbias_women5percentmore" ~ "women 5 percent more",
      "genderbias_women10percentmore" ~ "women 10 percent more",
      "genderbias_women100percentmore" ~ "women 100 percent more",
      "genderbias_men5percentmore" ~ "men 5 percent more",
      "genderbias_men10percentmore" ~ "men 10 percent more",
      "genderbias_men100percentmore" ~ "men 100 percent more",
      "ethnicity_caucasian5percentmore" ~ "caucasians 5 percent more",
      "ethnicity_caucasian10percentmore" ~ "caucasians 10 percent more",
      "ethnicity_caucasian100percentmore" ~ "caucasians 100 percent more",
      "ethnicity_asian5percentmore" ~ "asians 5 percent more",
      "ethnicity_asian10percentmore" ~ "asians 10 percent more",
      "ethnicity_asian100percentmore" ~ "asians 100 percent more",
      "partybias_democrat5percentmore" ~ "democrats 5 percent more",
      "partybias_democrat10percentmore" ~ "democrats 10 percent more",
      "partybias_democrat100percentmore" ~ "democrats 100 percent more",
      "partybias_republican5percentmore" ~ "republicans 5 percent more",
      "partybias_republican10percentmore" ~ "republicans 10 percent more",
      "partybias_republican100percentmore" ~ "republicans 100 percent more",
      "RW j = 0.15" ~ "RW, j = 0.15",
      .default = seed_selection
    ),
    seed_selection = factor(
      seed_selection,
      levels = c(
        "full network",
        "no bias",
        "seeds from data donation",
        "indegree bias linear",
        "women 5 percent more",
        "women 10 percent more",
        "women 100 percent more",
        "men 5 percent more",
        "men 10 percent more",
        "men 100 percent more",
        "caucasians 5 percent more",
        "caucasians 10 percent more",
        "caucasians 100 percent more",
        "asians 5 percent more",
        "asians 10 percent more",
        "asians 100 percent more",
        "democrats 5 percent more",
        "democrats 10 percent more",
        "democrats 100 percent more",
        "republicans 5 percent more",
        "republicans 10 percent more",
        "republicans 100 percent more",
        "RW, j = 0.15"
      )
    )
  ) %>% 
  arrange(seed_selection)


attribute_table <- attribute_summary %>%
  mutate(value = format_mean_se(mean, se, 3)) %>%
  select(seed_selection, attribute, value) %>%
  pivot_wider(
    names_from = attribute,
    values_from = value
  )

# attribute_table %>%
#   kbl(booktabs = TRUE,
#       caption = "Mean share (SE) of attributes by seed selection.") %>%
#   kable_styling(full_width = FALSE)

####
int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate  than democrats",
                   "Random Walk as silver standard, N = 10,373")

# part 1: gender and party 
kbl <- attribute_table %>%
  select(`Seed selection` = seed_selection, Female = female, Male = male, Democrat = democrat, Republican = republican) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_biases_attributes_part1}Means and standard errors for attribute shares (in percent) for gender and party affiliation for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 4, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


# part 2: ethnicity
kbl <- attribute_table %>%
  select(`Seed selection` = seed_selection, Caucasian = caucasian, Asian = asian, Hispanic = hispanic, `African American` = african_american) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_biases_attributes_part2}Means and standard errors for attribute shares (in percent) for ethnicity for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 4, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")



### Reweighted Attributes 

## k1

attributes_reweighted_k1 <- read.csv2("attributes_weights.csv", sep = ";") %>% 
  filter(sample_depth == "k1")


attributes_reweighted_long_k1 <- attributes_reweighted_k1 %>%
  select(-sample_size, -starts_with(c("X","other","unknown","nogender","independent"))) %>% 
  filter(sample_depth == "k1") %>% 
  pivot_longer(
    cols = -c(sample_num,sample_method,sample_depth),
    names_to = c("attribute", "weighting"),
    names_pattern = "([^_]+)_(.*)"
  ) %>%
  select(-sample_depth) %>% 
  rbind(., full_202109_attributes_long %>%  mutate(sample_num = 1,
                                                   sample_method = "whole network",
                                                   attribute = case_match(attribute,
                                                                          "african_american" ~ "afram", .default = attribute),
                                                   weighting = "whole network",
                                                   value = full_network / 911331) %>% 
          filter(attribute != "independent") %>% 
          select(sample_num,sample_method,attribute,weighting,value))

attributes_reweighted_long_k1 <- attributes_reweighted_long_k1 %>%
  mutate(
    weighting = case_match(weighting,
                           "naive_mean" ~ "unweighted",
                           "hajek_sampindegree" ~ "weighted by sample indegree",
                           "hajek_realindegree" ~ "weighted by real indegree"),
    weighting = factor(weighting, levels = rev(c("whole network","unweighted","weighted by sample indegree","weighted by real indegree"))),
    sample_method = case_match(sample_method,
                               "BFS_biasdegree_linear" ~ "BFS popularity bias",
                               "RW0.15" ~ "RW j = 0.15",
                               .default = sample_method),
    sample_method = factor(sample_method, levels = c("whole network","BFS","donated","BFS popularity bias","RW j = 0.15"))) %>% 

  select(sample_number = sample_num, weighting, attribute, value, sample_method)


summary_table <- attributes_reweighted_long_k1 %>%
  pivot_wider(names_from = attribute,
              values_from = value) %>% 
  group_by(sample_method, weighting) %>%
  summarise(
    across(
      c(female,male,democrat,republican,caucasian,asian,hispanic,afram),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  transmute(
    `Sampling method` = sample_method,
    Correction = weighting,
    Female = format_mean_se(female_mean, female_se, 3),
    Male = format_mean_se(male_mean, male_se, 3),
    Democrat = format_mean_se(democrat_mean, democrat_se, 3),
    Republican = format_mean_se(republican_mean, republican_se, 3),
    Caucasian = format_mean_se(caucasian_mean, caucasian_se, 3),
    Asian = format_mean_se(asian_mean, asian_se, 3),
    Hispanic = format_mean_se(hispanic_mean, hispanic_se, 3),
    `African American` = format_mean_se(afram_mean, afram_se, 3)
  ) %>% 
  mutate(`Sampling method` = factor(`Sampling method`, levels = c("whole network","BFS","donated","BFS popularity bias","RW j = 0.15")),
         Correction = factor(Correction, levels = c("whole network","unweighted","weighted by sample indegree","weighted by real indegree"))) %>% 
  arrange(`Sampling method`,Correction) %>% 
  select(-`Sampling method`)

int_headlines <- c("BFS, no bias",
                   "BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 10,373")

#part 1 : Female, Male, Democrat, Republican
kbl <- summary_table %>%
  select(Correction, Female, Male, Democrat, Republican) %>%  
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_reweighted_attributes_tbl_part1}Means and standard errors for reweighted for attribute shares (in percent) for gender and party affiliation. 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(1, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")



#part 2 : Ethnicity
kbl <- summary_table %>%
  select(Correction, Caucasian,Asian,Hispanic,`African American`) %>%  
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_reweighted_attributes_tbl_part2}Means and standard errors for reweighted for attribute shares (in percent) for ethnicity. 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(1, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")

## k2

attributes_reweighted_k2 <- read.csv2("attributes_weights.csv", sep = ";") %>% 
  filter(sample_depth == "k2")


attributes_reweighted_long_k2 <- attributes_reweighted_k2 %>%
  select(-sample_size, -starts_with(c("X","other","unknown","nogender","independent"))) %>% 
  filter(sample_depth == "k2") %>% 
  pivot_longer(
    cols = -c(sample_num,sample_method,sample_depth),
    names_to = c("attribute", "weighting"),
    names_pattern = "([^_]+)_(.*)"
  ) %>%
  select(-sample_depth) %>% 
  rbind(., full_202109_attributes_long %>%  mutate(sample_num = 1,
                                                   sample_method = "whole network",
                                                   attribute = case_match(attribute,
                                                                          "african_american" ~ "afram", .default = attribute),
                                                   weighting = "whole network",
                                                   value = full_network / 911331) %>% 
          filter(attribute != "independent") %>% 
          select(sample_num,sample_method,attribute,weighting,value))

attributes_reweighted_long_k2 <- attributes_reweighted_long_k2 %>%
  mutate(weighting = case_match(weighting,
                                "naive_mean" ~ "unweighted",
                                "hajek_sampindegree" ~ "weighted by sample indegree",
                                "hajek_realindegree" ~ "weighted by real indegree"),
         weighting = factor(weighting, levels = rev(c("whole network","unweighted","weighted by sample indegree","weighted by real indegree"))),
         sample_method = case_match(sample_method,
                                    "BFS_biasdegree_linear" ~ "BFS popularity bias",
                                    "RW0.15" ~ "RW j = 0.15",
                                    .default = sample_method),
         sample_method = factor(sample_method, levels = c("whole network","BFS","donated","BFS popularity bias","RW j = 0.15"))) %>% 
  select(sample_number = sample_num, weighting, attribute, value, sample_method) # %>% 


summary_table <- attributes_reweighted_long_k2 %>%
  pivot_wider(names_from = attribute,
              values_from = value) %>% 
  group_by(sample_method, weighting) %>%
  summarise(
    across(
      c(female,male,democrat,republican,caucasian,asian,hispanic,afram),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  transmute(
    `Sampling method` = sample_method,
    Correction = weighting,
    Female = format_mean_se(female_mean, female_se, 3),
    Male = format_mean_se(male_mean, male_se, 3),
    Democrat = format_mean_se(democrat_mean, democrat_se, 3),
    Republican = format_mean_se(republican_mean, republican_se, 3),
    Caucasian = format_mean_se(caucasian_mean, caucasian_se, 3),
    Asian = format_mean_se(asian_mean, asian_se, 3),
    Hispanic = format_mean_se(hispanic_mean, hispanic_se, 3),
    `African American` = format_mean_se(afram_mean, afram_se, 3)
  ) %>% 
  mutate(`Sampling method` = factor(`Sampling method`, levels = c("whole network","BFS","donated","BFS popularity bias","RW j = 0.15")),
         Correction = factor(Correction, levels = c("whole network","unweighted","weighted by sample indegree","weighted by real indegree"))) %>% 
  arrange(`Sampling method`,Correction) %>% 
  select(-`Sampling method`)

int_headlines <- c("BFS, no bias",
                   "BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 295,426")

#part 1 : Female, Male, Democrat, Republican
kbl <- summary_table %>%
  select(Correction, Female, Male, Democrat, Republican) %>%  
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_reweighted_attributes_tbl_part1}Means and standard errors for reweighted for attribute shares (in percent) for gender and party affiliation. 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(1, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")



#part 2 : Ethnicity
kbl <- summary_table %>%
  select(Correction, Caucasian,Asian,Hispanic,`African American`) %>%  
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_reweighted_attributes_tbl_part2}Means and standard errors for reweighted for attribute shares (in percent) for ethnicity. 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(1, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


## Step 21: Tables for Appendix: Network Statistics

########## k1 #############

samples_k1 <- read.csv("statistics_samples_k1.csv")

all_k1 <- samples_k1 %>% 
  select(sample_number, avg_indegree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic, seed_selection) %>% 
  rbind(., full_202109_stats %>% select(sample_number, avg_indegree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic) %>% mutate(seed_selection = "full network")) %>% 
  rbind(., donation_k1 %>% select(sample_number, avg_indegree = avg_degree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic)%>% mutate(seed_selection = "seeds from data donation")) %>%
  rbind(., rwj015_k1 %>% select(sample_number = sample_num, avg_indegree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic,seed_selection)) %>%
  mutate(
    seed_selection = case_match(
      seed_selection,
      "full network" ~ "full network",
      "no bias" ~ "no bias",
      "seeds from data donation" ~ "seeds from data donation",
      "biasdegree_linear" ~ "indegree bias linear",
      "genderbias_women5percentmore" ~ "women 5 percent more",
      "genderbias_women10percentmore" ~ "women 10 percent more",
      "genderbias_women100percentmore" ~ "women 100 percent more",
      "genderbias_men5percentmore" ~ "men 5 percent more",
      "genderbias_men10percentmore" ~ "men 10 percent more",
      "genderbias_men100percentmore" ~ "men 100 percent more",
      "ethnicity_caucasian5percentmore" ~ "caucasians 5 percent more",
      "ethnicity_caucasian10percentmore" ~ "caucasians 10 percent more",
      "ethnicity_caucasian100percentmore" ~ "caucasians 100 percent more",
      "ethnicity_asian5percentmore" ~ "asians 5 percent more",
      "ethnicity_asian10percentmore" ~ "asians 10 percent more",
      "ethnicity_asian100percentmore" ~ "asians 100 percent more",
      "partybias_democrat5percentmore" ~ "democrats 5 percent more",
      "partybias_democrat10percentmore" ~ "democrats 10 percent more",
      "partybias_democrat100percentmore" ~ "democrats 100 percent more",
      "partybias_republican5percentmore" ~ "republicans 5 percent more",
      "partybias_republican10percentmore" ~ "republicans 10 percent more",
      "partybias_republican100percentmore" ~ "republicans 100 percent more",
      "RW j = 0.15" ~ "RW, j = 0.15",
      .default = seed_selection
    ),
    seed_selection = factor(
      seed_selection,
      levels = c(
        "full network",
        "no bias",
        "seeds from data donation",
        "indegree bias linear",
        "women 5 percent more",
        "women 10 percent more",
        "women 100 percent more",
        "men 5 percent more",
        "men 10 percent more",
        "men 100 percent more",
        "caucasians 5 percent more",
        "caucasians 10 percent more",
        "caucasians 100 percent more",
        "asians 5 percent more",
        "asians 10 percent more",
        "asians 100 percent more",
        "democrats 5 percent more",
        "democrats 10 percent more",
        "democrats 100 percent more",
        "republicans 5 percent more",
        "republicans 10 percent more",
        "republicans 100 percent more",
        "RW, j = 0.15"
      )
    )
  )

format_mean_se <- function(mean, se, digits = 3) {
  ifelse(
    is.na(se),
    sprintf(paste0("%.", digits, "f"), mean),
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f)"),
      mean, se
    )
  )
}


summary_table <- all_k1 %>%
  group_by(seed_selection) %>%
  summarise(
    across(
      c(avg_indegree,
        density,
        clustering_coeff,
        clustering_coeff_global_vector,
        assortativity_ethnic),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  transmute(
    `Seed selection` = seed_selection,
    `Avg. indegree` = format_mean_se(avg_indegree_mean, avg_indegree_se,3),
    Density = format_mean_se(density_mean, density_se,5),
    `Avg. local clustering` =
      format_mean_se(clustering_coeff_mean, clustering_coeff_se,3),
    `Global clustering` =
      format_mean_se(clustering_coeff_global_vector_mean,
                     clustering_coeff_global_vector_se,5),
    `Ethnic assortativity` =
      format_mean_se(assortativity_ethnic_mean,
                     assortativity_ethnic_se,3)
  )

####
int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate  than democrats",
                   "Random Walk as silver standard, N = 10,373")

# part 1: Density, `Avg. local clustering`, `Ethnic assortativity`
kbl <- summary_table %>%
  select(`Seed selection`, Density, `Avg. local clustering coef.` = `Avg. local clustering`, `Ethnic assortativity`) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_biases_statistics_part1}Means and standard errors for network statistics for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 4, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")

### part 2: `Global clustering`, `Avg. indegree`
kbl <- summary_table %>%
  select(`Seed selection`, `Avg. indegree`,`Global clustering coef.` = `Global clustering`) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_biases_statistics_part2}Means and standard errors for network statistics for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


########## k2 #############

samples_k2 <- read.csv("statistics_samples_k2.csv")

all_k2 <- samples_k2 %>% 
  select(sample_number, avg_indegree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic, seed_selection) %>% 
  rbind(., full_202109_stats %>% select(sample_number, avg_indegree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic) %>% mutate(seed_selection = "full network")) %>% 
  rbind(., donation_k2 %>% select(sample_number, avg_indegree = avg_degree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic)%>% mutate(seed_selection = "seeds from data donation")) %>%
  rbind(., rwj015_k2 %>% select(sample_number = sample_num, avg_indegree, density, clustering_coeff,clustering_coeff_global_vector,assortativity_ethnic,seed_selection)) %>%
  mutate(
    seed_selection = case_match(
      seed_selection,
      "full network" ~ "full network",
      "no bias" ~ "no bias",
      "seeds from data donation" ~ "seeds from data donation",
      "biasdegree_linear" ~ "indegree bias linear",
      "genderbias_women5percentmore" ~ "women 5 percent more",
      "genderbias_women10percentmore" ~ "women 10 percent more",
      "genderbias_women100percentmore" ~ "women 100 percent more",
      "genderbias_men5percentmore" ~ "men 5 percent more",
      "genderbias_men10percentmore" ~ "men 10 percent more",
      "genderbias_men100percentmore" ~ "men 100 percent more",
      "ethnicity_caucasian5percentmore" ~ "caucasians 5 percent more",
      "ethnicity_caucasian10percentmore" ~ "caucasians 10 percent more",
      "ethnicity_caucasian100percentmore" ~ "caucasians 100 percent more",
      "ethnicity_asian5percentmore" ~ "asians 5 percent more",
      "ethnicity_asian10percentmore" ~ "asians 10 percent more",
      "ethnicity_asian100percentmore" ~ "asians 100 percent more",
      "partybias_democrat5percentmore" ~ "democrats 5 percent more",
      "partybias_democrat10percentmore" ~ "democrats 10 percent more",
      "partybias_democrat100percentmore" ~ "democrats 100 percent more",
      "partybias_republican5percentmore" ~ "republicans 5 percent more",
      "partybias_republican10percentmore" ~ "republicans 10 percent more",
      "partybias_republican100percentmore" ~ "republicans 100 percent more",
      "RW j = 0.15" ~ "RW, j = 0.15",
      .default = seed_selection
    ),
    seed_selection = factor(
      seed_selection,
      levels = c(
        "full network",
        "no bias",
        "seeds from data donation",
        "indegree bias linear",
        "women 5 percent more",
        "women 10 percent more",
        "women 100 percent more",
        "men 5 percent more",
        "men 10 percent more",
        "men 100 percent more",
        "caucasians 5 percent more",
        "caucasians 10 percent more",
        "caucasians 100 percent more",
        "asians 5 percent more",
        "asians 10 percent more",
        "asians 100 percent more",
        "democrats 5 percent more",
        "democrats 10 percent more",
        "democrats 100 percent more",
        "republicans 5 percent more",
        "republicans 10 percent more",
        "republicans 100 percent more",
        "RW, j = 0.15"
      )
    )
  )


summary_table <- all_k2 %>%
  group_by(seed_selection) %>%
  summarise(
    across(
      c(avg_indegree,
        density,
        clustering_coeff,
        clustering_coeff_global_vector,
        assortativity_ethnic),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  transmute(
    `Seed selection` = seed_selection,
    `Avg. indegree` = format_mean_se(avg_indegree_mean, avg_indegree_se,3),
    Density = format_mean_se(density_mean, density_se,5),
    `Avg. local clustering` =
      format_mean_se(clustering_coeff_mean, clustering_coeff_se,3),
    `Global clustering` =
      format_mean_se(clustering_coeff_global_vector_mean,
                     clustering_coeff_global_vector_se,5),
    `Ethnic assortativity` =
      format_mean_se(assortativity_ethnic_mean,
                     assortativity_ethnic_se,3)
  )

####
int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate  than democrats",
                   "Random Walk as silver standard, N = 295,426")

# part 1: Density, `Avg. local clustering`, `Ethnic assortativity`
kbl <- summary_table %>%
  select(`Seed selection`, Density, `Avg. local clustering coef.` = `Avg. local clustering`, `Ethnic assortativity`) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_biases_statistics_part1}Means and standard errors for network statistics for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 4, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")

### part 2: `Global clustering`, `Avg. indegree`
kbl <- summary_table %>%
  select(`Seed selection`, `Avg. indegree`,`Global clustering coef.` = `Global clustering`) %>% 
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_biases_statistics_part2}Means and standard errors for network statistics for the full network, BFS with donated seeds, 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(4, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


### Reweighted Statistics 

## k1

statistics_reweighted_k1 <- read.csv("", sep = ",")


statistics_reweighted_long_k1 <- statistics_reweighted_k1 %>%
  select(-starts_with("X"), -sample_size.y) %>%
  filter(sample_depth == "k1") %>% 
  pivot_longer(
    cols = -c(sample_num,sample_method,sample_depth,sample_size.x),
    names_to = c("statistic", "weighting"),
    names_pattern = "^((?:global|lcl)_[^_]+|[^_]+)_(.*)$"
  ) %>%
  left_join(., full_202109_stats_long %>%  mutate(full_network = full_network,
                                                  statistic = case_match(statistic,
                                                                         "assortativity_ethnic" ~"assortativity",
                                                                         "density" ~ "density",
                                                                         "clustering_coeff" ~ "lcl_coef",
                                                                         "clustering_coeff_global_vector" ~ "global_clc",
                                                                         .default = statistic)), by = "statistic") 

statistics_reweighted_long_donated_k1 <- statistics_reweighted_long_k1 %>% 
  filter(sample_method == "donated") 

statistics_reweighted_long_k1 <- statistics_reweighted_long_k1 %>%
  filter(sample_method != "donated") %>%
  left_join(.,statistics_reweighted_long_donated_k1 %>% select(donated = value,statistic,weighting), by = c("statistic","weighting")) %>% 
  mutate(sample_method = case_match(sample_method,
                                    "RW0.15" ~ "RW j = 0.15",
                                    "BFS_biasdegree_linear" ~ "BFS popularity bias",
                                    "BFS" ~ "BFS"),
         weighting = case_match(weighting,
                                "naive" ~ "unweighted",
                                "naive_mean" ~ "unweighted",
                                "hajek_naive" ~ "unweighted",
                                "hajek_sampindegree" ~ "weighted by sample indegree",
                                "hajek_realindegree" ~ "weighted by real indegree"),
         weighting = factor(weighting, levels = rev(c("unweighted","weighted by sample indegree","weighted by real indegree"))),
         sample_method = factor(sample_method, levels = unique(sample_method))) %>% 
  select(sample_number = sample_num, weighting, statistic, value, sample_method) %>% 
  rbind(., statistics_reweighted_long_donated_k1 %>% select(sample_number = sample_num, weighting, statistic, value, sample_method) %>% 
          mutate(weighting = case_match(weighting,
                                        "naive" ~ "unweighted",
                                        "naive_mean" ~ "unweighted",
                                        "hajek_naive" ~ "unweighted",
                                        "hajek_sampindegree" ~ "weighted by sample indegree",
                                        "hajek_realindegree" ~ "weighted by real indegree")))


summary_table <- statistics_reweighted_long_k1 %>%
  pivot_wider(names_from = statistic,
              values_from = value) %>% 
  group_by(sample_method, weighting) %>%
  summarise(
    across(
      c(density,
        lcl_coef,
        global_clc,
        assortativity),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  transmute(
    `Sampling method` = sample_method,
    Correction = weighting,
    Density = format_mean_se(density_mean, density_se,6),
    `Avg. local clustering` =
      format_mean_se(lcl_coef_mean, lcl_coef_se,3),
    `Global clustering` =
      format_mean_se(global_clc_mean,
                     global_clc_se,5),
    `Ethnic assortativity` =
      format_mean_se(assortativity_mean,
                     assortativity_se,3)
  ) %>% 
  rbind(full_202109_stats %>% select(`Sampling method` = sample_number,
                                     Correction  = sample_number,
                                     Density = density,
                                     `Avg. local clustering` = clustering_coeff,
                                     `Global clustering` = clustering_coeff_global_vector,
                                     `Ethnic assortativity` = assortativity_ethnic) %>% 
          mutate(Density = round(Density,6),
                 `Avg. local clustering` = round(`Avg. local clustering`,3),
                 `Global clustering` = round(`Global clustering`,5),
                 `Ethnic assortativity` = round(`Ethnic assortativity`,3))) %>% 
  mutate(`Sampling method` = factor(`Sampling method`, levels = c("whole network","BFS","donated","BFS popularity bias","RW j = 0.15")),
         Correction = factor(Correction, levels = c("whole network","unweighted","weighted by sample indegree","weighted by real indegree"))) %>% 
  arrange(`Sampling method`,Correction) %>% 
  select(-`Sampling method`)

int_headlines <- c("BFS, no bias",
                   "BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 10,373")


kbl <- summary_table %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k1_reweighted_statistics_tbl}Means and standard errors for reweighted network statistics. 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 1 layer via BFS.") 

pos <- seq(1, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


## k2

statistics_reweighted_k2 <- read.csv("", sep = ",")


statistics_reweighted_long_k2 <- statistics_reweighted_k2 %>%
  select(-starts_with("X"), -sample_size.y) %>%
  filter(sample_depth == "k2") %>% 
  pivot_longer(
    cols = -c(sample_num,sample_method,sample_depth,sample_size.x),
    names_to = c("statistic", "weighting"),
    names_pattern = "^((?:global|lcl)_[^_]+|[^_]+)_(.*)$"
  ) %>%
  left_join(., full_202109_stats_long %>%  mutate(full_network = full_network,
                                                  statistic = case_match(statistic,
                                                                         "assortativity_ethnic" ~"assortativity",
                                                                         "density" ~ "density",
                                                                         "clustering_coeff" ~ "lcl_coef",
                                                                         "clustering_coeff_global_vector" ~ "global_clc",
                                                                         .default = statistic)), by = "statistic")

statistics_reweighted_long_donated_k2 <- statistics_reweighted_long_k2 %>% 
  filter(sample_method == "donated") 

statistics_reweighted_long_k2 <- statistics_reweighted_long_k2 %>%
  filter(sample_method != "donated") %>%
  left_join(.,statistics_reweighted_long_donated_k2 %>% select(donated = value,statistic,weighting), by = c("statistic","weighting")) %>% 
  mutate(sample_method = case_match(sample_method,
                                    "RW0.15" ~ "RW j = 0.15",
                                    "BFS_biasdegree_linear" ~ "BFS popularity bias",
                                    "BFS" ~ "BFS"),
         weighting = case_match(weighting,
                                "naive" ~ "unweighted",
                                "naive_mean" ~ "unweighted",
                                "hajek_naive" ~ "unweighted",
                                "hajek_sampindegree" ~ "weighted by sample indegree",
                                "hajek_realindegree" ~ "weighted by real indegree"),
         weighting = factor(weighting, levels = rev(c("unweighted","weighted by sample indegree","weighted by real indegree"))),
         sample_method = factor(sample_method, levels = unique(sample_method))) %>% 
  select(sample_number = sample_num, weighting, statistic, value, sample_method) %>% 
  rbind(., statistics_reweighted_long_donated_k2 %>% select(sample_number = sample_num, weighting, statistic, value, sample_method) %>% 
          mutate(weighting = case_match(weighting,
                                        "naive" ~ "unweighted",
                                        "naive_mean" ~ "unweighted",
                                        "hajek_naive" ~ "unweighted",
                                        "hajek_sampindegree" ~ "weighted by sample indegree",
                                        "hajek_realindegree" ~ "weighted by real indegree"))) %>% 
  filter(statistic != "global_clc")


summary_table <- statistics_reweighted_long_k2 %>%
  pivot_wider(names_from = statistic,
              values_from = value) %>% 
  group_by(sample_method, weighting) %>%
  summarise(
    across(
      c(density,
        lcl_coef,
        assortativity),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        se = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  transmute(
    `Sampling method` = sample_method,
    Correction = weighting,
    Density = format_mean_se(density_mean, density_se,7),
    `Avg. local clustering` =
      format_mean_se(lcl_coef_mean, lcl_coef_se,3),
    `Ethnic assortativity` =
      format_mean_se(assortativity_mean,
                     assortativity_se,3)
  ) %>% 
  rbind(full_202109_stats %>% select(`Sampling method` = sample_number,
                                     Correction  = sample_number,
                                     Density = density,
                                     `Avg. local clustering` = clustering_coeff,
                                     `Ethnic assortativity` = assortativity_ethnic) %>% 
          mutate(Density = round(Density,7),
                 `Avg. local clustering` = round(`Avg. local clustering`,3),
                 `Ethnic assortativity` = round(`Ethnic assortativity`,3))) %>% 
  mutate(`Sampling method` = factor(`Sampling method`, levels = c("whole network","BFS","donated","BFS popularity bias","RW j = 0.15")),
         Correction = factor(Correction, levels = c("whole network","unweighted","weighted by sample indegree","weighted by real indegree"))) %>% 
  arrange(`Sampling method`,Correction) %>% 
  select(-`Sampling method`)

int_headlines <- c("BFS, no bias",
                   "BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 295,426")


kbl <- summary_table %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:k2_reweighted_statistics_tbl}Means and standard errors for reweighted network statistics. 100 samples samples per seed selection bis in BFS and 100 samples for RW. Sampling 2 layers via BFS.") 

pos <- seq(1, nrow(summary_table), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", 3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")

## Step 21: Tables for Appendix: NMSE


NMSE_attributes_biases <- read.csv("NMSE_attributes_biases_new.csv", sep = ",")

NMSE_attributes_biases <- NMSE_attributes_biases %>% 
  separate(sample_bias, into = c("bias_attribute", "bias"), sep = "_", extra = "merge") %>% 
  mutate(across(where(is.numeric),~ round(.x,3))) %>% 
  mutate(across(where(is.numeric),~ formatC(.x, format = "f", digits = 3)),
         across(where(is.character), ~ ifelse(.x == "0.000", ">0.001", .x)),
         across(where(is.character), ~ ifelse(grepl("NA",.x), "", .x)),
         across(ends_with("_se"), ~ ifelse(.x == "",.x,paste0(" (",.x,")")))) %>% 
  mutate(#NMSE
    NMSE_gender = paste0(gender_NMSE_av,gender_NMSE_se),
    NMSE_ethnicity = paste0(ethnicity_NMSE_av,ethnicity_NMSE_se),
    NMSE_party = paste0(party_NMSE_av,party_NMSE_se),
    bias = if_else(is.na(bias),bias_attribute,bias)) %>% 
  select(bias, sample_depth, starts_with("NMSE")) %>% 
  pivot_wider(
    id_cols = bias,
    names_from = sample_depth,
    values_from = c(NMSE_gender, NMSE_ethnicity, NMSE_party),
    names_sep = "."
  ) %>% 
  select(bias, ends_with("k1"),ends_with("k2")) %>% 
  arrange(match(bias, c("no bias","donated","linear",
                        "women5percentmore","women10percentmore","women100percentmore",
                        "men5percentmore","men10percentmore","men100percentmore",
                        "caucasian5percentmore","caucasian10percentmore","caucasian100percentmore",
                        "asian5percentmore","asian10percentmore", "asian100percentmore",
                        "democrat5percentmore","democrat10percentmore","democrat100percentmore",
                        "republican5percentmore","republican10percentmore","republican100percentmore",
                        "RW0.15"
  ))) %>% 
  mutate(bias = case_match(bias,
                           "no bias" ~ "no bias",
                           "linear" ~ "indegree bias linear",
                           "donated" ~"seeds from data donation",
                           "asian100percentmore" ~ "asians 100 percent more",
                           "asian10percentmore" ~ "asians 10 percent more",
                           "asian5percentmore" ~ "asians 5 percent more",
                           "caucasian100percentmore" ~ "caucasians 100 percent more",
                           "caucasian10percentmore" ~ "caucasians 10 percent more",
                           "caucasian5percentmore" ~ "caucasians 5 percent more",
                           "women100percentmore" ~ "women 100 percent more",
                           "women10percentmore" ~ "women 10 percent more",
                           "women5percentmore" ~ "women 5 percent more",
                           "men100percentmore" ~ "men 100 percent more",
                           "men10percentmore" ~ "men 10 percent more",
                           "men5percentmore" ~ "men 5 percent more",
                           "democrat100percentmore" ~ "democrats 100 percent more",
                           "democrat10percentmore" ~ "democrats 10 percent more",
                           "democrat5percentmore" ~ "democrats 5 percent more",
                           "republican100percentmore" ~ "republicans 100 percent more",
                           "republican10percentmore" ~ "republicans 10 percent more",
                           "republican5percentmore" ~ "republicans 5 percent more",
                           "RW0.15" ~ "RW, j = 0.15",
                           .default = NA))


# k1
NMSE_attributes_biases_k1 <- NMSE_attributes_biases %>% 
  select(bias, ends_with("k1"))

int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate  than democrats",
                   "Random Walk as silver standard, N = 10,373")

kbl <- NMSE_attributes_biases_k1 %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:NMSE_k1_biases_attributes}Average normalized mean square error and standard errors for nodal attribute shares over 100 samples for sampling 1 layer via BFS") 

pos <- seq(3, nrow(NMSE_attributes_biases_k1), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_attributes_biases_k1), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")

# k2
NMSE_attributes_biases_k2 <- NMSE_attributes_biases %>% 
  select(bias, ends_with("k2"))

int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate  than democrats",
                   "Random Walk as silver standard, N = 295,426")


kbl <- NMSE_attributes_biases_k2 %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:NMSE_k2_biases_attributes}Average normalized mean square error and standard errors for nodal attribute shares over 100 samples for sampling 2 layers via BFS")

pos <- seq(3, nrow(NMSE_attributes_biases_k2), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_attributes_biases_k2), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")

# NMSE network statistics
NMSE_statistics_biases <- read.csv("NMSE_statistics_biases.csv", sep = ",")

NMSE_statistics_biases <- NMSE_statistics_biases %>% 
  separate(sample_method_bias, into = c("bias_attribute", "bias"), sep = "_", extra = "merge") %>% 
  mutate(across(where(is.numeric),~ round(.x,3))) %>% 
  mutate(across(where(is.numeric),~ formatC(.x, format = "f", digits = 3)),
         across(where(is.character), ~ ifelse(.x == "0.000", ">0.001", .x)),
         across(where(is.character), ~ ifelse(grepl("NA",.x), "", .x)),
         across(ends_with("_se"), ~ ifelse(.x == "",.x,paste0(" (",.x,")")))) %>% 
  mutate(#NMSE
    NMSE_degree = paste0(avg_degree_NMSE_av,avg_degree_NMSE_se),
    NMSE_local_clustering_coef = paste0(avg_local_clustering_coeff_NMSE_av,avg_local_clustering_coeff_NMSE_se),
    NMSE_ethnic_assortativity = paste0(assortativity_ethnic_NMSE_av,assortativity_ethnic_NMSE_se),
    NMSE_density = paste0(density_NMSE_av,density_NMSE_se),
    NMSE_global_clustering_coef = paste0(clustering_coeff_global_full_NMSE_av,clustering_coeff_global_full_NMSE_se),
    bias = if_else(is.na(bias),bias_attribute,bias)) %>% 
  select(bias, sample_depth, starts_with("NMSE")) %>% 
  pivot_wider(
    id_cols = bias,
    names_from = sample_depth,
    values_from = c(NMSE_degree, NMSE_local_clustering_coef, NMSE_ethnic_assortativity, NMSE_density,NMSE_global_clustering_coef),
    names_sep = "."
  ) %>% 
  select(bias, ends_with("k1"),ends_with("k2")) %>% 
  arrange(match(bias, c("no bias","donated","linear",
                        "women5percentmore","women10percentmore","women100percentmore",
                        "men5percentmore","men10percentmore","men100percentmore",
                        "caucasian5percentmore","caucasian10percentmore","caucasian100percentmore",
                        "asian5percentmore","asian10percentmore", "asian100percentmore",
                        "democrat5percentmore","democrat10percentmore","democrat100percentmore",
                        "republican5percentmore","republican10percentmore","republican100percentmore",
                        "RW0.15"
  ))) %>% 
  mutate(bias = case_match(bias,
                           "no bias" ~ "no bias",
                           "linear" ~ "indegree bias linear",
                           "donated" ~"seeds from data donation",
                           "asian100percentmore" ~ "asians 100 percent more",
                           "asian10percentmore" ~ "asians 10 percent more",
                           "asian5percentmore" ~ "asians 5 percent more",
                           "caucasian100percentmore" ~ "caucasians 100 percent more",
                           "caucasian10percentmore" ~ "caucasians 10 percent more",
                           "caucasian5percentmore" ~ "caucasians 5 percent more",
                           "women100percentmore" ~ "women 100 percent more",
                           "women10percentmore" ~ "women 10 percent more",
                           "women5percentmore" ~ "women 5 percent more",
                           "men100percentmore" ~ "men 100 percent more",
                           "men10percentmore" ~ "men 10 percent more",
                           "men5percentmore" ~ "men 5 percent more",
                           "democrat100percentmore" ~ "democrats 100 percent more",
                           "democrat10percentmore" ~ "democrats 10 percent more",
                           "democrat5percentmore" ~ "democrats 5 percent more",
                           "republican100percentmore" ~ "republicans 100 percent more",
                           "republican10percentmore" ~ "republicans 10 percent more",
                           "republican5percentmore" ~ "republicans 5 percent more",
                           "RW0.15" ~ "RW, j = 0.15",
                           .default = NA))



NMSE_statistics_biases_k1_den_clc_asort <- NMSE_statistics_biases %>% 
  select(bias, NMSE_density.k1, NMSE_local_clustering_coef.k1, NMSE_ethnic_assortativity.k1)

NMSE_statistics_biases_k1_deg_glc <- NMSE_statistics_biases %>% 
  select(bias, NMSE_degree.k1, NMSE_global_clustering_coef.k1)

NMSE_statistics_biases_k2_den_clc_asort <- NMSE_statistics_biases %>% 
  select(bias, NMSE_density.k2, NMSE_local_clustering_coef.k2, NMSE_ethnic_assortativity.k2)

NMSE_statistics_biases_k2_deg_glc <- NMSE_statistics_biases %>% 
  select(bias, NMSE_degree.k2, NMSE_global_clustering_coef.k2)

int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate than democrats",
                   "Random Walk as silver standard, N = 10,373")

kbl <- NMSE_statistics_biases_k1_den_clc_asort %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:NMSE_k1_biases_statistics_part1}Average normalized mean square error and standard errors for network statistics over 100 samples for sampling 1 layer via BFS. For the real data donation only one sample exists.") %>%
  kable_styling(latex_options = "hold_position") 

pos <- seq(3, nrow(NMSE_statistics_biases_k1_den_clc_asort), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_statistics_biases_k1_den_clc_asort), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")




kbl <- NMSE_statistics_biases_k1_deg_glc %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:NMSE_k1_biases_statistics_part2}Average normalized mean square error and standard errors for network statistics over 100 samples for sampling 1 layer via BFS. For the real data donation only one sample exists.")

pos <- seq(3, nrow(NMSE_statistics_biases_k1_deg_glc), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_statistics_biases_k1_deg_glc), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


int_headlines <- c("women more likely to donate than men","men more likely to donate than women",
                   "caucasians more likely to donate than non-caucasians","asians more likely to donate than non-asians",
                   "democrats more likely to donate than republicans","republicans more likely to donate than democrats",
                   "Random Walk as silver standard, N = 295,426")



kbl <- NMSE_statistics_biases_k2_den_clc_asort %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:NMSE_k2_biases_statistics_part1}Average normalized mean square error and standard errors for network statistics over 100 samples for sampling 2 layers via BFS. For the real data donation only one sample exists.")

pos <- seq(3, nrow(NMSE_statistics_biases_k2_den_clc_asort), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_statistics_biases_k2_den_clc_asort), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")


kbl <- NMSE_statistics_biases_k2_deg_glc %>%
  kable(format = "latex", booktabs = TRUE, linesep = "",
        caption = "\\label{tab:NMSE_k2_biases_statistics_part2}Average normalized mean square error and standard errors for network statistics over 100 samples for sampling 2 layers via BFS. For the real data donation only one sample exists.")

pos <- seq(3, nrow(NMSE_statistics_biases_k2_deg_glc), by = 3)

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] ,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_statistics_biases_k2_deg_glc), "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
}

kbl %>%
  kable_styling(latex_options = "HOLD_position")




# NMSE weighted attributes

NMSE_attributes_reweighted <- read.csv("NMSE_attributes_biases_reweighted.csv")

NMSE_attributes_reweighted <- NMSE_attributes_reweighted %>% 
  mutate(across(where(is.numeric),~ round(.x,3))) %>% 
  mutate(across(where(is.numeric),~ formatC(.x, format = "f", digits = 3)),
         across(where(is.character), ~ ifelse(.x == "0.000", ">0.001", .x)),
         across(where(is.character), ~ ifelse(grepl("NA",.x), "", .x)),
         across(ends_with("_se"), ~ ifelse(.x == "",.x,paste0(" (",.x,")")))) %>% 
  mutate(#gender
    NMSE_gender_unweighted = paste0(gender_NMSE_naive_mean_av,gender_NMSE_naive_mean_se),
    NMSE_gender_sample_indegree = paste0(gender_NMSE_hajek_sampindegree_av,gender_NMSE_hajek_sampindegree_se),
    NMSE_gender_sample_realindegree = paste0(gender_NMSE_hajek_realindegree_av,gender_NMSE_hajek_realindegree_se),
    #ethnicity
    NMSE_ethnicity_unweighted = paste0(ethnicity_NMSE_naive_mean_av,ethnicity_NMSE_naive_mean_se),
    NMSE_ethnicity_sample_indegree = paste0(ethnicity_NMSE_hajek_sampindegree_av,ethnicity_NMSE_hajek_sampindegree_se),
    NMSE_ethnicity_sample_realindegree = paste0(ethnicity_NMSE_hajek_realindegree_av,ethnicity_NMSE_hajek_realindegree_se),
    #party
    NMSE_party_unweighted = paste0(party_NMSE_naive_mean_av,party_NMSE_naive_mean_se),
    NMSE_party_sample_indegree = paste0(party_NMSE_hajek_sampindegree_av,party_NMSE_hajek_sampindegree_se),
    NMSE_party_sample_realindegree = paste0(party_NMSE_hajek_realindegree_av,party_NMSE_hajek_realindegree_se)
  ) %>% 
  select(sample_method, sample_depth, starts_with("NMSE")) %>% 
  pivot_longer(
    cols = c(starts_with("NMSE")),
    names_to = "correction",
    values_to = "NMSE"
  ) %>% 
  mutate(
    attribute = str_extract(correction, "^(([^_]*_){1}[^_]*)"),
    correction  = str_remove(correction, "^(([^_]*_){1}[^_]*)_")
  ) %>% 
  pivot_wider(
    id_cols = c(sample_method,sample_depth, correction),
    names_from = attribute,
    values_from = NMSE
  ) %>% 
  filter(!(sample_method %in% c("RW0.15"))) %>% 
  mutate(sample_method = case_match(sample_method,
                                    "BFS" ~ "no bias",
                                    "BFS_biasdegree_linear" ~ "indegree bias linear",
                                    "donated" ~"seeds from data donation",
                                    "RW0.15" ~ "RW, j = 0.15",
                                    .default = NA),
         correction = case_match(correction,
                                 "unweighted" ~ "\\hspace{1em}unweighted",
                                 "sample_indegree" ~ "\\hspace{2em}... sample indegree",
                                 "sample_realindegree" ~ "\\hspace{2em}... real indegree"))


NMSE_attributes_reweighted_k1 <- NMSE_attributes_reweighted %>% 
  filter(sample_depth == "k1") %>% 
  arrange(match(sample_method, c("no bias","seeds from data donation","indegree bias linear",
                                 "RW, j = 0.15")))%>%
  group_by(sample_method) %>%
  mutate(n = n()) %>%   # count rows per group
  ungroup() 


NMSE_attributes_reweighted_k1 <- NMSE_attributes_reweighted_k1 %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "_", " "))) %>%
  rename_with(~ str_replace_all(.x, "_", " "))

NMSE_attributes_reweighted_k2 <- NMSE_attributes_reweighted %>% 
  filter(sample_depth == "k2")  %>% 
  arrange(match(sample_method, c("no bias","seeds from data donation","indegree bias linear",
                                 "RW, j = 0.15")))%>% 
  group_by(sample_method) %>%
  mutate(n = n()) %>%   # count rows per group
  ungroup() 

NMSE_attributes_reweighted_k2 <- NMSE_attributes_reweighted_k2 %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "_", " "))) %>%
  rename_with(~ str_replace_all(.x, "_", " ")) %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "_", " "))) %>%
  rename_with(~ str_replace_all(.x, "_", " "))

int_headlines <- c("BFS, no bias","BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 10,373")


kbl <- NMSE_attributes_reweighted_k1 %>%
  select(-`sample method`, -`sample depth`, -n) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE,
        caption = "\\label{tab:reweighted_attributes_k1}Average normalized mean square error and standard errors for nodal attribute shares over 100 samples for sampling 1 layer via BFS corrected by sample indegree, sample indegree without seeds and real indegree from the complete dataset. For the real data donation only one sample exists.") 

pos <- c(0,2,4,6)
pos2 <- pos +1

offset <- 0

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] + offset,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_attributes_reweighted_k1)-3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
  
  offset <- offset + 1
  
  kbl <- kbl %>%
    row_spec(
      pos[i] +  offset,
      extra_latex_after = paste0(
        "\\multicolumn{", ncol(NMSE_attributes_reweighted_k1)-3, "}{l}{\\hspace{2em}\\textit{\\hspace{1em}weighted by ...}} \\\\"
      )
    )
  offset <- offset + 1
}


kbl %>%
  kable_styling(latex_options = "HOLD_position")


int_headlines <- c("BFS, no bias","BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 295,426")

kbl <- NMSE_attributes_reweighted_k2 %>%
  select(-`sample method`, -`sample depth`,-n) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE,
        caption = "\\label{tab:reweighted_attributes_k2}Average normalized mean square error and standard errors for nodal attribute shares over 100 samples for sampling 2 layers via BFS corrected by sample indegree, sample indegree without seeds and real indegree from the complete dataset. For the real data donation only one sample exists.")

offset <- 0

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] + offset,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_attributes_reweighted_k2)-3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
  
  offset <- offset + 1
  
  kbl <- kbl %>%
    row_spec(
      pos[i] +  offset,
      extra_latex_after = paste0(
        "\\multicolumn{", ncol(NMSE_attributes_reweighted_k2)-3, "}{l}{\\hspace{2em}\\textit{\\hspace{1em}weighted by ...}} \\\\"
      )
    )
  offset <- offset + 1
}


kbl %>%
  kable_styling(latex_options = "HOLD_position")


# NMSE weighted network statistics

NMSE_statistics_reweighted <- read.csv("NMSE_statistics_biases_reweighted.csv", sep = ",")

NMSE_statistics_reweighted <- NMSE_statistics_reweighted %>% 
  mutate(across(where(is.numeric),~ round(.x,3))) %>% 
  mutate(across(where(is.numeric),~ formatC(.x, format = "f", digits = 3)),
         across(where(is.character), ~ ifelse(.x == "0.000", ">0.001", .x)),
         across(where(is.character), ~ ifelse(grepl("NA",.x), "", .x)),
         across(ends_with("_se"), ~ ifelse(.x == "",.x,paste0(" (",.x,")")))) %>% 
  mutate(#density
    NMSE_density_unweighted = paste0(density_NMSE_naive_mean_av,density_NMSE_naive_mean_se),
    NMSE_density_sample_indegree = paste0(density_NMSE_hajek_sampindegree_av,density_NMSE_hajek_sampindegree_se),
    NMSE_density_sample_realindegree = paste0(density_NMSE_hajek_realindegree_av,density_NMSE_hajek_realindegree_se),
    #lcl  coef
    NMSE_lcl_unweighted = paste0(lcl_coef_NMSE_naive_mean_av,lcl_coef_NMSE_naive_mean_se),
    NMSE_lcl_sample_indegree = paste0(lcl_coef_NMSE_hajek_sampindegree_av,lcl_coef_NMSE_hajek_sampindegree_se),
    NMSE_lcl_sample_realindegree = paste0(lcl_coef_NMSE_hajek_realindegree_av,lcl_coef_NMSE_hajek_realindegree_se),
    #ethnic_assortativity
    NMSE_ethassort_unweighted = paste0(assortativity_ethnic_NMSE_naive_mean_av,assortativity_ethnic_NMSE_naive_mean_se),
    NMSE_ethassort_sample_indegree = paste0(assortativity_ethnic_NMSE_hajek_sampindegree_av,assortativity_ethnic_NMSE_hajek_sampindegree_se),
    NMSE_ethassort_sample_realindegree = paste0(assortativity_ethnic_NMSE_hajek_realindegree_av,assortativity_ethnic_NMSE_hajek_realindegree_se),
    #global clustering coefficient
    NMSE_gcl_unweighted = paste0(clustering_coeff_global_NMSE_naive_mean_av,clustering_coeff_global_NMSE_naive_mean_se),
    NMSE_gcl_sample_indegree = paste0(clustering_coeff_global_NMSE_hajek_sampindegree_av,clustering_coeff_global_NMSE_hajek_sampindegree_se),
    NMSE_gcl_sample_realindegree = paste0(clustering_coeff_global_NMSE_hajek_realindegree_av,clustering_coeff_global_NMSE_hajek_realindegree_se)
  ) %>%
  select(sample_method, sample_depth, starts_with("NMSE")) %>% 
  pivot_longer(
    cols = c(starts_with("NMSE")),
    names_to = "correction",
    values_to = "NMSE"
  ) %>% 
  mutate(
    attribute = str_extract(correction, "^(([^_]*_){1}[^_]*)"),
    correction  = str_remove(correction, "^(([^_]*_){1}[^_]*)_")
  ) %>% 
  pivot_wider(
    id_cols = c(sample_method,sample_depth, correction),
    names_from = attribute,
    values_from = NMSE
  ) %>% 
  mutate(sample_method = case_match(sample_method,
                                    "BFS" ~ "no bias",
                                    "BFS_biasdegree_linear" ~ "indegree bias linear",
                                    "donated" ~"seeds from data donation",
                                    "RW0.15" ~ "RW, j = 0.15",
                                    .default = NA),
         correction = case_match(correction,
                                 "unweighted" ~ "\\hspace{1em}unweighted",
                                 "sample_indegree" ~ "\\hspace{2em}... sample indegree",
                                 "sample_realindegree" ~ "\\hspace{2em}... real indegree")) %>% 
  arrange(match(sample_method, c("no bias","seeds from data donation","indegree bias linear",
                                 "RW, j = 0.15")))



NMSE_statistics_reweighted_k1 <- NMSE_statistics_reweighted %>% 
  filter(sample_depth == "k1") %>% 
  arrange(match(sample_method, c("no bias","seeds from data donation","indegree bias linear",
                                 "RW, j = 0.15"))) %>% 
  group_by(sample_method) %>%
  mutate(n = n()) %>%   # count rows per group
  ungroup() %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "_", " "))) %>%
  rename_with(~ str_replace_all(.x, "_", " ")) %>% 
  rename(`NMSE local clustering coef` = `NMSE lcl` ,
         `NMSE ethnic assortativity` = `NMSE ethassort`)

int_headlines <- c("BFS, no bias","BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 10,373")


kbl <- NMSE_statistics_reweighted_k1 %>%
  select(-`sample method`, -`sample depth`, -n) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE,
        caption = "\\label{tab:reweighed_statistics_k1_part1}Average normalized mean square error and standard errors for network statistics over 100 samples for sampling 1 layer via BFS corrected by sample indegree, sample indegree without seeds and real indegree from the complete dataset. For the real data donation only one sample exists.") 

pos <- c(0,1,2,3) #seq(0, 20, by = 4)
pos2 <- pos +1

offset <- 0

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] + offset,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_statistics_reweighted_k1)-3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
  
  offset <- offset + 1
  
  kbl <- kbl %>%
    row_spec(
      pos[i] +  offset,
      extra_latex_after = paste0(
        "\\multicolumn{", ncol(NMSE_statistics_reweighted_k1)-3, "}{l}{\\hspace{2em}\\textit{\\hspace{1em}weighted by ...}} \\\\"
      )
    )
  offset <- offset + 1
}


kbl %>%
  kable_styling(latex_options = "HOLD_position")

#k2
NMSE_statistics_reweighted_k2 <- NMSE_statistics_reweighted %>% 
  filter(sample_depth == "k2")  %>% 
  arrange(match(sample_method, c("no bias","seeds from data donation","indegree bias linear",
                                 "RW, j = 0.15")))%>%
  group_by(sample_method) %>%
  mutate(n = n()) %>%   # count rows per group
  ungroup() %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "_", " "))) %>%
  rename_with(~ str_replace_all(.x, "_", " ")) %>% 
  mutate(across(where(is.character), ~ str_replace_all(., "_", " "))) %>%
  rename_with(~ str_replace_all(.x, "_", " ")) %>% 
  rename(`NMSE local clustering coef` = `NMSE lcl` ,
         `NMSE ethnic assortativity` = `NMSE ethassort`)


int_headlines <- c("BFS, no bias","BFS, seeds from data donation",
                   "BFS, indegree bias linear",
                   "Random Walk (j = 0.15) as silver standard, N = 295,426")

kbl <- NMSE_statistics_reweighted_k2 %>%
  select(-`sample method`, -`sample depth`, -n) %>%
  kable(format = "latex", booktabs = TRUE, linesep = "", escape = FALSE,
        caption = "\\label{tab:reweighed_statistics_k2}Average normalized mean square error and standard errors for network statistics over 100 samples for sampling 2 layers via BFS corrected by sample indegree and real indegree from the complete dataset. For the real data donation only one sample exists.") 

pos <- c(0,1,2,3) 
pos2 <- pos +1

offset <- 0

for (i in seq_along(pos)) {
  kbl <- kbl %>%
    row_spec(
      pos[i] + offset,#+ (i - 1),
      extra_latex_after = paste0(
        "\\addlinespace\n",
        "\\multicolumn{", ncol(NMSE_statistics_reweighted_k2)-3, "}{l}{\\textit{", int_headlines[i], "}} \\\\"
      )
    )
  
  offset <- offset + 1
  
  kbl <- kbl %>%
    row_spec(
      pos[i] +  offset,
      extra_latex_after = paste0(
        "\\multicolumn{", ncol(NMSE_statistics_reweighted_k2)-3, "}{l}{\\hspace{2em}\\textit{\\hspace{1em}weighted by ...}} \\\\"
      )
    )
  offset <- offset + 1
}


kbl %>%
  kable_styling(latex_options = "HOLD_position")






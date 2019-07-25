#### Header ####
# Stand-alone program to run the Spatial Wallinga Teunis R effective method
# Requires:
# Case data
# Weight Matrix formatted as three columns: id_source, id_destination, weight
# Written by Corey Peak (peak@mail.harvard.edu)

# Chiefdom is the Admin 3 poplation unit at which level the case count data are aggregated by day.

#### Helper Function: make_arrays ####
make_arrays <- function(id.var, date.var){
  
  # Numerator array dimensions are row:infected, column:infector, z:chiefdom
  numerator_array <- array(rep(0, (length(unique(date.var)))^2*length(unique(id.var))), dim=c(length(unique(date.var)), length(unique(date.var)), length(unique(id.var))))
  
  # Denominator array dimensions are row:infected, column:nothing, z:chiefdom
  denominator_array <- array(rep(0, (length(unique(date.var)))*length(unique(id.var))), dim=c(length(unique(date.var)), 1, length(unique(id.var))))
  
  # Quotient array dimensions are row:infected, column:infector, z:chiefdom
  quotient_array <- numerator_array
  
  # Make a key for each chiefdom matrix
  id_array_key <- data.frame(id = unique(id.var), key = seq(1:length(unique(id.var))))
  
  return(list(numerator_array = numerator_array,
              denominator_array = denominator_array, 
              quotient_array = quotient_array, 
              id_array_key = id_array_key))
}

#### Calculate quotient array ####
WT_Spatial <- function(id.var, date.var, case.var, weights, generation_interval, subset_ids = "all"){
  require(dplyr)
  
  Break = 0
  # Troubleshoot
  if (length(id.var) != length(date.var)){Break = Break + 1; cat("\nError 1")}
  if (length(id.var) != length(case.var)){Break = Break + 1; cat("\nError 2")}
  if (length(unique(id.var))^2 != nrow(weights)){Break = Break + 1; cat("\nError 3")}
  if (sum(id.var %in% weights[,1]) != length(id.var)){Break = Break + 1; cat("\nError 4")}
  if (sum(id.var %in% weights[,2]) != length(id.var)){Break = Break + 1; cat("\nError 5")}
  
  generation_interval_cumsum <- cumsum(generation_interval)
  max_GT <- which(generation_interval_cumsum > 0.99)[1]
  generation_interval <- generation_interval[1:(max_GT+1)]
  
  if (Break == 0){
    
    df <- data.frame(id = id.var, date = date.var, case = case.var)
    df$id <- as.character(id.var)
    names(weights) <- c("id_source", "id_destination", "weight")
    
    # Remove units with no cases ever
    df_total <- df %>% group_by(id) %>% summarize(total_cases = sum(case)) %>% filter(total_cases > 0)
    df <- df %>% filter(id %in% df_total$id)
    
    if (("all" %in% subset_ids) != 1){
      df <- df[df$id %in% subset_ids,]
      id.var <- id.var[id.var %in% subset_ids]
      if (length(subset_ids) > 1){
        weights <- weights[weights$id_source %in% subset_ids & weights$id_destination %in% subset_ids,]
      } else { # If only looking at a single location, then you can compare the model performance to the cholera_df Wallinga-Teunis output
        weights <- data.frame(id_source = subset_ids, id_destination = subset_ids, weight = 1)
      }
    }
    
    # Make arrays (See Helper Function)
    arrays <- make_arrays(id.var, date.var)
    numerator_array <- arrays$numerator_array
    denominator_array <- arrays$denominator_array
    quotient_array <- arrays$quotient_array
    id_array_key <- arrays$id_array_key
    
    # Keep track of internal Reff
    internal_Reff_array <- denominator_array
    
    first_day <- min(date.var)
    
    # Check to make sure the weights file is in the right format
    if (sum(weights$id_source != rep(id_array_key$id, each = length(id_array_key$id)))>0){"Error with weight matrix"}
    if (sum(weights$id_destination != rep(id_array_key$id, length(id_array_key$id)))>0){"Error with weight matrix"}
    
    for (infected_id in unique(df$id)){
      cat("\n", infected_id)
      infected_key <- id_array_key[id_array_key$id == infected_id,"key"]
      
      # First work through one unit
      subset = df %>% filter(id %in% infected_id) 
      weights_subset = weights %>% filter(id_destination %in% infected_id)
      
      for (i in 1:nrow(subset)) {                  # For each day in that unit
        if (subset[i, "case"] > 0){                # If someone in that unit was infected on that day
          min_j <- max(1, (i-max_GT))  # Only iterate through days with meaningful generation interval probabilities
          for (j in min_j:i) {                         # For each day when that(those) infection(s) could have occured
            infector_candidates_on_day_j = df[df$date == subset[j,"date"] & df$case > 0,]
            if (nrow(infector_candidates_on_day_j) > 0){
              candidate_ids <- as.character(infector_candidates_on_day_j[, "id"])
              candidate_keys <- as.numeric(sapply(candidate_ids, function(x) which(id_array_key$id == x)))
              generation_interval_value <- as.numeric(generation_interval[subset[i,"date"]-subset[j,"date"]+1] )
                # The "money" given by somone on day i to people on day j in candidate_id
                numerator_array[i,j,candidate_keys] <- infector_candidates_on_day_j$case * generation_interval_value * weights_subset[candidate_keys,"weight"]
            }
          } #Now we've gone through each candidate infector
          denominator_array[i,1,infected_key] <- sum(numerator_array[i,,]) # The total "money" given by the people infected on day i in id
          if (denominator_array[i,1,infected_key] > 0){
            quotient_array[i,,] = quotient_array[i,,] + numerator_array[i,,] * subset[i, "case"]  / denominator_array[i,1,infected_key]  # "money" given to candidates for infecting people on day i in id
            internal_Reff_array[,,infected_key] <- internal_Reff_array[,,infected_key] + numerator_array[i,,infected_key] * subset[i, "case"]  / denominator_array[i,1,infected_key]
          }
          numerator_array[,,] <- 0 
          cat(".")
        }
      }
    }
  }
  
  df_out = calculate_Reff(df, quotient_array, internal_Reff_array, id_array_key)
  return(df_out)
}

#### Calculate Reff from quotient_array ####
calculate_Reff <- function(df, quotient_array, internal_Reff_array, id_array_key){
  cat("\nCalculating R effective")
  #initialize R-effective
  df$Reff <- NA
  df$Reff_internal <- NA
  # df$Reff_imported <- NA
  df$Reff_external <- NA
  
  for (id in df$id){
    id_key <- id_array_key[id_array_key$id == id,"key"]
    
    # Calculate Reff
    df[df$id == id, "Reff"] <- apply(quotient_array[,,id_key], 2, sum ) # Total infections attributed to that day
    df[df$id == id, "Reff"] <- df[df$id == id, "Reff"] / df[df$id == id, "case"]
    df[df$id == id & is.na(df$Reff)==1, "Reff"] <- 0
    
    # Calculate Reff_internal
    df[df$id == id, "Reff_internal"] <- internal_Reff_array[,,id_key] / df[df$id == id, "case"]
    df[df$id == id & is.na(df$Reff_internal)==1, "Reff_internal"] <- 0
    
    # cat("\nID Number", id)
  }
  # Calculate Reff_external
  df[,"Reff_external"]  <- round(df[,"Reff"] -  df[,"Reff_internal"], 6)
  
  df$exported_cases <- df$case * df$Reff_external
  
  return(df)
}

#### Summary Function ####
WT_Spatial_Summary <- function(df){
  df %>% group_by(id) %>%
    summarize(days_Reff_over_1 = sum(Reff > 1),
              total_cases = sum(case),
              total_cases_generated = sum(case*Reff),
              total_exported_cases = sum(exported_cases),
              total_internal_cases_generated = sum(case*Reff_internal),
              total_imported_cases = total_cases - total_internal_cases_generated,
              export_import_ratio = total_exported_cases / (0.001+total_imported_cases))
}

#### To practice, here are some play data ####
days <- 1:6
id.var_play <- rep(c("A", "B"), each = length(days))
date.var_play <- rep(days, length(unique(id.var_play)))
case.var_play <- c(0,20,40,80,20,0,
                   0,5,0,1,0,0)

tau_play <- c(0, 0.3, 0.3, 0.2, 0.1, 0.1, 0)
plot(tau_play, type = "l", xlab = "Days", ylab = "Generation Interval Density")

weight_play <- data.frame(id_source = c("A","A","B","B"), id_destination = c("A", "B", "A", "B"), weight = c(1,0.5,0.5,1))
weight_play

play_output <- WT_Spatial(id.var = id.var_play, date.var = date.var_play, case.var = case.var_play, weights = weight_play, generation_interval = tau_play)
play_output

WT_Spatial_Summary(play_output)

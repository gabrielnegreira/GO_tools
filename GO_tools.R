#get required libraries
library(httr)
library(jsonlite)
library(xml2)
library(rlist)
library(dplyr)
library(XML)
library(readr)

#clean_GO_terms####
#this function will take a vector of GO_terms and use the quickGO api to retrieve updated terms for them.
clean_GO_terms <- function(GO_terms){
  #remove NAs and duplicates
  GO_terms <- GO_terms[!is.na(GO_terms)]
  if(length(GO_terms) != length(unique(GO_terms))){
    warning("found multiple occurences of identical terms. Only one occurence will be considered.")
    GO_terms <- unique(GO_terms)
  }
  
  #create the dataframe where the new information will be stored
  final_data <- data.frame(input_term = GO_terms)
  #setting row names to the input term makes it easier to map results later to the correct rows in the dataframe.
  rownames(final_data) <- final_data$input_term
  
  #create a single string containing the GO terms in the expected format for the URL request
  string_to_search <- gsub("\\:", "%3A", GO_terms)
  
  #the quick GO API accepts a maximum of 525 terms per query, so we'll break the vector to run multiple queries.
  breaks <- seq(0, length(string_to_search), 525)
  breaks <- unique(c(breaks, length(string_to_search)))
  res_list <- list()
  for(i in c(2:length(breaks))){
    #create the request URL for the complete information about the GO terms
    indexes <- c(breaks[i-1]:breaks[i])
    requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", paste0(string_to_search[indexes], collapse = "%2C"), "/complete")
    #make the request
    r <- GET(requestURL, accept("application/json"))
    #check for errors
    stop_for_status(r)
    #convert the JSON formatted data to a list
    res <- content(r)$results
    #append the results to the final list with all information per term
    indexes <- c((length(res_list)+1):(length(res_list)+length(res)))
    res_list[indexes] <- res
  }
  
  #set the names of the list objects to the GO term id
  names(res_list) <- sapply(res_list, function(x) x$id)
  
  #update the results dataframe with the information about whether the term is in the quickGO database
  final_data$is_in_db <- final_data$input_term %in% unique(c(names(res_list), as.character(unlist(sapply(res_list, function(x) x$secondaryIds)))))
  #now we have to find which terms are secondary ids
  final_data$is_secondary <- final_data$input_term %in% as.character(unlist(sapply(res_list, function(x) x$secondaryIds)))
  
  #find the primary ids of the secondary terms
  final_data$primary_term <- NA
  for(secondary_term in filter(final_data, is_secondary)$input_term){
    index <- which(grepl(secondary_term, sapply(res_list, function(x) x$secondaryIds)))
    final_data[secondary_term,]$primary_term <- unique(names(res_list[index]))
  }
  
  #defined the optimal term for each input term
  final_data <- final_data %>%
    mutate(db_term = ifelse(is.na(primary_term), input_term, primary_term),
           db_term_name = as.character(sapply(res_list, function(x) x$name)[db_term]),
           is_obsolete = as.logical(as.character(sapply(res_list, function(x) x$isObsolete)[db_term])),
           db_term_aspect = as.character(sapply(res_list, function(x) x$aspect)[db_term]),
           replaced_by = ifelse(is_obsolete, as.character(sapply(res_list, function(x){
             x <- x$replacements
             x <- list.filter(x, type == "replaced_by")
             x <- sapply(x, function(x) as.character(x$id))
           })[db_term]), NA),
           replaced_by = ifelse(replaced_by == "list()", NA, replaced_by),
           best_term = ifelse(is_obsolete & !is.na(replaced_by), replaced_by, 
                              ifelse(is_in_db, db_term, input_term)))
  
  #now that we have defined the best terms, we need some information about these terms which are not present in the input_term. Thus, we make a new query.
  GO_terms <- final_data %>%
    filter(best_term != db_term) %>%
    pull("best_term")
  
  #create a single string containing the GO terms in the expected format for the URL request
  string_to_search <- gsub("\\:", "%3A", GO_terms)
  
  #the quick GO API accepts a maximum of 525 terms per query, so we'll break the vector to run multiple queries.
  breaks <- seq(0, length(string_to_search), 525)
  breaks <- unique(c(breaks, length(string_to_search)))
  res_list <- list()
  for(i in c(2:length(breaks))){
    #create the request URL for the complete information about the GO terms
    indexes <- c(breaks[i-1]:breaks[i])
    requestURL <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", paste0(string_to_search[indexes], collapse = "%2C"), "/complete")
    #make the request
    r <- GET(requestURL, accept("application/json"))
    #check for errors
    stop_for_status(r)
    #convert the JSON formatted data to a list
    res <- content(r)$results
    #append the results to the final list will all information per term
    indexes <- c((length(res_list)+1):(length(res_list)+length(res)))
    res_list[indexes] <- res
  }
  
  #set the names of the list objects to the GO term id
  names(res_list) <- sapply(res_list, function(x) x$id)
  
  #append the new information
  final_data <- final_data %>%
    mutate(best_term_name = ifelse(best_term %in% names(res_list), as.character(sapply(res_list, function(x) x$name)[best_term]), db_term_name),
           best_term_aspect = ifelse(best_term %in% names(res_list), as.character(sapply(res_list, function(x) x$aspect)[best_term]), db_term_aspect))
  
  return(final_data)
}

#revigo_query####
#this function uses the API of revigo (http://revigo.irb.hr/FAQ) to summarize a list of GO terms into simpler terms. 
#it accepts as input both a vector of GOterms, or a data frame with two columns, the first being the GOterms and the second being the value (like pvalue for instance)
revigo_query <- function(GO_terms, cutoff = "0.7", valueType = c("PValue", "Higher", "Lower", "HigherAbsolute", "HigherAbsLog2"), speciesTaxon = "0", measure = c("SIMREL", "LIN", "RESNIK", "JIANG"), removeObsolete = TRUE){
  
  #check function parameters
  valueType <- match.arg(valueType)
  measure <- match.arg(measure)
  
  #first convert the provided table to a temporary tsv file
  write_tsv(as.data.frame(GO_terms), file = "temp.tsv")
  #now load the tsv file as is
  input <- readChar("temp.tsv",file.info("temp.tsv")$size)
  file.remove("temp.tsv")
  #now use the RESTful API to query the results (for more information check: http://revigo.irb.hr/FAQ)
  print("Making request to RESTful API of revigo (this might take a while) ...")
  results <- httr::POST(
    url = "http://revigo.irb.hr/Revigo",
    body = list(
      cutoff = cutoff,
      valueType = valueType,
      speciesTaxon = speciesTaxon,
      measure = measure,
      goList = input,
      removeObsolete = tolower(as.character(removeObsolete))
    ),
    # application/x-www-form-urlencoded
    encode = "form"
  )
  print("done!")
  #convert the html output to a dataframe
  dat <- httr::content(results, encoding = "UTF-8")
  dat <- readHTMLTable(htmlParse(dat))
  
  #fix the go_type names to match other databases such as the quickGO
  names(dat) <- c("MolecularFunction" = "molecular_function", "CellularComponent" = "cellular_component", "BiologicalProcess" = "biological_process")[names(dat)]
  
  #at this point dat is a list with 3 data frames, one for biological processes, one for cellular component, and one for molecular function.
  #here it will merge these 3 dataframes into a single one. 
  dat_df <- c()
  for(i in c(1:length(dat))){
    df <- dat[[i]]
    df$go_type <- names(dat)[i]
    dat_df <- rbind(dat_df, df)
  }
  #formating the column names
  colnames(dat_df) <- colnames(dat_df) %>%
    gsub(" ", "_", .) %>%
    tolower()
  
  #finally it will, for each GO term, add a column indicating its representative term
  #OBS: parental terms have a "null" string in the 'Representative' column.
  dat_df$repr_id <- dat_df$term_id
  dat_df$repr_name <- dat_df$name
  
  for(i in c(1:nrow(dat_df))){
    if(tolower(dat_df$representative[i]) == "null"){
      next
    }else{
      dat_df$repr_id[i] <- dat_df$repr_id[i-1]
      dat_df$repr_name[i] <- dat_df$repr_name[i-1]
    }
  }
  rownames(dat_df) <- dat_df$term_id
  return(dat_df)
}

#GO_enrich####
#function to perform a simple GO enrichment analysis using hypergeometric distribution statistics.
#it expects two vectors, one containing the GO_terms found among the differential expressed genes (set_terms).
#and another vector containing all the GO terms found in the reference genome (ref_terms).
GO_enrich <- function(set_terms, ref_terms, test_type = c("enrichment", "depletion", "both")){
  test_type <- match.arg(test_type)
  test_type <- c("enrichment" = "greater", "depletion" = "less", "both" = "two.sided")[test_type]
  results <- as.data.frame(sort(table(set_terms), decreasing = TRUE))
  results$pvalue <- NA
  results$count_in_ref <- NA
  for(i in c(1:nrow(results))){
    #get the term as well as the times it is found in the set and in the reference
    term <- results$set_terms[i]
    x <- sum(set_terms == term)
    m <- sum(ref_terms == term)
    k <- length(set_terms)
    n <- length(ref_terms) - m
    
    #build the contingency matrix for fisher test
    cont_matrix <- matrix(ncol = 2, nrow = 2)
    cont_matrix[1,1] <- x
    cont_matrix[2,1] <- k-x
    cont_matrix[1,2] <- m-x
    cont_matrix[2,2] <- n-(k-x)
    
    #perform_fisher test
    results$pvalue[i] <- fisher.test(cont_matrix, alternative = test_type)$p.value
    
    #add extra info
    results$count_in_ref[i] <- m
  }
  colnames(results)[c(1:2)] <- c("term", "count_in_set")
  results <- results %>%
    mutate(prop_in_set = count_in_set/sum(count_in_set),
           prop_in_ref = count_in_ref/length(ref_terms),
           log2enrichment = log2(prop_in_set/prop_in_ref))
  return(results)
}

#load tidyverse libraries
#install.packages(c("httr", "jsonlite", "xml2", "rlist", "dplyr", "XML", "readr")) #uncoment if needed.
library(dplyr)

#download GO tools from github repo
#unlink("GO_tools", recursive = TRUE)
#system("git clone https://github.com/gabrielnegreira/GO_tools.git")
#source it
source("GO_tools/GO_tools.R")


#import example file
go_df <- read.table("GO_tools/example.txt", comment.char = "%", header = FALSE, sep = "", col.names = c("GO_ID", "p_value"), stringsAsFactors = FALSE)
go_df %>%
  head() %>%
  knitr::kable(format = "markdown", row.names = FALSE)

#run `clean_GO_terms()` and print only the terms where the suggested terms differ from the input terms.
cleaned_terms <- clean_GO_terms(go_df$GO_ID) 
cleaned_terms %>%
  filter(input_term != best_term) %>%
  knitr::kable(format = "markdown", row.names = FALSE)

go_df <- go_df %>%
  left_join(cleaned_terms %>% 
              rename(GO_ID = input_term, cleaned_GO_ID = best_term) %>% 
              select(GO_ID, cleaned_GO_ID), 
            by = "GO_ID")

go_df %>%
  filter(GO_ID != cleaned_GO_ID) %>%
  knitr::kable(format = "markdown", row.names = FALSE)

#run `revigo_query()`
reduced_terms <- revigo_query(go_df$cleaned_GO_ID)
reduced_terms  %>%
  filter(term_id != repr_id) %>%
  head() %>%
  knitr::kable(format = "markdown", row.names = FALSE)

#replace original terms by the reduced ones
go_df <- go_df %>%
  left_join(reduced_terms %>% rename(cleaned_GO_ID = term_id, reduced_GO_ID = repr_id) %>% select(cleaned_GO_ID, reduced_GO_ID), by = "cleaned_GO_ID")

#arrange by reduced terms
go_df <- go_df %>%
  filter(cleaned_GO_ID != reduced_GO_ID) %>%
  group_by(reduced_GO_ID) %>%
  mutate(n_terms = n()) %>%
  arrange(desc(n_terms)) %>%
  select(-n_terms)

#print table  
go_df %>%
  head() %>%
  knitr::kable(format = "markdown", row.names = FALSE)

#run get_GO_data
go_meta <- get_GO_data(unique(go_df$reduced_GO_ID))
go_meta %>%
  head() %>%
  knitr::kable(format = "markdown", row.names = FALSE)
  
#append name to go_df
go_df <- go_df %>%
  left_join(go_meta %>%
              rename(reduced_GO_ID = input_term, reduced_GO_name = name, reduced_GO_type = aspect) %>%
              select(reduced_GO_ID, reduced_GO_name, reduced_GO_type),
            by = "reduced_GO_ID")

#print table  
go_df %>%
  head() %>%
  knitr::kable(format = "markdown", row.names = FALSE)

#perform GO enrichment analysis
set_terms <- go_df$reduced_GO_ID[1:6]
  
ref_terms <- go_df$reduced_GO_ID

GO_enrich(set_terms, ref_terms)%>%
  knitr::kable(format = "markdown", row.names = FALSE)


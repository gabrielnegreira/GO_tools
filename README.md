# GO_tools.R

The `GO_tools.R` script provides a suite of functions for processing and analyzing Gene Ontology (GO) terms in R. It offers utilities to clean, update, and retrieve detailed metadata for GO terms via the [QuickGO API](https://www.ebi.ac.uk/QuickGO/), summarize GO terms via the [ReviGO API](http://revigo.irb.hr/), and perform GO enrichment analysis using hypergeometric statistics.

## Usage
1) Install required packages:
```r
install.packages(c("httr", "jsonlite", "xml2", "rlist", "dplyr", "XML", "readr"))
```
2) Clone the repository and source the `GO_tools.R` script from this repository:
```r
#download GO tools from github repo
system("git clone https://github.com/gabrielnegreira/GO_tools.git")
#source it
source("GO_tools/GO_tools.R")
```
We can then append metadata of interest to our main dataframe:
```r
```

This will load the functions to the environment.
## Overview of Functions

### 1. `clean_GO_terms()`

This function takes a vector of GO terms and uses the [QuickGO API](https://www.ebi.ac.uk/QuickGO/) to look for updated terms. It identifies whether each term is primary or secondary, and for secondary terms, it determines the corresponding primary term. It also flags obsolete terms and suggest replacements for them when applicable. It returns a dataframe with these information, including the `best_term` column which contains a suggested replacement for the input term. 

**Usage Example:**
1) `clean_GO_terms()` only requires a vector of GO terms as input:
```r
#inport example file
go_df <- read.table("GO_tools/example.txt", comment.char = "%", header = FALSE, sep = "", col.names = c("GO_ID", "p_value"), stringsAsFactors = FALSE)

head(go_df)
```
|GO_ID      |   p_value|
|:----------|---------:|
|GO:0040017 | 0.0001894|
|GO:0051272 | 0.0001894|
|GO:0016049 | 0.0019016|
|GO:0008361 | 0.0022678|
|GO:0007091 | 0.0023467|
|GO:0030335 | 0.0029797|

```r
cleaned_terms <- clean_GO_terms(go_df$GO_ID) 
cleaned_terms %>%
  filter(input_term != best_term)
```
The result looks like this:
|input_term |is_in_db |is_secondary |primary_term |db_term    |db_term_name                                       |is_obsolete |db_term_aspect     |replaced_by |best_term  |best_term_name                                     |best_term_aspect   |
|:----------|:--------|:------------|:------------|:----------|:--------------------------------------------------|:-----------|:------------------|:-----------|:----------|:--------------------------------------------------|:------------------|
|GO:0006465 |TRUE     |FALSE        |NA           |GO:0006465 |obsolete signal peptide processing                 |TRUE        |biological_process |GO:0051604  |GO:0051604 |protein maturation                                 |biological_process |
|GO:0043193 |TRUE     |TRUE         |GO:0045893   |GO:0045893 |positive regulation of DNA-templated transcription |FALSE       |biological_process |NA          |GO:0045893 |positive regulation of DNA-templated transcription |biological_process |
|GO:0030384 |TRUE     |TRUE         |GO:0046488   |GO:0046488 |phosphatidylinositol metabolic process             |FALSE       |biological_process |NA          |GO:0046488 |phosphatidylinositol metabolic process             |biological_process |
|GO:0015893 |TRUE     |TRUE         |GO:0042908   |GO:0042908 |xenobiotic transport                               |FALSE       |biological_process |NA          |GO:0042908 |xenobiotic transport                               |biological_process |
|GO:0046489 |TRUE     |TRUE         |GO:0006661   |GO:0006661 |phosphatidylinositol biosynthetic process          |FALSE       |biological_process |NA          |GO:0006661 |phosphatidylinositol biosynthetic process          |biological_process |
|GO:0005830 |TRUE     |TRUE         |GO:0022626   |GO:0022626 |cytosolic ribosome                                 |FALSE       |cellular_component |NA          |GO:0022626 |cytosolic ribosome                                 |cellular_component |
|GO:0004289 |TRUE     |FALSE        |NA           |GO:0004289 |obsolete subtilase activity                        |TRUE        |molecular_function |GO:0004252  |GO:0004252 |serine-type endopeptidase activity                 |molecular_function |
|GO:0008022 |TRUE     |FALSE        |NA           |GO:0008022 |obsolete protein C-terminus binding                |TRUE        |molecular_function |GO:0005515  |GO:0005515 |protein binding                                    |molecular_function |

We can then replace the terms by the suggested ones:
```R
go_df <- go_df %>%
  left_join(cleaned_terms %>% rename(GO_ID = input_term) %>% select(GO_ID, best_term), by = "GO_ID")

go_df %>%
  filter(GO_ID != best_term)
```
|GO_ID      |   p_value|cleaned_GO_ID |
|:----------|---------:|:-------------|
|GO:0006465 | 0.0049359|GO:0051604    |
|GO:0043193 | 0.0066019|GO:0045893    |
|GO:0030384 | 0.0142193|GO:0046488    |
|GO:0015893 | 0.0184365|GO:0042908    |
|GO:0046489 | 0.0425017|GO:0006661    |
|GO:0005830 | 0.0338005|GO:0022626    |
|GO:0004289 | 0.0108298|GO:0004252    |
|GO:0008022 | 0.0265540|GO:0005515    |

### 2. `revigo_query()`

This function summarizes GO terms into simpler, higher representative terms using the [ReviGO API](http://revigo.irb.hr/). This function can take either a vector of GO terms or a data frame with GO terms and associated values (e.g., p-values).

**Parameters:**
- `cutoff`: Similarity cutoff for grouping terms (default is "0.7").
- `valueType`: Type of the provided value (default "PValue"). Other options include "Higher", "Lower", "HigherAbsolute", and "HigherAbsLog2".
- `speciesTaxon`: NCBI taxon ID (default "0").
- `measure`: Similarity measure to use; options include "SIMREL", "LIN", "RESNIK", "JIANG".
- `removeObsolete`: Logical flag indicating whether to remove obsolete terms (default is TRUE).

See the [Revigo documentation](http://revigo.irb.hr/FAQ) for more details.

**Usage**
```r
#run `revigo_query()`
reduced_terms <- revigo_query(go_df$cleaned_GO_ID)
reduced_terms  %>%
  filter(term_id != repr_id) %>%
  head()
```
The result looks like this:
|term_id    |name                                      |frequency |value |uniqueness |dispensability |eliminated |representative |go_type            |repr_id    |repr_name                              |
|:----------|:-----------------------------------------|:---------|:-----|:----------|:--------------|:----------|:--------------|:------------------|:----------|:--------------------------------------|
|GO:0030334 |regulation of cell migration              |0.355%    |NaN   |0.73       |0.87           |True       |40013          |biological_process |GO:0040013 |negative regulation of locomotion      |
|GO:0030335 |positive regulation of cell migration     |0.216%    |NaN   |0.70       |0.91           |True       |40013          |biological_process |GO:0040013 |negative regulation of locomotion      |
|GO:0030336 |negative regulation of cell migration     |0.064%    |NaN   |0.75       |0.76           |True       |40013          |biological_process |GO:0040013 |negative regulation of locomotion      |
|GO:0040017 |positive regulation of locomotion         |0.240%    |NaN   |0.71       |0.91           |True       |40013          |biological_process |GO:0040013 |negative regulation of locomotion      |
|GO:0006661 |phosphatidylinositol biosynthetic process |0.457%    |NaN   |0.91       |0.85           |True       |46488          |biological_process |GO:0046488 |phosphatidylinositol metabolic process |
|GO:0001710 |mesodermal cell fate commitment           |0.005%    |NaN   |0.83       |0.85           |True       |48333          |biological_process |GO:0048333 |mesodermal cell differentiation        |

To use it to replace terms by their reduced alternatives, we can do this:
```R
#arrange by reduced terms
go_df <- go_df %>%
  filter(cleaned_GO_ID != reduced_GO_ID) %>%
  group_by(reduced_GO_ID) %>%
  mutate(n_terms = n()) %>%
  arrange(desc(n_terms)) %>%
  select(-n_terms)

#print table  
go_df %>%
  head()
```

|GO_ID      |   p_value|cleaned_GO_ID |reduced_GO_ID |
|:----------|---------:|:-------------|:-------------|
|GO:0040017 | 0.0001894|GO:0040017    |GO:0040013    |
|GO:0030335 | 0.0029797|GO:0030335    |GO:0040013    |
|GO:0030334 | 0.0140864|GO:0030334    |GO:0040013    |
|GO:0030336 | 0.0321750|GO:0030336    |GO:0040013    |
|GO:0004861 | 0.0002106|GO:0004861    |GO:0004867    |
|GO:0004860 | 0.0108788|GO:0004860    |GO:0004867    |

Notice how the `GO:0040017`, `GO:0030335`, `GO:0030334`, and `GO:0030336` terms were colapsed to the `GO:0040013` term. 

### 3. `get_GO_data`
  
This function fetches additional metadata for a given vector of GO terms by querying the QuickGO API.

**Usage**
```r
go_meta <- get_GO_data(unique(go_df$reduced_GO_ID))
go_meta %>%
  head()
```
Result:
|input_term |is_in_db |is_secondary |primary_term |db_term    |name                                         |is_obsolete |aspect             |replaced_by |
|:----------|:--------|:------------|:------------|:----------|:--------------------------------------------|:-----------|:------------------|:-----------|
|GO:0040013 |TRUE     |FALSE        |NA           |GO:0040013 |negative regulation of locomotion            |FALSE       |biological_process |NA          |
|GO:0004867 |TRUE     |FALSE        |NA           |GO:0004867 |serine-type endopeptidase inhibitor activity |FALSE       |molecular_function |NA          |
|GO:0050678 |TRUE     |FALSE        |NA           |GO:0050678 |regulation of epithelial cell proliferation  |FALSE       |biological_process |NA          |
|GO:0030134 |TRUE     |FALSE        |NA           |GO:0030134 |COPII-coated ER to Golgi transport vesicle   |FALSE       |cellular_component |NA          |
|GO:0005254 |TRUE     |FALSE        |NA           |GO:0005254 |chloride channel activity                    |FALSE       |molecular_function |NA          |
|GO:0016209 |TRUE     |FALSE        |NA           |GO:0016209 |antioxidant activity                         |FALSE       |molecular_function |NA          |


We can then use it to append metadata of interest to our main dataframe:
```R
#append name to go_df
go_df <- go_df %>%
  left_join(go_meta %>%
              rename(reduced_GO_ID = input_term, reduced_GO_name = name, reduced_GO_type = aspect) %>%
              select(reduced_GO_ID, reduced_GO_name, reduced_GO_type),
            by = "reduced_GO_ID")

go_df %>%
  head()
```
|GO_ID      |   p_value|cleaned_GO_ID |reduced_GO_ID |reduced_GO_name                              |reduced_GO_type    |
|:----------|---------:|:-------------|:-------------|:--------------------------------------------|:------------------|
|GO:0040017 | 0.0001894|GO:0040017    |GO:0040013    |negative regulation of locomotion            |biological_process |
|GO:0030335 | 0.0029797|GO:0030335    |GO:0040013    |negative regulation of locomotion            |biological_process |
|GO:0030334 | 0.0140864|GO:0030334    |GO:0040013    |negative regulation of locomotion            |biological_process |
|GO:0030336 | 0.0321750|GO:0030336    |GO:0040013    |negative regulation of locomotion            |biological_process |
|GO:0004861 | 0.0002106|GO:0004861    |GO:0004867    |serine-type endopeptidase inhibitor activity |molecular_function |
|GO:0004860 | 0.0108788|GO:0004860    |GO:0004867    |serine-type endopeptidase inhibitor activity |molecular_function |

### 4. `GO_enrich`
  
This function performs a GO enrichment analysis using hypergeometric statistics and Fisher's exact test. It compares a set of GO terms (e.g., from differentially expressed genes) against a reference set (e.g., the entire genome).

**Parameters:**
- `set_terms`: A vector of GO terms from the test set.
- `ref_terms`: A vector of GO terms from the reference set.
- `test_type`: Specifies the alternative hypothesis. Options are:
- `"enrichment"` (uses "greater"),
- `"depletion"` (uses "less"),
- `"both"` (uses "two.sided").

**Usage:**  
```R
#perform GO enrichment analysis
set_terms <- go_df$reduced_GO_ID[1:6]
ref_terms <- go_df$reduced_GO_ID
GO_enrich(set_terms, ref_terms)

```
|term       | count_in_set|    pvalue| count_in_ref| prop_in_set| prop_in_ref| log2enrichment|
|:----------|------------:|---------:|------------:|-----------:|-----------:|--------------:|
|GO:0040013 |            4| 0.0020506|            4|   0.6666667|   0.1818182|       1.874469|
|GO:0004867 |            2| 0.1688312|            3|   0.3333333|   0.1363636|       1.289507|


## Notes

- **Internet Connection:**  
  Both the QuickGO and Revigo APIs require an active internet connection.
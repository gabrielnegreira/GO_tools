# GO_tools.R

The `GO_tools.R` script provides a suite of functions for processing and analyzing Gene Ontology (GO) terms. It offers utilities to clean and update GO term information via the QuickGO API, summarize GO terms using the Revigo API, retrieve detailed metadata for GO terms, and perform GO enrichment analysis using hypergeometric statistics.

## Requirements

Before using this script, ensure that the following R packages are installed:

- **httr**
- **jsonlite**
- **xml2**
- **rlist**
- **dplyr**
- **XML**
- **readr**

You can install these packages using the following command in R:

```r
install.packages(c("httr", "jsonlite", "xml2", "rlist", "dplyr", "XML", "readr"))
```

## Overview of Functions

### 1. `clean_GO_terms`

**Purpose:**  
Cleans a vector of GO terms by removing missing values and duplicates, then retrieves updated term information from the QuickGO API. It identifies whether each term is primary or secondary, and for secondary terms, it determines the corresponding primary term.

**Key Steps:**
- Removes NA values and duplicate GO terms.
- Splits queries into chunks (max 525 terms) for the QuickGO API.
- Retrieves detailed term information (e.g., name, aspect, obsolete status).
- Updates each term to its optimal representation.

**Usage Example:**

### 2. `revigo_query`

**Purpose:**  
Summarizes a list of GO terms into simpler, representative terms using the Revigo API. This function can take either a vector of GO terms or a data frame with GO terms and associated values (e.g., p-values).

**Parameters:**
- `cutoff`: Similarity cutoff for grouping terms (default is "0.7").
- `valueType`: Type of the provided value (default "PValue"). Other options include "Higher", "Lower", "HigherAbsolute", and "HigherAbsLog2".
- `speciesTaxon`: NCBI taxon ID (default "0").
- `measure`: Similarity measure to use; options include "SIMREL", "LIN", "RESNIK", "JIANG".
- `removeObsolete`: Logical flag indicating whether to remove obsolete terms (default is TRUE).

### 3. `get_GO_data`

**Purpose:**  
Fetches additional metadata for a given vector of GO terms by querying the QuickGO API. Similar to `clean_GO_terms`, it handles missing values and duplicates, and returns detailed information such as term name, aspect, and obsolete status.

### 4. `GO_enrich`

**Purpose:**  
Performs a GO enrichment analysis using hypergeometric statistics and Fisher's exact test. It compares a set of GO terms (e.g., from differentially expressed genes) against a reference set (e.g., the entire genome).

**Parameters:**
- `set_terms`: A vector of GO terms from the test set.
- `ref_terms`: A vector of GO terms from the reference set.
- `test_type`: Specifies the alternative hypothesis. Options are:
- `"enrichment"` (uses "greater"),
- `"depletion"` (uses "less"),
- `"both"` (uses "two.sided").

**Returns:**  
A data frame with:
- Term frequencies in both the test and reference sets.
- Proportions in both the test and reference sets.
- Log2 enrichment ratios.
- P-values from Fisher's test.

## Notes

- **Internet Connection:**  
  Both the QuickGO and Revigo APIs require an active internet connection.

- **API Limitations:**  
  The functions are designed to handle API limitations by chunking large queries (maximum 525 terms per query).

- **Temporary Files:**  
  The `revigo_query` function writes a temporary TSV file (`temp.tsv`) to disk. Ensure you have appropriate write permissions.

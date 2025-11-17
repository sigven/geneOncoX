# Function that retrieves geneOncoX data from Google Drive

Function that retrieves geneOncoX data from Google Drive

## Usage

``` r
get_gox_data(cache_dir = NA, force_download = FALSE, db = "alias")
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

- db:

  type of dataset to be retrieved

## Value

pre-processed gene/transcript annotation records and metadata

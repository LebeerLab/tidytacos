# Can create test taco

    Code
      tt
    Output
      $counts
      # A tibble: 6 x 3
        count sample_id taxon_id
        <dbl> <chr>     <chr>   
      1  1500 s1        t1      
      2   280 s2        t1      
      3   456 s3        t1      
      4  1300 s1        t2      
      5   356 s2        t2      
      6   678 s3        t2      
      
      $samples
      # A tibble: 3 x 2
        sample  sample_id
        <chr>   <chr>    
      1 sample1 s1       
      2 sample2 s2       
      3 sample3 s3       
      
      $taxa
      # A tibble: 2 x 2
        taxon  taxon_id
        <chr>  <chr>   
      1 taxon1 t1      
      2 taxon2 t2      
      
      attr(,"class")
      [1] "tidytacos"

# Complain about missing packages

    Code
      force_optional_dependency("dadjokeapi")
    Condition
      Error in `force_optional_dependency()`:
      ! The dadjokeapi package must be installed to use this function. 


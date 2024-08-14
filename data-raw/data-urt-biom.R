library(tidytacos)

urt <- read_tidytacos("tidytacos/urt")
urt %>% to_biom("biom/urt.biom")
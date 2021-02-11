#SARS-CoV-2 research on Ghana
# 11th February, 2021 by OKOH Olayinka Sunday
#Anchor University, Lagos.


library(readr) #to read in tsv data
library(dplyr) #to manipulate data
#library(read)
#raw_data <- read_tsv(file.choose()) #This could have been used interactively. But to avoid manual selection each time
#read in the data
raw_data<- read_tsv("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/updatedSARCoV2_project/metadata_2021-01-06_18-26.tsv.zip", col_names = TRUE,
                    col_types = cols(strain ="f", virus = "c", gisaid_epi_isl = "-",
                                     genbank_accession = "-", date = "D", region = "f",
                                     country = "c",  division = "c", location = "c",
                                     region_exposure = "c", country_exposure = "c", division_exposure = "c",
                                     segment = "c", host = "f", age = "i", sex = "f",
                                     Nextstrain_clade = "f", pangolin_lineage = "f", GISAID_clade = "f",
                                     originating_lab = "c", submitting_lab = "c", authors = "-", url = "-", title = "-",
                                     paper_url = "-",date_submitted = "D" 
                    ))


#filter Ghana data
ghana_gisaid <- raw_data %>% filter(country == "Ghana") %>% select(-c(virus, region, country, region_exposure, country_exposure, division_exposure, segment, host, originating_lab, submitting_lab)) %>% as.data.frame()

#Exploratory analysis
dim(ghana_gisaid) # 70 sequence submitted
table(ghana_gisaid$sex) # 28 males, 11 females, 31 persons undeclared sex
sum(is.na(ghana_gisaid$age)) # 37 persons without age
range(na.omit(ghana_gisaid$age)) # 16 - 83

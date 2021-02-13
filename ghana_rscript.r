#SARS-CoV-2 research on Ghana
# 11th February, 2021 by OKOH Olayinka Sunday
#Anchor University, Lagos.


library(readr) #to read in tsv data
library(dplyr) #to manipulate data
library(ggplot2)

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
summary(ghana_gisaid$age) #Min is 16, 1st 29, Median 37, Mean 38.88, 3rd qd 49.00, Max 83.00, NAs is 37
mean_age <- mean(na.omit(ghana_gisaid$age)) #38.88
range(ghana_gisaid$length) #28,878 - 29,899 is the range of the sequence length
####
ggplot(ghana_gisaid, aes(x = age)) + geom_density(fill = "green") +
  geom_vline(xintercept = (mean_age)) 
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_agehist.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_agehist.pdf", height =  120, width = 200, units = "mm")

###

#get the number of sequences submitted from each state
ghana_sequences <- ghana_gisaid %>% select(division, strain) %>% group_by(division) %>% summarize (sequences = n()) %>% mutate(Percentage = round((sequences/sum(sequences)*100),1))

length(unique(ghana_gisaid$pangolin_lineage)) #13 lineages circulating in Ghana 
pango_lineages <- ghana_gisaid %>% group_by(pangolin_lineage) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1))
Nextstrain_clade <- ghana_gisaid %>% group_by(Nextstrain_clade) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1))
GISAID_clade <- ghana_gisaid %>% group_by(GISAID_clade) %>% summarize(N = n()) %>% mutate (Percentage = round((N/sum(N))*100,1))


#Number of Pango lineages in ghana
#Demography
#sex
sex <- c("Male", "Female", "NA")
patients <- as.data.frame(sex)
patients$N <- c(28, 11, 31)
patients <- patients %>% mutate(Percentage = round((N/sum(N))*100, 1))

ggplot(patients, aes(x = reorder(sex, N), y = N)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5)+
  labs(title = "Sex of SARS-CoV-2 patients in Ghana", y = "Number of patients", x = "Sex of patients")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_sex_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_sex_bar.pdf", height =  120, width = 200, units = "mm")

##Regions that submitted sequences
ggplot(ghana_sequences, aes(y = reorder(division, sequences), x = sequences)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "SARS-CoV-2 sequences submitted in GISAID from Ghana", y = "Region", x = "Sequences")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_sequences_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_sequences_bar.pdf", height =  120, width = 200, units = "mm")

##Pango lineages circulating in Ghana
ggplot(pango_lineages, aes(y = reorder(pangolin_lineage, N), x = N)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "SARS-CoV-2 lineages circulating in Ghana", y = "lineages", x = "Number of sequences", caption = "Source: GISAID, 2021")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_lineages_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_lineages_bar.pdf", height =  120, width = 200, units = "mm")

##NextStrain clades diversity 
ggplot(na.omit(Nextstrain_clade), aes(y = reorder(Nextstrain_clade, N), x = N)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "NextStrain clade diversity in Ghana", y = "NextStrain clades", x = "Number of sequences", caption = "Source: GISAID, 2021")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_Nextstrain_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_Nextstrain_bar.pdf", height =  120, width = 200, units = "mm")




##GISAID Clade diversity
ggplot(GISAID_clade, aes(y = reorder(GISAID_clade, N), x = N)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_dodge(width = 1), vjust = 0.5) +
  labs(title = "GISAID clade diversity in Ghana", y = "GISAID clades", x = "Number of sequences", caption = "Source: GISAID, 2021")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_GISAID_Clade_bar.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_GISAID_clade_bar.pdf", height =  120, width = 200, units = "mm")




###To plot lock down date
ghana_date <- ghana_gisaid %>% select(strain, date) %>% group_by(date) %>% summarize(N = n())

#Plot a chart 
ggplot(ghana_date, aes(x = date, y = cumsum(N))) +
  geom_line(color = 'red') +
  geom_vline(xintercept = as.Date(c("2020-03-30","2020-07-01"))) +
  labs(title = "Cumulative SARS-CoV-2 sequences in Ghana", y = "Number of sequences", x = "Date", caption = "Source: 2021")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_lockdown.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_lockdown.pdf", height =  120, width = 200, units = "mm")




##To plot incidence graph of the top 5 lineages

Top_lineages <- ghana_gisaid %>% select(pangolin_lineage, date)  %>% filter(pangolin_lineage == c("A", "B.1.1", "A.1.1", "B.1", "B.1.1.307")) %>% group_by(date, pangolin_lineage) %>% summarise (N = n()) %>% mutate(percentage = N/sum(N))

ggplot(Top_lineages, aes(x = date, y = percentage, fill = pangolin_lineage)) +
                          geom_area(size = 0.3, color = 'white') + 
  labs(title = "Incidence of top pango lineages circulating in Ghana", y = "Percentage", x = "Date", caption =  "Source: GISAID, 2021")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_incidence.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_incidence.pdf", height =  120, width = 200, units = "mm")

###Incidence of all the lineages in Ghana
lineages_incidence <- ghana_gisaid %>% select(pangolin_lineage, date)  %>% group_by(date, pangolin_lineage) %>% summarise (N = n()) %>% mutate(percentage = N/sum(N))
ggplot(lineages_incidence, aes(x = date, y = percentage, fill = pangolin_lineage)) +
  geom_area(size = 0.3, color = 'white') + 
  labs(title = "Incidence Pango lineages circulating in Ghana", y = "Percentage", x = "Date", caption =  "Source: GISAID, 2021")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_incidence_all.jpg", height =  120, width = 200, units = "mm")
ggsave("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/plots/ghana_incidence_all.pdf", height =  120, width = 200, units = "mm")






##### CASES ARE REQUIRED FOR DATA SCRIPT BELOW
#read in data from ncdc.gov.ng
#updated_ghana_ncdc <- read_excel("C:/Users/YinkaOkoh/Desktop/Bioinformatics_Data/SARS-CoV2_Project_2021/SARS-CoV-2_Ghana2020/updated_ghana_covid19_ncdc.xlsx")
#ghana_ncdc <- as.data.frame(updated_ghana_ncdc)
#ghana_cases_and_sequences <- left_join(ghana_ncdc, ghana_sequences, by = "State")

#Manipulate to include some calculated data
#ghana_cases_and_sequences <- ghana_cases_and_sequences %>% mutate(Cases_1M = Confirmed/Population, Test_1M = Testing/Population, Positive_1KTests = (Confirmed/Testing)*1000, Sequences_1HCases = (sequences/Confirmed)*100, Recovery_1HCases = (Recoveries/Confirmed)*100, Deaths_1HCases = (Deaths/Confirmed)*100)

#for GISAID data set = ghana_gisaid
#for ncdc data set and sequences = ghana_cases_and_sequences

#State without sequences in gisaid 
#sum(is.na(ghana_cases_and_sequences$sequences)) # 15























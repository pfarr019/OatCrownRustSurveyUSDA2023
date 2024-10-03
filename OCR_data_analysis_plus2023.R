
# R code for 2023 USDA Oat Crown Rust Virulence Survey --------------------


# set working directory
setwd("C:/Users/ERIN.MOREAU/OneDrive - USDA/Cloud documents/Oat research/OCR manuscript/APS 2024/Reanalysis with 2023 data")
# setwd("")

# in the original spreadsheet, I changed the 0-4 scale codes of 2015-2022 to the following scale:
# 0=0
# 0;=1
# ;=2
# ;C=3
# 1;=4
# 1=5
# 2=6
# 3=7
# 3+=8
# 4=9

#this conversion is based on miller et al 2021

library(tidyverse) #several packages including ggplot2 and dplyr

#read in data
raw_OCR_survey <- read.csv("Supplementary File S1 Raw data 1993-2023.csv")


# Transforming and cleaning the data --------------------------------------

# check data
str(raw_OCR_survey)
table(raw_OCR_survey$HiFi)
table(raw_OCR_survey$TAM.O.405)
table(raw_OCR_survey$H548)

# change the integer type to numeric
raw_OCR_survey$TAM.O.405 <- as.numeric(raw_OCR_survey$TAM.O.405)

# clean up the isolate names to the get rid of "-" between state and number for some of the years (95LA-026 -> 95LA026)
raw_OCR_survey <- raw_OCR_survey %>%
  mutate(isolate= gsub("([A-Za-z])-", "\\1", isolate)) #only deleting "-" that follow a letter

# Doing a race clone correction for the data BEFORE changing from more comprehensive scale 0-4 or 0-9 to binary avir (0) or vir (1)
# if two or more single pustule isolates from the same collection (usually submitted as Pca leaf or single plant) ex. 16TX31-1 and 16TX31-2, have the
# exact same virulence profile then that is an artificial oversampling of that race.

# First, make a secondary ID column that only has the initial ID before the "-"
raw_OCR_survey_1 <- raw_OCR_survey %>% separate(isolate, c("new_id", "second"), sep = "-", remove=FALSE) %>% select(-second)

str(raw_OCR_survey_1)

# dropping all duplicate rows based on new_id and the differential line scores
raw_OCR_survey_1 <- raw_OCR_survey_1 %>% distinct(new_id,Pc14,Pc35,Pc36,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,Pc91,Pc94,Pc96,MARVELOUS,H548,IAB605Xsel.,WIX4361.9,TAM.O.405,Belle,HiFi,Leggett,Stainless, .keep_all = TRUE)

# dropping the new_id column
raw_OCR_survey <- raw_OCR_survey_1[,-2]

# split dataframe into two based on the year, exclude 1992 because the data was collected by a different person and only includes accessions from MN, so not really part of the national survey
OCR_before2015 <- filter(raw_OCR_survey, year < 2015 & year >=1993)
OCR_after2015 <- filter(raw_OCR_survey, year >= 2015)

table(OCR_before2015$Pc14)
table(OCR_after2015$Pc14)

# for before 2015, take columns 8-47 with the rating data and change to 1 if virulent (>=3) or 0 if avirulent
OCR_before2015[8:47] <- lapply(OCR_before2015[8:47], function(x) ifelse(x>=3, 1, 0))

# for 2015 and after, take columns 8-47 with the rating data and change to 1 if virulent (>=7) or 0 if avirulent
OCR_after2015[8:47] <- lapply(OCR_after2015[8:47], function(x) ifelse(x>=7, 1, 0))

# recombine the two datasets
OCR_survey <- rbind(OCR_before2015,OCR_after2015)

str(OCR_survey)

OCR_survey %>% count(region)

# for each isolate, calculate the average amount of virulences and the virulence count, append to the end of the dataframe
OCR_survey$meanvirulence = rowMeans(OCR_survey[,c(8:47)], na.rm = TRUE)
OCR_survey$countvirulence = rowSums(OCR_survey[,c(8:47)], na.rm = TRUE)


#filter out 2013 data because only two samples and add in later
OCR_wout_2013<- OCR_survey %>% filter(year!=2013)

# these two observations must have typos for values Pc91, Pc94,Pc96, IAB605Xsel., WIX4361.9, TAM.O.405, Belle, HiFi, Leggett, and Stainless
# because they were not being tested for at that time.
# replace with NAs to clean up the data set
OCR_wout_2013 %>%
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="01TX010", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="01TX010", NA, .))) -> OCR_wout_2013_clean



# Summary measures ---------------------------------


# checking the isolates with the most and least virulences in 2023
OCR_wout_2013_clean %>% 
  filter(year==2023) %>% 
  arrange(countvirulence)

OCR_wout_2013_clean %>% 
  filter(year==2023) %>% 
  arrange(desc(countvirulence))

# determining which isolates are virulent to Pc14, Pc40, Pc91, and Pc94
OCR_wout_2013_clean %>% 
  filter(year==2023) %>% 
  filter(Pc14==1) %>% 
  filter(Pc40==1) %>% 
  filter(Pc91==1) %>% 
  filter(Pc94==1)

# determining how the mean and median countvirulence compares to previous survey years
OCR_wout_2013_clean %>%
  group_by(year) %>% 
  summarise(mean(countvirulence)) %>% 
  arrange(desc(year))

OCR_wout_2013_clean %>%
  group_by(year, region) %>% 
  summarise(mean(countvirulence)) %>% 
  arrange(desc(year))

  

# Pca race analysis -------------------------------------------------------


## for all study years only using the 30 original lines and excluding the isolates with missing data --------

library(data.table)

OCR_strains_allminusNAs <- OCR_wout_2013_clean %>%
  select(-c(Belle, Stainless,Leggett,HiFi,IAB605Xsel.,WIX4361.9,TAM.O.405,Pc91,Pc94,Pc96)) %>% 
  na.omit()
setDT(OCR_strains_allminusNAs)[,group :=.GRP,by = .(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548,MARVELOUS)]

OCR_n_strains_allminusNAs <-OCR_strains_allminusNAs %>% 
  count(group) %>% 
  arrange(-n)

strain_counts_all <-OCR_n_strains_allminusNAs %>% 
  count(n, name= "frequency") %>% 
  mutate(totalNofstrains= n*frequency) %>%
  mutate(percentage= totalNofstrains/sum(totalNofstrains)) %>% 
  mutate(percentageofstrains= frequency/sum(frequency))
sum(strain_counts_all$frequency)
sum(strain_counts_all$totalNofstrains)

isolate1688_all <-OCR_strains_allminusNAs %>% 
  filter(group==1688)
isolate1688_all %>% count(region, year)

isolate2443_all <-OCR_strains_allminusNAs %>% 
  filter(group==2443)
isolate2443_all %>% count(region, year)

isolate561_all <-OCR_strains_allminusNAs %>% 
  filter(group==561)
isolate561_all %>% count(region, year)

isolate1676_all <-OCR_strains_allminusNAs %>% 
  filter(group==1676)
isolate1676_all %>% count(region, year)

isolate26_all <-OCR_strains_allminusNAs %>% 
  filter(group==26)
isolate26_all %>% count(region, year)

# Sorting out and printing the isolates with the least amount of virulences
head(OCR_strains_allminusNAs[order(OCR_strains_allminusNAs$countvirulence),], n=20L)

# determining how many races are exclusive to each region when there are at least 5 isolates
# sort out any races with less than 5 isolates
OCR_5isolates_plus <- OCR_n_strains_allminusNAs[OCR_n_strains_allminusNAs$n >= 5]
OCR_5isolates_plus_data <-dplyr::filter(OCR_strains_allminusNAs, group %in% OCR_5isolates_plus$group)
str(OCR_5isolates_plus_data)
OCR_5isolates_plus_data %>% count(group,region)
by_year <-OCR_5isolates_plus_data %>% count(group,year)


# visualize with a bar chart (sqrt scale) with numbers
ggplot(strain_counts_all, aes(x=n, y=frequency)) + 
  geom_bar(stat='identity', fill="light blue")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=frequency), vjust=-0.3, size=3)+
  # scale_y_sqrt()+
  coord_trans(y= 'sqrt')+
  scale_y_continuous(limits=c(0,3000), expand=c(0,0))+
  scale_x_continuous(limits=c(-0.5,47), expand=c(0,0))+
  xlab("Number of Isolates per Race")+
  ylab("Number of Races \n(Square root axis transformation)")+
  ggtitle("Oat Crown Rust Races")

ggsave("Histogram races all years.tiff", scale=1, dpi=600, width=150, height=100, units="mm")


## 2023 Pca race analysis --------------------------------------------------

# for 2023 only, using all 40 lines, excluding missing data

OCR_strains_2023minusNAs <- OCR_survey %>%
  filter(year == 2023) %>%
  na.omit()
setDT(OCR_strains_2023minusNAs)[,group :=.GRP,by = .(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548,MARVELOUS, Belle, Stainless,Leggett,HiFi,IAB605Xsel.,WIX4361.9,TAM.O.405,Pc91,Pc94,Pc96)]

str(OCR_strains_2023minusNAs)
tail(OCR_strains_2023minusNAs)

OCR_n_strains_2023minusNAs <-OCR_strains_2023minusNAs %>% 
  count(group) %>% 
  arrange(-n)

strain_counts_2023 <-OCR_n_strains_2023minusNAs %>% 
  count(n, name= "frequency") %>% 
  mutate(totalNofstrains= n*frequency) %>%
  mutate(percentage= totalNofstrains/sum(totalNofstrains)) %>% 
  mutate(percentageofstrains= frequency/sum(frequency))
sum(strain_counts_2023$frequency)
sum(strain_counts_2023$totalNofstrains)
strain_counts_2023

isolate2_all <-OCR_strains_2023minusNAs %>% 
  filter(group==2)
isolate2_all %>% count(region, year)

isolate1_all <-OCR_strains_2023minusNAs %>% 
  filter(group==1)
isolate1_all %>% count(region, year)

isolate116_all <-OCR_strains_2023minusNAs %>% 
  filter(group==116)
isolate116_all %>% count(region, year)

isolate123_all <-OCR_strains_2023minusNAs %>% 
  filter(group==123)
isolate123_all %>% count(region, year)

isolate31_all <-OCR_strains_2023minusNAs %>% 
  filter(group==31)
isolate31_all %>% count(region, year)

# Plotting the heat map of observations per state -------------------------

## Only 2023 ---------------------------------------------------------------

states_us <- ne_states(country="United States of America",returnclass = 'sf')
class(states_us)

#working from the original dataframe, need to rename state column to postal so that it merges correctly with the sf object 
str(OCR_survey)
OCR_survey_2023 <- OCR_survey %>% 
  filter(year == 2023)
str(OCR_survey_2023)
samples_per_state <- OCR_survey_2023 %>% count(state, sort = FALSE, name="rust_count") %>% rename(postal="state")
#add region to the samples_per_state by matching columns from OCR survey
# samples_per_state <-left_join(samples_per_state, OCR_survey, by= "postal")
# merge(samples_per_state, OCR_survey[, c("state", "region")], by="state")
# %>% rename(postal="state")

#merge the sf and the samples_per_state data frame
states_us <- merge(states_us, samples_per_state, all = TRUE)

#read in a simple CSV with the regions N vs S defined and merge with states_us
regions_key <- read.csv("Postal codes and regions.csv")
states_us <-merge(states_us, regions_key, by="postal")

#statistics on the percentages for each state for the results section
sum(samples_per_state$rust_count)
samples_per_state$perc_per_state =(samples_per_state$rust_count/sum(samples_per_state$rust_count)*100)
arrange(samples_per_state, rust_count)
# Stats on North vs South 
str(OCR_survey_2023)
nrow(OCR_survey_2023[OCR_survey_2023$region == 'South', ])
nrow(OCR_survey_2023[OCR_survey_2023$region == 'North', ])

(mainland <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt", name="Sample count")+
    geom_sf_label(aes(label = rust_count), size = 2, fontface = "bold", label.padding = unit(0.15, "lines")) +
    geom_sf(fill = "transparent", color = "black", linewidth=0.3, #adding a thick border around the two regions
            data = . %>% group_by(rust_region) %>% summarise()) +
    coord_sf(crs = st_crs(2163), xlim = c(-2500000, 2500000), ylim = c(-2300000, 730000)))

(alaska <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt")+
    theme(legend.position = "none")+
    geom_sf_text(aes(label = rust_count), size = 2, fontface = "bold") +
    geom_sf_label(aes(label = rust_count), size = 2, fontface = "bold", label.padding = unit(0.15, "lines")) +
    coord_sf(crs = st_crs(3467), xlim = c(-2400000, 1600000), ylim = c(200000, 2500000), expand = FALSE, datum = NA))+
  xlab("") + ylab("")

(hawaii  <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt")+
    theme(legend.position = "none")+
    coord_sf(crs = st_crs(4135), xlim = c(-161, -154), ylim = c(18, 
                                                                23), expand = FALSE, datum = NA))

mainland + xlab("Longitude") + ylab("Latitude") +
  ggtitle("2023 Rust survey samples by state")+
  annotation_custom(
    grob = ggplotGrob(alaska + theme(axis.title.x=element_blank(), axis.title.y=element_blank())),
    xmin = -2750000,
    xmax = -2750000 + (1600000 - (-2400000))/2.5,
    ymin = -2450000,
    ymax = -2450000 + (2500000 - 200000)/2.5
  ) +
  annotation_custom(
    grob = ggplotGrob(hawaii),
    xmin = -1250000,
    xmax = -1250000 + (-154 - (-161))*120000,
    ymin = -2450000,
    ymax = -2450000 + (23 - 18)*120000
  )

ggsave("Figure 1 Rust sample map 2023.tiff", scale=1, dpi=600, width=178, units="mm")



library("rnaturalearth")
library("sf")

## All years, from 1993-2023 -----------------------------------------------

states_us <- ne_states(country="United States of America",returnclass = 'sf')
class(states_us)

#working from the original dataframe, need to rename state column to postal so that it merges correctly with the sf object 
str(OCR_survey)
samples_per_state <- OCR_survey %>% count(state, sort = FALSE, name="rust_count") %>% rename(postal="state")
#add region to the samples_per_state by matching columns from OCR survey
# samples_per_state <-left_join(samples_per_state, OCR_survey, by= "postal")
# merge(samples_per_state, OCR_survey[, c("state", "region")], by="state")
# %>% rename(postal="state")

#merge the sf and the samples_per_state data frame
states_us <- merge(states_us, samples_per_state, all = TRUE)

#read in a simple CSV with the regions N vs S defined and merge with states_us
regions_key <- read.csv("Postal codes and regions.csv")
states_us <-merge(states_us, regions_key, by="postal")

#statistics on the percentages for each state for the results section
sum(samples_per_state$rust_count)
samples_per_state$perc_per_state =(samples_per_state$rust_count/sum(samples_per_state$rust_count)*100)
arrange(samples_per_state, rust_count)
# Stats on North vs South 
str(OCR_survey)
nrow(OCR_survey[OCR_survey$region == 'South', ])
nrow(OCR_survey[OCR_survey$region == 'North', ])

(mainland <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt", name="Sample count")+
    geom_sf_label(aes(label = rust_count), size = 2, fontface = "bold", label.padding = unit(0.15, "lines")) +
    geom_sf(fill = "transparent", color = "black", linewidth=0.3, #adding a thick border around the two regions
            data = . %>% group_by(rust_region) %>% summarise()) +
    coord_sf(crs = st_crs(2163), xlim = c(-2500000, 2500000), ylim = c(-2300000, 730000)))

(alaska <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt")+
    theme(legend.position = "none")+
    geom_sf_text(aes(label = rust_count), size = 2, fontface = "bold") +
    geom_sf_label(aes(label = rust_count), size = 2, fontface = "bold", label.padding = unit(0.15, "lines")) +
    coord_sf(crs = st_crs(3467), xlim = c(-2400000, 1600000), ylim = c(200000, 2500000), expand = FALSE, datum = NA))+
  xlab("") + ylab("")

(hawaii  <- ggplot(data = states_us) +
    geom_sf(data=states_us, aes(fill = rust_count)) +
    scale_fill_gradient(low="light blue", high="dark blue", trans= "sqrt")+
    theme(legend.position = "none")+
    coord_sf(crs = st_crs(4135), xlim = c(-161, -154), ylim = c(18, 
                                                                23), expand = FALSE, datum = NA))

mainland + xlab("Longitude") + ylab("Latitude") +
  ggtitle("1993-2023 Rust survey samples by state")+
  annotation_custom(
    grob = ggplotGrob(alaska + theme(axis.title.x=element_blank(), axis.title.y=element_blank())),
    xmin = -2750000,
    xmax = -2750000 + (1600000 - (-2400000))/2.5,
    ymin = -2450000,
    ymax = -2450000 + (2500000 - 200000)/2.5
  ) +
  annotation_custom(
    grob = ggplotGrob(hawaii),
    xmin = -1250000,
    xmax = -1250000 + (-154 - (-161))*120000,
    ymin = -2450000,
    ymax = -2450000 + (23 - 18)*120000
  )

ggsave("Rust sample map 1993-2023.tiff", scale=1, dpi=600, width=178, units="mm")


# Boxplot count virulence -------------------------------------------------


# create boxplot showing countvirulence by year with the alpha (transparency) corresponding to the relative amount of samples

# counting up the number of observations per year and dividing by the largest to get values between 0 and 1
samples_per_year <- OCR_wout_2013_clean %>% count(year, sort = FALSE) %>% mutate(n_adjusted=n/max(n))

typeof(samples_per_year)
str(samples_per_year)


boxplot1 <-ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2023), y=countvirulence)) + 
  geom_boxplot(fill="dark orange", alpha=samples_per_year$n_adjusted, outlier.alpha = 1, outlier.size=0.5, notch=TRUE)+
  xlab("Year")+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2023), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 12.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 12.5, xmax = 17.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot 
  theme(text=element_text(size=10))
  
boxplot1

# manually print out a legend color scale for the side because unable to print out a legend using alpha

library(ggpubr)

boxplot_legend <-ggplot(samples_per_year, aes(x=factor(year, levels=1993:2023), y= n_adjusted, fill=n)) + 
  geom_point()+
  scale_fill_gradient(low=rgb(1,0.549,0,0), high=rgb(1,0.549,0,1), limit = c(0,544), name="Number of \nIsolates",
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+ #manually specifying dark orange with alpha 1 and alpha 0
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8.5))

boxplot_legend
boxplot_legend <- get_legend(boxplot_legend) #extracting the legend
boxplot_legend <- as_ggplot(boxplot_legend) #making the legend a ggplot object
boxplot_legend

boxplotwlegend <- ggarrange(boxplot1, boxplot_legend, widths=c(1,0.12),
                          nrow=1)
boxplotwlegend
ggsave("Boxplot all years with legend.tiff", scale=1, dpi=600, width=150, height=113, units="mm", path=".", bg="white")



# 2023 data heatmap -------------------------------------------------------

# starting with the unconverted 0-9 data
# converting 0-9 scale to 0-1 scale with three levels, for visualization
# 0-3=0
# 4-6=0.2 (intermediate resistance, 1 or greater (some sporulation))
# 7-9=1

# The intermediate value is less than 0.5 so that the clustering algorithm puts more weight into S virulence rating when clustering
# Take columns 8-47 with the rating data and change to above rating system
str(raw_OCR_survey)
OCR_2023_heatmap <- raw_OCR_survey %>% 
  filter(year == 2023) %>% 
  mutate(across(c(8:47),
                ~case_when(. <4 ~ 0,
                           . ==4 ~ 0.2,
                           . ==5 ~ 0.2,
                           . ==6 ~ 0.2,
                           . >6 ~ 1))) %>% 
  select(!(MARVELOUS))
str(OCR_2023_heatmap)

#load required packages
library(ComplexHeatmap)
library(circlize)
#library(RColorBrewer) #another alternative package for colors if so desired
library(colorRamps)

#Reassigning the isolate variable (first column) to row names
OCR_2023_heatmap_named_rows <- OCR_2023_heatmap[,-1]
rownames(OCR_2023_heatmap_named_rows) <- OCR_2023_heatmap[,1]

# removing supplementary data columns
OCR_2023_heatmap_matrix <- as.matrix(OCR_2023_heatmap_named_rows[  ,c(7:ncol(OCR_2023_heatmap_named_rows))])
str(OCR_2023_heatmap_matrix)

# preparing data for percentage virulent row annotations for each differential line
str(OCR_wout_2013_clean)
heatmap_data_2023 <- OCR_wout_2013_clean %>% group_by(year) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) %>% 
  filter(year==2023) %>% 
  select(!(MARVELOUS))

# change data into long format
heatmap_data_long_2023 <-gather(heatmap_data_2023, differential_line, percent_virulent, Pc14:Stainless) %>% 
  select(-year)
str(heatmap_data_long_2023)
#converting into a matrix
matrix_percent <- heatmap_data_long_2023 %>% 
  select(-differential_line) %>% 
  as.matrix()

# manually assigning a vector of colors based on the percentage of virulence for each accession using colorRampPalette to simulate the YlOrRd color palette
# Load the RColorBrewer package to check the original YlOrRd palette (optional) 
# install.packages("RColorBrewer") 
library("RColorBrewer")
# View the original YlOrRd colors (optional) 
display.brewer.pal(9, "YlOrRd") 
# Display the 9-color version for reference 
# Manually specify colors similar to the YlOrRd palette 
ylorrd_colors <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026") 
# Create a color ramp function from these colors 
ylorrd_ramp <- colorRampPalette(ylorrd_colors) 
# Generate 100 colors from this palette 
custom_palette <- ylorrd_ramp(101) 
# Plot to visualize the recreated palette (optional) 
image(1:100, 1, as.matrix(1:100), col = custom_palette, xlab = "YlOrRd", ylab = "", yaxt = "n", bty = "n")
palette_func <-colorRampPalette(c("white", "yellow","red")) 
colors=palette_func(101)
custom_palette
# creating a vector of colors based on the percentage values of virulence to each differential line
color_vector <- custom_palette[as.numeric(cut(matrix_percent,breaks=101, labels=FALSE))]
color_vector

# for a simpler palette, but not as good at showing differences over a wide range of values
# palette_func <-colorRampPalette(c("white", "yellow","red")) 
# colors=palette_func(101)
# color_vector <- colors[as.numeric(cut(matrix_percent,breaks=101, labels=FALSE))] #could also use
# color_vector

# Manually adding a zero to the beginning of the one digit numbers so they stack up correctly
matrix_percent<-round(matrix_percent, digits=0)
matrix_percent <- matrix(sprintf("%02d", matrix_percent), nrow(matrix_percent))
# adding leading and trailing spaces to characters so that they have space in the 
# the just parameter does not work correctly so one has to do this manual work around 
matrix_percent <- matrix(paste(" ", matrix_percent, " ", sep = ""), nrow = nrow(matrix_percent)) 

# setting up the row annotation
row_ha = rowAnnotation("Percent\nVirulent"= anno_text(matrix_percent, gp=gpar(fill=color_vector, fontsize=10, just="right"), show_name = TRUE))


# preparing data for the column annotation N vs S
Region1 <- select(OCR_2023_heatmap, region)
Region1 <- unlist(Region1)
Region1 <- unname(Region1)
typeof(Region1)
is.character(Region1) #should return TRUE
tOCR_2023_heatmap_matrix <- t(OCR_2023_heatmap_matrix) #transpose the matrix
#adopting colors from ColorBrewer 2.0 for similarity to purple=north and green=south in the "Difference in percentage of virulent isolates per year by region" graph, although I made the purple is a "rung" darker than the green on the scale here (https://colorbrewer2.org/#type=diverging&scheme=PRGn&n=8) so they are easier to distinguish
column_ha = HeatmapAnnotation(Region=Region1, col=list(Region= c("North" ="#762a83", "South" = "#5aae61"))) 

# making a custom color palette for the heatmap, using a 5 color scale without the 3 and 4 intermediate values to emphasize the virulent isolates while still showing a difference between the HR and MR
yellowred = colorRampPalette(c("lightyellow","red3"), space="rgb")(4)
yellowred =yellowred[-c(2)]

# using Heatmap from the ComplexHeatmap package
# Different clustering methods with hclust
# The default method for clustering_distance_rows is "euclidean" and the default method for clustering_method_rows is "complete".
scoring_plot <- Heatmap(tOCR_2023_heatmap_matrix,
                        show_row_names = TRUE,
                        show_column_names = FALSE,
                        row_names_side = "left",
                        border_gp = gpar(col="black"),
                        column_names_gp = gpar(cex=0.8),
                        column_title = "Isolates",
                        row_title = "Differential Lines",
                        top_annotation = column_ha,
                        left_annotation= row_ha,
                        heatmap_legend_param = list(title = "Differential\nResponse", color_bar = "discrete", at = c(0, 0.2, 1), labels = c("HR", "MR", "S")), col = yellowred)

scoring_plot


# Save the plots to a file
png("Heatmap 2023 isolates.tiff", width = 3500, height = 2800, units="px", res=300)
#pdf("scoring_plot_.pdf")
scoring_plot
dev.off()

# Violin Plot N vs S ------------------------------------------------------


# making a violin plot to compare N and S
library(cowplot)
OCR_2023 <- filter(OCR_wout_2013_clean, year == 2023)
# Perform a wilcox test
compare_means(countvirulence ~ region, data = OCR_2023)
my_comparisons <- list(c("North", "South"))
violin_plot <- ggplot(OCR_2023, aes(x=region, y=countvirulence, fill=region)) + 
  geom_violin(trim=TRUE)+
  theme_minimal_hgrid()+
  ylim(0,40)+
  coord_cartesian(expand=FALSE)+
  scale_fill_manual(values=c("#762a83", "#5aae61"))+
  labs(title="Distribution of isolate virulence by region", y="Number of virulences per isolate", x="Region")+
  theme(legend.position="none")+
  geom_boxplot(width=0.1, fill="white", color="light gray", outlier.size=3)+
  stat_compare_means( aes(label = paste0(..method.., "\n", "p = ", ..p.format..)), label.x = 1.35, label.y = 35)
violin_plot

ggsave("N vs S violin plot.tiff", scale=1, units="mm")


# Heatmap of virulences per differential line over time -------------------

library("colorspace")

str(OCR_wout_2013_clean)
heatmap_data <- OCR_wout_2013_clean %>% group_by(year) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100))

# change data into long format
heatmap_data_long <-gather(heatmap_data, differential_line, percent_virulent, Pc14:Stainless)
str(heatmap_data_long)

# createheatmap
heatmap_year <-ggplot(heatmap_data_long, aes(factor(year), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 3)+
  ylab("Differential Line")+
  ggtitle("Percentage of virulent isolates per year")
heatmap_year

# Order the heatmap rows by hierarchical clustering so that patterns are easier to see

heatmap_data_matrix <- as.matrix(heatmap_data[, -1]) # -1 to omit categories from matrix
# Cluster based on euclidean distance
clust <- hclust(dist(t(heatmap_data_matrix)))
rust_clust <- as.dendrogram(hclust(dist(t(heatmap_data_matrix))))

library("ggdendro")

# Create dendrogram
dendro_data_for_plot <- dendro_data(rust_clust, type = "rectangle")
dendro_rust <- ggplot(segment(dendro_data_for_plot)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() +
  theme_dendro()+
  theme(plot.margin = unit(c(5.5,0,5.5,5.5), "points"), text=element_text(size=7.5))+
  scale_y_reverse(expand = c(0, 0))+
  scale_x_continuous(expand=c(0.015, 0.015))+
  annotate("text", x=29.0, y=170, label= "A", color="red", cex=3)+ #manually designate the groups of isolates that I discuss in the article
  annotate("text", x=13.3, y=230, label= "B", color="red", cex=3)+
  annotate("text", x=4.7, y=230, label= "C", color="red", cex=3)+
  annotate("text", x=24.2, y=169, label= "D", color="red", cex=3)

dendro_rust

#ordered heatmap using dendrogram
heatmap_year_ordered <-ggplot(heatmap_data_long, aes(factor(year), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, size=11), plot.margin = unit(c(5.5,5.5,5.5,0), "points"), text=element_text(size=9))+ #decrease plot margin on left side so lines up closer to the dendrogram
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2.5)+
  # ylab("Differential Line")+
  ggtitle("Percentage of virulent isolates per year \n")+
  annotate("segment", x = 20.5, xend = 20.5, y = 40.5, yend = 0.5, colour = "gray", linewidth=1, alpha=0.4)
heatmap_year_ordered

# create summary figure grouped by decade
# 1993-2002, 2003-2012, 2013-2022
# Use the clean data with the two observations from 2013
OCR_survey %>%
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="96WI064", NA, .))) %>% 
  mutate(across("Pc91":"Pc96", ~ifelse(isolate=="01TX010", NA, .))) %>% 
  mutate(across("IAB605Xsel.":"TAM.O.405", ~ifelse(isolate=="01TX010", NA, .))) -> OCR_survey_clean

# now subset data into three decades and take the mean
OCR_survey_clean %>%
  filter(year<=2002) %>% 
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_1993_to_2002
OCR_survey_clean %>%
  filter(year>2002, year<2013) %>% 
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_2003_to_2012
OCR_survey_clean %>%
  filter(year>=2013) %>% 
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_2013_to_2023
# bind the means together
Mean_by_decade <- bind_rows("1993-2002"= Mean_1993_to_2002, "2003-2012"= Mean_2003_to_2012, "2013-2022"= Mean_2013_to_2023, .id="Decade")

# change data into long format
Mean_by_decade_long <-gather(Mean_by_decade, differential_line, percent_virulent, Pc14:Stainless)
str(Mean_by_decade_long)

# ordered heatmap by decade
heatmap_decade_ordered <- ggplot(Mean_by_decade_long, aes(factor(Decade), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  # scale_fill_gradient2(low = "light blue", mid = "white", high = "red", midpoint = 0.5)+
  # scale_fill_gradient2(low = "white", high = "red")+
  # scale_fill_distiller(palette = "YlGnBu")+
  # scale_fill_viridis_c(option="I")+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), text=element_text(size=7.5))+
  scale_x_discrete(labels=c("1993-2002" = "1993-\n2002", "2003-2012" = "2003-\n2012",
                            "2013-2023" = "2013-\n2023"))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2)+
  ggtitle("Percentage of \nvirulent isolates \nper decade")
heatmap_decade_ordered

#arrange the two ordered plots together with the dendrogram
arrangedplotsordered <- ggarrange(dendro_rust, heatmap_year_ordered, heatmap_decade_ordered,
                                  labels=c("a", "","b"),
                                  nrow=1,
                                  align="h",
                                  widths=c(0.08, 1, 0.24))
arrangedplotsordered

annotate_figure(arrangedplotsordered,
                left = text_grob("Differential Line", rot = 90))

#just the year plot with the dendrogram
arrangedplotsordered <- ggarrange(dendro_rust, heatmap_year_ordered,
                                  nrow=1,
                                  align="h",
                                  widths=c(0.08, .95))
arrangedplotsordered

annotate_figure(arrangedplotsordered,
                left = text_grob("Differential Line", rot = 90))


## Linear regression for each pc gene over time ----------------------------
# https://www.datacamp.com/tutorial/linear-regression-R
# lm([target] ~ [predictor / features], data = [data source])

# doing all the lines together actually ignores the years that don't have observations for all of the lines (only takes last 10 years)
# do separately for the three sets of lines depending on how long they have been used

lmOG =lm(cbind(Pc14,Pc35,Pc36,Pc38,Pc39,Pc40,Pc45,Pc46,Pc48,Pc50,Pc51,Pc52,Pc53,Pc54,Pc55,Pc56,Pc57,Pc58,Pc59,Pc60,Pc61,Pc62,Pc63,Pc64,Pc67,Pc68,Pc70,Pc71,H548, MARVELOUS)~year, data=heatmap_data)
summary(lmOG)

lm2005 =lm(cbind(Pc91,Pc94,Pc96,IAB605Xsel.,WIX4361.9,TAM.O.405)~year, data=heatmap_data)
summary(lm2005)

lm2010 =lm(cbind(Belle,HiFi,Leggett,Stainless)~year, data=heatmap_data)
summary(lm2010)

# Was not able to efficiently copy out the slope, p value, and r2 with code so did manually, results are in lmdata.csv

lmdata <- read.csv("lmdata.csv")
lmdata <- lmdata %>% 
  mutate(slope= case_when(Pvalue <= 0.05 ~ year)) %>%
  mutate(signif= case_when(Pvalue <=0.05 & Pvalue > 0.01 ~ "*",
                           Pvalue <=0.01 & Pvalue > 0.001 ~ "**",
                           Pvalue <=0.001 ~ "***")) %>%
  mutate(slope_sig=paste(as.character(round(slope,1)),signif, sep=" "))%>% #combining the slope and significance into one string
  mutate(slope_sig = na_if(slope_sig,"NA NA")) %>% #turning "Na Na" into a real <NA>
  mutate(y=1) #adding a y=1 for plotting in the heatmap

### creating an ordered heatmap to add to the final heatmap to communicate the slope and significance
heatmap_slope <-ggplot(lmdata, aes(factor(y), factor(gene, levels=rev(unique(gene))), fill=year))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  # scale_fill_gradient2(low = "light blue", mid = "white", high = "red", midpoint = 0.5)+
  # scale_fill_gradient2(low = "white", high = "red")+
  # scale_fill_distiller(palette = "YlGnBu")+
  # scale_fill_viridis_c(option="I")+
  scale_fill_continuous_diverging(palette = "Blue-Red 3", rev = FALSE, na.value="white")+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        plot.title = element_blank(), 
        plot.margin = unit(c(5.5,5.5,5.5,0), "points"), #decrease plot margin on left side so lines up closer to the dendrogram
        text=element_text(size=9))+ 
  scale_x_discrete(labels =c("1" = "change\nper year"))+
  geom_text(aes(label = slope_sig), color = "black", size = 2.5)+
  # ylab("Differential Line")+
  ggtitle("Slope")
heatmap_slope

#arrange all the plots together
arrangedplotsordered_slope <- ggarrange(dendro_rust, heatmap_year_ordered, heatmap_slope, heatmap_decade_ordered,
                                        labels=c("a", "","","b"),
                                        nrow=1,
                                        align="h",
                                        widths=c(0.08, 1, 0.06, 0.22))
arrangedplotsordered_slope

annotate_figure(arrangedplotsordered_slope,
                left = text_grob("Differential Line", rot = 90, size=8))


ggsave("Figure 4 heatmap 450dpi.tiff", scale=1, dpi=450, width=246, height=178, units="mm", path=".", bg="white")


#all plots minus decade
arrangedplotsordered_slope <- ggarrange(dendro_rust, heatmap_year_ordered, heatmap_slope,
                                        nrow=1,
                                        align="h",
                                        widths=c(0.08, 1, 0.06))
arrangedplotsordered_slope

annotate_figure(arrangedplotsordered_slope,
                left = text_grob("Differential Line", rot = 90, size=13))


ggsave("Supplementary Figure heatmap 30 years.tiff", scale=1, dpi=450, width=246, height=178, units="mm", path=".", bg="white")



# Comparing the S to N accessions wtih a Heatmap --------------------------

str(OCR_wout_2013_clean)
heatmap_data_byregion <- OCR_wout_2013_clean %>% group_by(region, year) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100))
str(heatmap_data_byregion)
head(heatmap_data_byregion, n=20L)
tail(heatmap_data_byregion, n=20L)

# Taking a closer look at the differences in numbers between northern and southern isolates in the same year to see if I can compare percentages
samples_per_year_perRegion <- OCR_wout_2013_clean %>% count(year, region, sort = FALSE)
samples_per_year_perRegion


ggplot(samples_per_year_perRegion, aes(x=year, y=n, fill=region)) +
  geom_bar(stat='identity', position='dodge')

heatmap_data_byregion
#arrange by year
heatmap_data_byregion <- arrange(heatmap_data_byregion, year)


# subtract the north by the south for each year by creating northern and southern data tibbles and extracting them from one another
heatmap_data_north <- heatmap_data_byregion %>%
  filter(region=="North")

heatmap_data_south <- heatmap_data_byregion %>%
  filter(region=="South")

#subtract north from south
heatmap_data_diff <-(heatmap_data_north[,c(3:42)])-(heatmap_data_south[,c(3:42)])

#re-add the year to the data set for visualization
year_vec <- pull(heatmap_data_north, year)
heatmap_data_diff<- heatmap_data_diff %>%
  mutate(year=year_vec,
         .before=Pc14)

### change data into long format
heatmap_data_diff_long <-gather(heatmap_data_diff, differential_line, percent_virulent, Pc14:Stainless)

#heatmap
heatmap_NvsS <- ggplot(heatmap_data_diff_long, aes(factor(year), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile()+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  xlab("Year")+
  labs(fill = "North Vir Perc minus \nSouth Vir Perc   \n")+
  # ylab("Differences in virulence percentages from North to South")+
  scale_fill_continuous_diverging(palette = "Purple_Green", rev = TRUE, limits=c(-100,100))+
  theme(plot.title = element_text(hjust = 0.5, size=11), axis.title.y=element_blank(), axis.title.x=element_blank(), legend.position="bottom", text=element_text(size=9))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2.5)+
  ggtitle("Difference in percentage of virulent isolates per year by region")+
  annotate("segment", x = 20.5, xend = 20.5, y = 40.5, yend = 0.5, colour = "gray", linewidth=1, alpha=0.6)
heatmap_NvsS


#######Decade heatmap comparing N vs S########

OCR_survey_clean %>%
  filter(region=="North") %>%
  filter(year<=2002) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_North_1993_to_2002
OCR_survey_clean %>%
  filter(region=="South") %>%
  filter(year<=2002) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_South_1993_to_2002

OCR_survey_clean %>%
  filter(region=="North") %>%
  filter(year>2002, year<2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_North_2003_to_2012
OCR_survey_clean %>%
  filter(region=="South") %>%
  filter(year>2002, year<2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_South_2003_to_2012

OCR_survey_clean %>%
  filter(region=="North") %>%
  filter(year>=2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_North_2013_to_2023
OCR_survey_clean %>%
  filter(region=="South") %>%
  filter(year>=2013) %>%
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100)) ->Mean_South_2013_to_2023

# bind the means together
South_decade_mean <- bind_rows("1993-2002"= Mean_South_1993_to_2002, "2003-2012"= Mean_South_2003_to_2012, "2013-2023"= Mean_South_2013_to_2023, .id="Decade")
North_decade_mean <- bind_rows("1993-2002"= Mean_North_1993_to_2002, "2003-2012"= Mean_North_2003_to_2012, "2013-2023"= Mean_North_2013_to_2023, .id="Decade")
str(North_decade_mean)

# subtract mean of south from mean of north
heatmap_data_mean_diff <- North_decade_mean[,c(2:41)] - South_decade_mean[,c(2:41)]
str(heatmap_data_mean_diff)

# Re-add the decade names to dataset
decade_vec <- pull(South_decade_mean, Decade)
heatmap_data_mean_diff<- heatmap_data_mean_diff %>%
  mutate(decade=decade_vec,
         .before=Pc14)
str(heatmap_data_mean_diff)

#put into long form
heatmap_data_mean_long <-gather(heatmap_data_mean_diff, differential_line, percent_virulent, Pc14:Stainless)
str(heatmap_data_mean_long)

#ordered heatmap
heatmap_mean_diff <- ggplot(heatmap_data_mean_long, aes(factor(decade), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  # labs(fill = "North Vir Perc\n      minus \nSouth Vir Perc")+
  scale_fill_continuous_diverging(palette = "Purple_Green", rev = TRUE, limits=c(-100,100))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5), text=element_text(size=8))+
  scale_x_discrete(labels=c("1993-2002" = "1993-\n2002", "2003-2012" = "2003-\n2012",
                            "2013-2023" = "2013-\n2023"))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 2.2)+
  ggtitle("Differences \nby decade")
heatmap_mean_diff

# Manually creating dendrogram groups as a graphic
line_coord <- data.frame(x1 = c(0,0,0,0,0,0,0,0,1,1,1,1), x2 = c(1,1,1,1,1,1,1,1,1,1,1,1), 
                         y1 = c(1,8,9,20,22,25,26,37,1,9,22,26), y2 = c(1,8,9,20,22,25,26,37,8,20,25,37))

dendro_rust_groups <- ggplot(line_coord, aes(x = x1, y = y1, xend = x2, yend = y2)) + 
  geom_segment() +
  # coord_flip() +
  theme_dendro()+
  theme(plot.margin = unit(c(5.5,0,5.5,5.5), "points"))+
  scale_y_continuous(expand = c(0, 0), limits=c(0.5,40.5))+ #forsome reason have to manually adjust the y axis limits so that it lines up properly on the combined graphic 
  scale_x_reverse(expand=c(0.015, 0.015), limits=c(3.5,0))+
  annotate("text", x=2.5, y=(37-26)/2+26, label= "A", color="red", cex=3)+
  annotate("text", x=2.5, y=(25-22)/2+22, label= "D", color="red", cex=3)+
  annotate("text", x=2.5, y=(20-9)/2+9, label= "B", color="red", cex=3)+
  annotate("text", x=2.5, y=(8-1)/2+1, label= "C", color="red", cex=3)

dendro_rust_groups

# Arranging all three graphs together with the dendrogram groups clearly marked
arrangeddiffplotswgroups <- ggarrange(dendro_rust_groups, heatmap_NvsS, heatmap_mean_diff,
                               labels=c("", "a", "b"),
                               nrow=1,
                               align="h",
                               widths=c(0.035,1,0.22),
                               common.legend = TRUE,
                               legend="bottom")
arrangeddiffplotswgroups <- annotate_figure(arrangeddiffplotswgroups,
                left = text_grob("Differences in virulence percentages from North to South", rot = 90, size=10))

ggsave("Figure 6 N vs S heatmap 450dpi.tiff", scale=1, dpi=450, width=246, height=178, units="mm", path=".", bg="white")

# cutting out the decade for the poster
arrangeddiffplotswgroups <- ggarrange(dendro_rust_groups, heatmap_NvsS,
                                      nrow=1,
                                      align="h",
                                      widths=c(0.035,0.95),
                                      common.legend = TRUE,
                                      legend="bottom")
arrangeddiffplotswgroups <- annotate_figure(arrangeddiffplotswgroups,
                                            left = text_grob("Differences in virulence percentages from North to South", rot = 90, size=10))

ggsave("Figure 6 N vs S heatmap 450dpi.tiff", scale=1, dpi=450, width=246, height=178, units="mm", path=".", bg="white")


# Boxplot comparing the differences in virulences between N and S  --------


str(OCR_wout_2013_clean)
NvS_boxplot <-ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2023), y=countvirulence, fill=region, group = interaction(year, region))) + 
  # geom_boxplot(alpha=samples_per_year$n_adjusted, outlier.alpha = 1, outlier.size=0.5)+
  geom_boxplot(outlier.size=0.5)+
  xlab("Year")+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2023), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 13.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 13.5, xmax = 18.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot

NvS_boxplot

## Wilcox Test -------------------------------------------------------------

library(rstatix)

OCR_wout_2013_long_ttest <- OCR_wout_2013_clean %>%
  select(isolate, year, region, countvirulence)

wilcox.stat.test <- OCR_wout_2013_long_ttest %>%
  group_by(year) %>%
  wilcox_test(countvirulence ~ region) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
wilcox.stat.test %>% print(n=30)

wilcox.stat.test <- wilcox.stat.test %>%
  add_xy_position(x = "year")%>%
  select(-x) %>% #remove x because improperly calculated and will not show up in the right place on the graph
  add_column(x=1, .before = 'xmin') %>%
  mutate(x= year-1992)%>%
  mutate(xmin= x-0.2) %>%
  mutate(xmax= x+0.2)
wilcox.stat.test

# Add the significances to the graph
str(OCR_wout_2013_clean)
ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2023), y=countvirulence, fill=region, group = interaction(year, region))) +
  geom_boxplot(outlier.size=0.5)+
  xlab("Year")+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2023), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 13.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 13.5, xmax = 18.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  stat_pvalue_manual(wilcox.stat.test, label = "p.adj.signif", inherit.aes = FALSE, bracket.nudge.y= -0.5, hide.ns=TRUE)

ggsave("Sup Figure S2 NvsS.tiff", scale=1, dpi=300, width=2500, height=1800, units="px", path=".")

# #now make the same graph but for combining with barplot below
NvS_boxplot <-ggplot(OCR_wout_2013_clean, aes(x=factor(year, levels=1993:2023), y=countvirulence, fill=region, group = interaction(year, region))) +
  # geom_boxplot(alpha=samples_per_year$n_adjusted, outlier.alpha = 1, outlier.size=0.5)+
  geom_boxplot(outlier.size=0.5)+
  ylab("Number of virulences per rust isolate")+
  ggtitle("Virulences per rust isolate")+
  theme_bw()+
  theme(legend.position="none", axis.title.x=element_blank())+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2023), drop=FALSE)+
  scale_y_continuous(limits = c(0,40), expand = c(0, 0))+
  # geom_jitter(color="black", size=0.2, alpha=0.3) + #not as much of a fan of jitter for this plot
  annotate("rect", xmin = 0.4, xmax = 13.5, ymin = 30, ymax = 40,fill = "grey")+ #adding grey rectangles because didn't have all 40 lines in the beginning
  annotate("rect", xmin = 13.5, xmax = 18.5, ymin = 36, ymax = 40,fill = "grey")+
  annotate("point", x = 21, y = 16, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  annotate("point", x = 21, y = 24, size=0.5) + #manually adding two observation dots to 2013 because it originally tried to turn them into a boxplot
  stat_pvalue_manual(wilcox.stat.test, label = "p.adj.signif", inherit.aes = FALSE, bracket.nudge.y= -0.5, hide.ns=TRUE)
NvS_boxplot

# Graphing the data for north and south in order to compare the sample sizes by year
samples_per_year_perRegion <- OCR_survey_clean %>% count(year, region, sort = FALSE)
samples_per_year_perRegion


Sample_Number_by_region_bargraph <-ggplot(samples_per_year_perRegion, aes(x=factor(year, levels=1993:2023), y=n, fill=region)) +
  geom_bar(stat='identity', position='dodge', width= 0.65)+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), breaks=factor(1993:2023), drop=FALSE)+
  xlab("Year")+
  ylab("Number of isolates by region")
Sample_Number_by_region_bargraph


arrangedNvsS <- ggarrange(NvS_boxplot, Sample_Number_by_region_bargraph,
                          labels=c("a", "b"),
                          nrow=2,
                          align="v")
arrangedNvsS
ggsave("Supplementary Figure NvsSarranged.tiff", scale=1, dpi=300, width=4000, height=2300, units="px", path=".")


# Comparing N vs S four years of data -------------------------------------

str(OCR_wout_2013_clean)
heatmap_data <- OCR_wout_2013_clean %>% group_by(year, region) %>% #using two variables for grouping
  summarise(across("Pc14":"Stainless", ~ mean(.x, na.rm = TRUE)*100), .groups="keep") %>% 
  filter(year>2019)

### change data into long format
heatmap_data_long <-gather(heatmap_data, differential_line, percent_virulent, Pc14:Stainless) %>% 
  unite("year_region", year:region, sep= " ", remove=FALSE)
str(heatmap_data_long)

#createheatmap
heatmap_year <-ggplot(heatmap_data_long, aes(factor(year_region), factor(differential_line, levels=rev(unique(differential_line))), fill=percent_virulent))+
  geom_tile(show.legend = FALSE)+
  scale_fill_continuous_sequential(palette = "YlOrRd", rev = TRUE)+
  scale_y_discrete(limits = colnames(heatmap_data_matrix)[clust$order])+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = round(percent_virulent, 1)), color = "black", size = 3)+
  ylab("Differential Line")+
  ggtitle("Percent virulent isolates per year by region")+
  annotate("segment", x = 2.5, xend = 2.5, y = 40.5, yend = 0.5, colour = "black", linewidth=1, alpha=0.4)+
  annotate("segment", x = 4.5, xend = 4.5, y = 40.5, yend = 0.5, colour = "black", linewidth=1, alpha=0.4)+
  annotate("segment", x = 6.5, xend = 6.5, y = 40.5, yend = 0.5, colour = "black", linewidth=1, alpha=0.4)+
  scale_x_discrete(labels=c("2020 North" = "2020\nNorth", "2020 South" = "2020\nSouth", "2021 North" = "2021\nNorth", "2021 South" = "2021\nSouth",
                            "2022 North" = "2022\nNorth", "2022 South" = "2022\nSouth", "2023 North" = "2023\nNorth","2023 South" = "2023\nSouth"))
heatmap_year

dendro_rust_groups

arrangedplotsordered <- ggarrange(dendro_rust_groups, heatmap_year,
                                  nrow=1,
                                  align="h",
                                  widths=c(0.1, 1))

annotate_figure(arrangedplotsordered,
                left = text_grob("Differential Line", rot = 90))

ggsave("Supplemental Figure heatmapNvS 2020-2023.tiff")

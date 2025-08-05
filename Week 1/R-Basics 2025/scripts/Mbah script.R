# Creating working Directory and R project ----
getwd()                                     #To check current working directory
setwd('provide path to the desired folder') # Selecting the desired working directory
setwd("C:/Users/Admin/OneDrive/Desktop/R-Basics 2025")  #I have select R-basics in my desktop as my working directory
# create a new project at the top right and choose the already existing directory (R-Basics 2025) in the desktop
#Click on create. Reload the R-data.


# Organizing Project folders ----
# create sub-folders including data, scripts and results where files would be stored
dir.create('Data') 
dir.create('scripts')
dir.create('results')


#Set repositories and install packages to load data ----
setRepositories()                   #set repository to be able to install packages
installed.packages()                #to download and install R packages
library()                           #to load packages


# Import data from different file formats ----
library(readr)                      #load the package required to import csv or tsv files into R
library(readxl)                     #load the package required to load excel file into R
patient_info <- read.csv(file.choose()) #open and select the csv data
View(patient_info)
data <- read_tsv(file.choose())     #open and select the tsv file format
code <- read.delim(file.choose())   #read a file (either csv,tsv) into a tibble
?(package)                          # Opens Help page to know how a package is used




# Data set structure and Conversions----
str(patient_info)                   #Check the structure of the object

#As factor
gender <- as.factor(patient_info $gender) 
str(gender)
print (patient_info$gender)
diagnosis <- factor(patient_info$diagnosis, levels = c('Cancer', 'Normal'), labels = c(1, 2))
str(diagnosis)
smoker <- as.factor(patient_info$smoker)
smoker_binary <- ifelse(patient_info$smoker == 'yes', 1, 0)
smoker_binary <- as.factor(smoker_binary)

#As list
my_list <- list(patient_info$age, patient_info$gender)
?list

#As Numeric
bmi <- as.numeric(patient_info$bmi)


# we can evaluate data structures with logical expressions
myLogical <- is.numeric(patient_info$diagnosis)
class(diagnosis)


#Saving Output----
# To save R script use keyboard shortcut 
# ctrl + s 
# Or go to the file menu & select 'Save As' option

# To open R script use keyboard shortcut 
# ctrl + Enter 
# Or go to the file menu & select 'Open File' option

# Save file  as CSV in results folder
write.csv(patient_info, file = "results/patient_data.csv")

# Save the entire R workspace
save.image(file = "Mbah_chinedu_Class_1b.RData")

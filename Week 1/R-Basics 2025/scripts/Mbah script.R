# Creating working Directory and R project ----
getwd()                                     #To check current working directory
setwd('provide path to the desired folder') # Selecting the desired working directory
setwd("C:/Users/Admin/OneDrive/Desktop/MyBioinformaticsJourney2025/Week 1/R-Basics 2025")  #I have select R-basics in my desktop as my working directory
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
data <- read.csv(file.choose()) #open and select the csv data
View(patient_info)
data <- read_tsv(file.choose())     #open and select the tsv file format
code <- read.delim(file.choose())   #read a file (either csv,tsv) into a tibble
?(package)                          # Opens Help page to know how a package is used




# Data set structure and Conversions----
# before analyzing any dataset, you need to ensure that they are stored in the right format

str(data)                   #Check the structure of the object
clean.data <- data

#As factor (Convert gender, diagnosis, and smoker to factor)
clean.data$gender <- as.factor(clean.data$gender)
print(levels(clean.data$gender))
print(class(clean.data$gender))
str(clean.data$gender)

clean.data$gender <- factor(clean.data$gender)  #convert to factors with levels
str(clean.data$gender)

clean.data$gender_num <- as.numeric(clean.data$gender) #convert from factors to numeric

# FOR DIAGNOSIS
clean.data$diagnosis <- factor(clean.data$diagnosis, levels = c("Cancer", "Normal"))
str(clean.data$diagnosis)

# FOR SMOKER
clean.data$smoker <- as.factor(clean.data$smoker)
str(clean.data$smoker)

clean.data$gender_num <- NULL #helps to remove column

smoker_binary <- ifelse(clean.data$smoker == 'yes', 1, 0)
smoker_binary <- as.factor(smoker_binary)
str(clean.data)

#As list
my_list <- list(clean.data$age, clean.data$gender)
?list

#As Numeric
clean.data$bmi <- as.numeric(clean.data$bmi)


# we can evaluate data structures with logical expressions
myLogical <- is.numeric(clean.data$diagnosis)
class(clean.data$diagnosis)


#Saving Output----
# To save R script use keyboard shortcut 
# ctrl + s 
# Or go to the file menu & select 'Save As' option

# To open R script use keyboard shortcut 
# ctrl + Enter 
# Or go to the file menu & select 'Open File' option

# Save file  as CSV in results folder

write.csv(clean.data, file = "results/transformed.data.csv")

# Save the entire R workspace

save.image(file = "results/Mbah_chinedu_Class_1b.RData")

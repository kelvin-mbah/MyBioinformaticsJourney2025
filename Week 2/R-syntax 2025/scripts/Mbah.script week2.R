#Operators in R ----
#operators are used to store variables in R (<- and =)
ose_height <- c(5.6, 5.8, 6.3, 6.5, 6.7)
ose_weight <- c(50, 55, 53, 59, 54)
allergy_status <- c("yes", "no", 'yes', 'no', 'yes')

#Arithmetic
#Symbols are used to perform arithmetic operations 

bmi <- ose_weight/(ose_height)^2

#comparisons
# Comparison operators ask TRUE/FALSE questions about values
# They do not calculate a number, they return a logical output
ose_weight > 50
ose_height >= 5.4

allergy_status == 'no'
allergy_status != 'no'

data <- data.frame(ose_height, ose_weight, allergy_status, bmi)

 data[, 4] <- 0

 

#Data Pre-processing/Cleaning/Assessment ----
 #Real world data set are always or sometimes messy
 #so we have to clean or transform such data set
 #Note: always check the structure of the data set before analyzing
 
 str(data)            #Provides a preview of the data (column names and data types)
 head(data)           #checks the first six rows (observations)
 tail(data, n = 4)    #shows the last four observation (rows) since we specified using n
 dim(data)            #provides the dimensions (Row, Column) of the data
 
 #To show the names of variables
 colnames(data)
 rownames(data)
 names(data)
 
 #we can access data variables using: dollar symbol ($) or index []
 data$ose_weight     #extract a single column or variable from the entire dataset
data[, 2] 
data [1, c(1:3)]
data$ose_language <- c('English', 'Spanish', 'Pidgin', 'yoruba', 'igbo') #add new variables or columms

#Missing Values
is.na(data)         #Identify missing values
sum(is.na(data))    #Gives total missing values
colSums(is.na(data))    #Missing values in column
rowSums(is.na(data))    #Missing values in row

#Handling/Removing Na: 
Data <- data[, colSums(is.na(data)) == 0]   

# replace NA value with 0
Data1 <- data
Data1[is.na(Data1)] <- 0

# replace NA value with mean
Data2 <- data
Data[is.na(Data2)] <- mean(data$bmi, na.rm = TRUE)


#FUNCTIONS IN R ----

#Function is like a reusable tool used to do a specific job quickly
# You create it once and let it do the work whenever you need it
# HOW TO CREATE A FUNCTION:
        # Name: This is the name given to your function
        # Ingredients: This are the input (argument) you give the function to work on
        # Steps: This contains the instructions that tell the function what to do
        # Results : This is the output(return value) after the function does its job

#example:

classify_gene <- function(logFC, padj){
  if(logFC >1 && padj < 0.05) {
    return('upregulated')
  } else if(logFC < -1 && padj < 0.05) {
    return('downregulated')
  } else {
    return('Not Significant')
  }
}
  
  data <- read.csv(file.choose())
str(data)

(classify_gene(1.5, 0.01))
(classify_gene(-1.5, 0.01))



#LOOPS ----

#Instead of repeating steps for each file, we use loops.

# Typical loop workflow:

# Define the input folder (where raw data files are stored) and the output folder (where results will be saved).
# Specify the list of files that need to be processed.
# Prepare an empty list in R to store results for later use.
# For each file in the list:
#          Import the data into R.
#          Check and handle missing values (NA).
#          Calculate BMI using the calculate_BMI function.
#          Save the processed results both as a CSV file and in the R results list.

# Define input and output folders
input_dir <- 'data'
output_dir <- 'result'

#Create output folder if not created
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

#Specify the list of files that need to be processed
#files name must tally with the files you want to work on
raw_files <- c('DEGs_Data_1.csv', 'DEGs_Data_2.csv')


# prepare an empty list
result_list <- list()


for(file_names in raw_files) {
  cat('processing', file_names, ': \n')
  input_file_path <- file.path(input_dir, file_names)

  # import the dataset
  data <- read.csv(input_file_path, header = TRUE)  
  cat('files imported: checking for missing values...\n')
  
  #Handling missing values in padj
  if('padj' %in% names(data)) {
    missing_count <- sum(is.na(data$padj))
   
     cat('missing value in padj:', missing_count, '...\n')
    data$padj[is.na(data$padj)]  <- 1
  }
  
  #Handling missing values in logFC
  if('logFC' %in% names(data)) {
    missing_count <- sum(is.na(data$logFC))
    
    cat('missing value in logFC:', missing_count, '...\n')
    data$logFC[is.na(data$logFC)]  <- 1
  }
  
  #Classify Gene using Function
  data$status <- mapply(classify_gene, data$logFC, data$padj, SIMPLIFY = TRUE)
  cat('Genes have been classified successfully:...\n')
  
  #save result in R
  result_list[[file_names]] <- data
  
  #Save result in output directory
  output_file_path <- file.path(output_dir, paste0('DEG_result', file_names))
  write.csv(data, output_file_path, row.names = FALSE)
  cat('Results saved to', output_file_path, '...\n')
}


# The loop repeats until all files are processed.

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]

# save FILE
save.image(file = "result/Mbah_chinedu_Class_2b.RData")

# Summary----

# Loops automate repetitive work 
# making your workflow faster 
# consistent, and reproducible
###############################################################################
# Choose method to specify which peptide should be used in parent mass list:
# "all"         : don't discard peptides
# "STY-phospho" : keep only peptides with phosphorylation on S,T, or Y
# "K-acetyl"    : keep only peptides with acetylation on K
#
# ChargeFilter defines the method that is applied for charge state filtering:
#   - "optimal" : only charge states are allowed that are likely to be detected
#                 according to the number of basic amino acids and
#                 phosphorylations (empiric by dataset 15-137)
#   - "oneMax"  : only the one most likely charge State is accepted
#   - "noFilter": all charge states are accepted
#
# MSX defines the number of MSX used in your PRM method (0 mean disabled), MSX
# groups are defined by RT if present, otherwise by m/z. RT will be deleted.
#   - "1-10"    : peptides are gouped in x MSX IDs
#   - "0"       : no MSX IDS
###############################################################################


# load packages

if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

# functions
countCharOccurance <- function(char, text, allowNA = TRUE) {
  text2 <- gsub(pattern = char, replacement = substr(x = char,
                                                     start = 1,
                                                     stop = (nchar(char)-1)),
                x = text,
                fixed = TRUE)
  return (nchar(text) - nchar(text2))
}

print("test1")

# empirical charge state reference table (Data: 15-137)
basicSites <- rep(seq(0, 8),5)

phosphoSites <- c(rep(0, 9), rep(1,9), rep(2,9), rep(3,9), rep(4,9))

optimalStart <- c(rep(c(2,2,2,2,3,4,4,5,5),2),
                       rep(c(rep(2,6), rep(3,3)),3))

optimalStop <- c(rep(c(2,2,3,4,4,5,6,6,6),2),
                      c(2,2,3,4,4,4,5,5,5),
                      c(2,2,3,3,4,4,5,5,5),
                      c(2,2,3,3,3,3,4,4,4))

oneCharge <- c(c(2,2,3,3,3,4,5,5,6),
               rep(c(2,2,2,3,3,4,4,5,5), 2),
               rep(c(2,2,2,2,3,4,4,5,5),2))

ChargeRefTable <- data.frame(basicSites,
                             phosphoSites,
                             optimalStart,
                             optimalStop,
                             oneCharge)


print("test2")


# Input
inputFile <- file.choose()

method <- select.list(choices = c("all", "STY-phospho", "K-acetyl", "K-GlyGly"),
                      graphics = TRUE,
                      title = "Which peptides to keep:")

ChargeFilter <- select.list(choices = c("optimal", "oneMax", "noFilter"),
                            graphics = TRUE,
                            title = "Choose Charge state filter:")
MSX <- select.list(choices = c("0","1", "2", "3", "4", "5",
                               "6", "7", "8", "9", "10"),
                   graphics = TRUE,
                   title = "Choose MSX count (0 = off):")

MSX <- as.integer(MSX)


# load file
input <- read.csv(file = inputFile)

input$Comment <- as.character(input$Comment)

# filter by method
if(method == "STY-phospho"){
  
  # keep only "[+80]"
  input <- input[grep(pattern = "[+80.0]",
                      x = input$Comment,
                      ignore.case = FALSE,
                      fixed = TRUE), ]
  
} else {
  
  if(method == "K-acetyl"){
    
    # keep only "[+42.0]"
    input <- input[grep(pattern = "[+42.0]",
                        x = input$Comment,
                        ignore.case = FALSE,
                        fixed = TRUE), ]
    
    # delete c-terminal modifications
    input$NewComment <- substr(x = input$Comment,
                            start = 1,
                            stop = regexpr(pattern = "\\(",
                                           text = input$Comment) - 2)
    
    input <- input[!grepl(pattern = "[+42.0]",
                          x = substr(input$NewComment,
                                     start = nchar(input$NewComment) - 6,
                                     stop = nchar(input$NewComment)),
                          fixed = TRUE),
                   - dim(input)[2]]
    
  } else {
    
    if(method == "K-GlyGly"){
      
      # keep only "[+114.0]"
      input <- input[grep(pattern = "[+114.0]",
                          x = input$Comment,
                          ignore.case = FALSE,
                          fixed = TRUE), ]
      
      # delete c-terminal modifications
      input$NewComment <- substr(x = input$Comment,
                                 start = 1,
                                 stop = regexpr(pattern = "\\(",
                                                text = input$Comment) - 2)
      
      input <- input[!grepl(pattern = "[+114.0]",
                            x = substr(input$NewComment,
                                       start = nchar(input$NewComment) - 7,
                                       stop = nchar(input$NewComment)),
                            fixed = TRUE),
                     - dim(input)[2]]
      
    } else {
      
      if(method != "all"){
        
        print("Undefined Method has been chosen. Method was set to all.")
        
        method <- "all"
        
      }
    }
  }
}

# get basic sites
input$basicSites <- countCharOccurance(char = "R",
                                       text = input$Comment) +
  countCharOccurance(char = "H", text = input$Comment) +
  countCharOccurance(char = "K", text = input$Comment) -
  countCharOccurance(char = "K[+42.0]", text = input$Comment) -
  countCharOccurance(char = "K[+114.0]", text = input$Comment) 

# get number of PhosphoSites
input$phosphoSites <- str_count(pattern = "80.0",
                           string = input$Comment)

# get optimal charge states from ReferenceTable
isolationList <- merge(x = input,
                       y = ChargeRefTable,
                       by = c("basicSites", "phosphoSites"))

# filter charge states according to ChargeFilter
if(ChargeFilter == "optimal"){
  
  isolationList <- isolationList[
    isolationList$CS..z. >= isolationList$optimalStart &
      isolationList$CS..z. <= isolationList$optimalStop, ]
  
} else {
  
  if(ChargeFilter == "oneMax"){
    
    isolationList <- isolationList[
      isolationList$CS..z. == isolationList$oneCharge, ]
    
  } else{
    
    if(ChargeFilter != "noFilter"){
      
      print("No defined charge filter has been chosen. Was set to noFilter")
      
      ChargeFilter <- "noFilter"
      
    }
  }
}

# filter out double m/z entries and sort by  
isolationList <- isolationList[!duplicated(isolationList$Mass..m.z.), ] 

isolationList <- isolationList[order(isolationList$Mass..m.z.), ]

# import file with original column names
input2 <- read.csv(file = inputFile, check.names = FALSE)

# Use optimal MSX IDs

if(is.integer(MSX) & MSX > 0){
  
  # check if RT are present and sort by them
  if(!any(is.na(isolationList$Start..min.))){
    
    isolationList <- isolationList[order(isolationList$Start..min.), ]
  }
  # give MSX IDs optimal to RT or m/z
  MSXIDs <- rep(1:ceiling(dim(isolationList)[1]/MSX), MSX)
  
  isolationList$MSXID <- MSXIDs[1:dim(isolationList)[1]]
  
  # clean data frame
  isolationList <- isolationList[ ,c(-1,
                                   -2,
                                   -(dim(isolationList)[2]-1),
                                   -(dim(isolationList)[2]-2),
                                   -(dim(isolationList)[2]-3))]
  
  # reorder data frame
  isolationList <- isolationList[ , c(1,2,3,4,5,6,7,8,10,9)]
  
  # delete RTs
  isolationList$Start..min. <- NA
  
  isolationList$End..min. <- NA
  
  # check if MSX ID was already in original file and name accordingly
  if(!any(grepl(pattern = "MSX", x = names(input2), ignore.case = TRUE))){
    names(isolationList)[c(1,2,3,4,5,6,7,8,10)] <- names(input2)
    names(isolationList)[9] <- "MSX ID"
    
  } else {
    
    names(isolationList) <- names(input2)
    
  }
  
} else {
  
  # if MSX was set 0, clean and rename data frame
  isolationList <- isolationList[, c(-1,
                                   -2,
                                   -(dim(isolationList)[2]),
                                   -(dim(isolationList)[2]-1),
                                   -(dim(isolationList)[2]-2))]
  
  names(isolationList) <- names(input2)
  
}

# write Isolation file
write.csv(x = isolationList[ , - dim(input)[2]],
          file = paste0(dirname(inputFile),
                        "/",
                        substr(x = basename(inputFile),
                               start = 1,
                               stop = nchar(basename(inputFile))-4),
                        "_",
                        method,
                        "_",
                        ChargeFilter,
                        "_IsolatioList.csv"),
          row.names = FALSE)
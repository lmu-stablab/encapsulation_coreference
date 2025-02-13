#------------------------------
# DESCRIPTION
#------------------------------

# This script contains code for the pre-processing 
# of and descriptive analysis of the data.
# The processed data is saved in the folder 0_data/2_preparedData.

#-----------HEADER-------------
# load packages
library(knitr)        # output table editing etc.
library(readxl)       # reading Excel files
library(tidyverse)    # data management and ggplot2
library(DPKogHelpers) # Helper functions for the DPKog project
theme_set(theme_bw()) # set global theme for ggplot2 plots

# define paths
path_rawData  <- "0_data/1_rawData/"
path_prepData <- "0_data/2_preparedData/"

# response variables
response_vec <- c("TRT","FRT","RRT")

# list all raw data Excel files
files <- list.files(path_rawData)
files <- as.list(files[grepl("xls", tolower(files))])
# read data
dat_list <- lapply(files, 
                   function(file) readxl::read_excel(paste0(path_rawData, file), 
                                                     sheet = 1, na = ".")) 
#------------------------------

#------------------------------ DIMENSIONS
# Dimension of the data
strings <- lapply(seq_len(length(dat_list)), function(i) {
  paste0("\n", files[i], ": ", paste(dim(dat_list[[i]]), collapse = " x "))
})
string_dataDim <- paste(strings, collapse = "\n")
print(string_dataDim)

# Distribution of AOIs
tab_list <- lapply(seq_len(length(dat_list)), function(i) {
  dat <- as.data.frame(table(dat_list[[i]]$`IA_ID`)) %>% 
    dplyr::rename(AOI = Var1) %>% mutate(AOI = as.character(AOI))
  colnames(dat)[2] <- paste0("h_file",i)
  return(dat)
})
tab <- Reduce(left_join, tab_list)
tab[is.na(tab)] <- 0 # for if one AOI is not appearing in one dataset
kable(tab)
#------------------------------

#------------------------------ DATA PREPARATION
dat_list[[2]]$TRIAL_LABEL <- as.numeric(substr(dat_list[[2]]$TRIAL_LABEL, 8, 10))

# define observation and trials
for (i in seq_len(length(dat_list))) {
  dat_list[[i]]$Participant <- paste0(i, dat_list[[i]]$DATA_FILE)
  dat_list[[i]]$trial <- gsub("_(\\d+)_", "_", dat_list[[i]]$`CONDITION+AOI+SET`) 
  dat_list[[i]]$Topic <- substr(dat_list[[i]]$`CONDITION+AOI+SET`, 7, 8)
  dat_list[[i]]$observation <- paste0(dat_list[[i]]$trial, "_", dat_list[[i]]$Participant)
  dat_list[[i]]$condition <- substr(dat_list[[i]]$AOI_LABEL, 1, nchar(dat_list[[i]]$AOI_LABEL) - 2)
}

# create one big dataset
dat <- dplyr::bind_rows(dat_list)

# set missing FRT & RRT for TRT == 0 to 0
dat[dat$AOI_TOTAL_READING_TIME == 0 
    & is.na(dat$AOI_FIRST_PASS_READING_TIME) 
    & is.na(dat$AOI_REREADING_TIME), 
    c("AOI_FIRST_PASS_READING_TIME", "AOI_REREADING_TIME")] <- 0

# rename variables and define letters per word and AOI.conditions
dat <- dat %>%
  dplyr::rename("TRT"      = "AOI_TOTAL_READING_TIME",
                "FRT"      = "AOI_FIRST_PASS_READING_TIME",
                "RRT"      = "AOI_REREADING_TIME", 
                "nLetters" = "CHARACTER_COUNT",
                "nWords"   = "WORD_COUNT",
                "AOI"      = "IA_ID",
                "AOI.condition" = "AOI_LABEL") %>%
  mutate(TRT.WD = TRT / nWords,
         FRT.WD = FRT / nWords,
         RRT.WD = RRT / nWords,
         nLetters.WD = nLetters / nWords)

# ensure that all reading time variables are of class numeric
vars <- paste0(rep(response_vec, times = 2), rep(c("",".WD"), each = length(response_vec)))
for (var in vars)
  dat[,var] <- as.numeric(unlist(dat[,var]))
#------------------------------

#------------------------------ CONSISTENCY CHECKS
error_messages <- c()

# CHECK 1: missing values in the response
for (var in response_vec) {
  if (any(is.na(dat[,var]))) {
    error_messages <- append(error_messages, paste0("Missing values exist in ",var,"!"))
  }
}

paste0(sum(is.na(dat$FRT)), " observations with missing FRT and RRT are removed.")
dat <- dat %>% filter(!is.na(FRT), !is.na(RRT))
error_messages <- c()

# Check pre-calculated reading times per word
for (var in response_vec) {
  if (any(dat[,var]/dat$nWords != dat[,paste0(var,".WD")])) {
    error_messages <- append(error_messages,
                             c(paste0("The values of ",var,".WD differ from ",var,"/nWords!"),
                               paste0("  -> Calculate values of ",var,".WD anew.")))
    dat[,paste0(var,".WD")] <- dat[,var]/dat$nWords
  }
}

# Final message
if (length(error_messages) == 0)
  error_messages <- "No inconsistencies found."

print(error_messages)

# Distribution of AOI.conditions
tab <- table(dat$AOI.condition) %>% as.data.frame() %>%
  dplyr::rename(AOI.condition = Var1, h = Freq)
kable(tab)
#------------------------------

#------------------------------ OUTLIER HANDLING
# 1. Any.first.skip: First Pass Reading Time in AOI 1 = 0ms 
# 2. Any.fast.reader (extremely low values): First Pass Reading Time in AOI 1 <80ms 
#    and Second Pass Reading Time in AOI 1 <80ms
# 3. Any.slow.reader (extremely high values): Total Reading Time in AOI 1 >800ms

# 1. Any.first.skip
messages <- c()
aois1 <- c(1) 
outlier1 <- lapply(aois1, function(aoi) { dat %>% filter(AOI == aoi, FRT == 0)})
for (i in seq_len(length(aois1))) {
  new_message <- paste0(nrow(outlier1[[i]])," outliers in AOI '",aois1[i],"'.",
                        "\n",length(unique(outlier1[[i]]$observation)),
                        " (",paste0(100*round(length(unique(outlier1[[i]]$observation)) / 
                                                length(unique(dat$observation)),4),"%"),
                        ") of all observations are affected.")
  messages <- append(messages, new_message)
}
print(messages)

# 2. Any.fast.reader
messages <- c()
aois2 <- c(1)
outlier2 <- lapply(aois2, function(aoi) { dat %>% filter(AOI == aoi, FRT.WD < 80, RRT.WD < 80)}) 
for (i in seq_len(length(aois2))) {
  new_message <- paste0(nrow(outlier2[[i]])," outliers in AOI '",aois2[i],"'.",
                        "\n",length(unique(outlier2[[i]]$observation)),
                        " (",paste0(100*round(length(unique(outlier2[[i]]$observation)) / 
                                                length(unique(dat$observation)),4),"%"),
                        ") of all observations are affected.")
  messages <- append(messages, new_message)
}
print(messages)

# 3. Any.slow.reader
messages <- c()
aois3 <- c(1)
outlier3 <- lapply(aois3, function(aoi) { dat %>% filter(AOI == aoi, TRT.WD > 800) }) 
for (i in seq_len(length(aois3))) {
  new_message <- paste0(nrow(outlier3[[i]])," outliers in AOI '",aois3[i],"'.",
                        "\n",length(unique(outlier3[[i]]$observation)),
                        " (",paste0(100*round(length(unique(outlier3[[i]]$observation)) / 
                                                length(unique(dat$observation)),4),"%"),
                        ") of all observations are affected.")
  messages <- append(messages, new_message)
}
print(messages)


# OUTLIER DELETION
outlier <- dplyr::bind_rows(outlier1, outlier2, outlier3)
dat_new <- dat %>% filter(!(observation %in% unique(outlier$observation)))

# DESCRIPTION
cat(sprintf(
  "Outliers were detected in %d of %d participants.\n\n",
  length(unique(outlier$Participant)),
  length(unique(dat$Participant))
))

cat(sprintf(
  "Overall, %d (i.e. %s) of %d observations were excluded.\n",
  length(unique(outlier$observation)),
  paste0(100 * round(length(unique(outlier$observation)) / length(unique(dat$observation)), 3), "%"),
  length(unique(dat$observation))
))

# Negative values
error_messages <- c()
for (var in response_vec) {
  if (any(dat[,paste0(var,".WD")] < 0))
    error_messages <- append(error_messages, paste0("--- Negative values exist in ",var,"!"))
}
if (length(error_messages) == 0)
  error_messages <- "No negative reading times found."
print(error_messages)


# Distribution after exclusion of outliers
dat_new$AOI <- as.character(dat_new$AOI)
vars <- c("AOI", paste0(response_vec, ".WD"))
plot_dat <- dat_new %>% select(!!vars) %>%
  gather("response", "value", -AOI)
ggplot(plot_dat, aes(x = AOI, y = value)) +
  geom_boxplot() +
  facet_wrap(~ response, nrow = length(response_vec), scales = "free")

# save prepared data as .csv
write.csv2(dat_new, file = paste0(path_prepData, "data_prepared.csv"), row.names = F, fileEncoding = "utf-8")
message <- paste0("The prepared data were successfully saved as '",path_prepData,"data_prepared.csv'.")















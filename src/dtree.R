################## In the name of Allah #######################
#                                                             #
#          Introduction to Bioinformatics, bahman 1400        #
#                    Sharif University                        #
#                                                             #
#  Name:    Mohammad Hosein Movasaghinia                      #
#  Uni-num: 400200919                                         #
#                                                             #
##############################################################*
#                       Decision Tree                         #
##############################################################*

library(GEOquery)
library(caTools)
library(party)
library(dplyr)
library(magrittr)

#### Set working directory ####
curD <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(sub(paste0("/",sub("(.+)/","",curD)),"",curD))

#### Set global variables ####
series <- "GSE48558"
platform <- "GPL6244"

#### Load data ####
gset <- getGEO(series, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#### Set group ####
gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")

sml <- strsplit(gsms, split="")[[1]]

#### filter by X and add Group to data ####
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
gs <- factor(sml)
groups <- make.names(c("normal", "test"))
levels(gs) <- groups
gset$group <- gs

sname <- gset$source_name_ch1
sname <- sub(" ", "", sname)
sname <- gsub("[+]", "p", sname)
sname.gs <- factor(sname)

# Expression of gset and scaleing
ex <- exprs(gset)
ex.scale <- t(scale(t(ex), scale = F))
pc <- prcomp(ex.scale)
plot(pc)
pcr <- data.frame(pc$r[,1:3], group = gset$group)


sample_data = sample.split(pcr, SplitRatio = 0.75)
train_data <- subset(pcr, sample_data == TRUE)
test_data <- subset(pcr, sample_data == FALSE)

model<- ctree(formula = group ~ ., data=train_data, controls = ctree_control(maxdepth = 10))
plot(model)


predict_model<-predict(model, test_data)

m_at <- table(test_data$group, predict_model)
m_at

ac_Test <- sum(diag(m_at)) / sum(m_at)
print(paste('Accuracy for test is found to be', ac_Test))


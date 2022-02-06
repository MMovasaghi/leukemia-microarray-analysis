################## In the name of Allah #######################
#                                                             #
#          Introduction to Bioinformatics, bahman 1400        #
#                    Sharif University                        #
#                                                             #
#  Name:    Mohammad Hosein Movasaghinia                      #
#  Uni-num: 400200919                                         #
#                                                             #
##############################################################*
#   SVM for Classification of Normal vs. AML Patient Cells    #
##############################################################*

install.packages('caTools')
install.packages('e1071')
library(e1071)
library(caTools)
library(GEOquery)

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
pcr <- data.frame(pc$r[,1:2], group = gset$group)

set.seed(123)
split = sample.split(pcr$group, SplitRatio = 0.75)

training_set = subset(pcr, split == TRUE)
test_set = subset(pcr, split == FALSE)


classifier = svm(formula = group ~ .,
                 data = training_set,
                 type = 'C-classification',
                 kernel = 'sigmoid',
                 degree=4,
                 coef0=1)

y_pred = predict(classifier, newdata = test_set)
cm = table(test_set[, 3], y_pred)
cm


# Plotting the training data set results
set = training_set
X1 = seq(min(set[, 1]) - 1, max(set[, 1]) + 1, by = 0.01)
X2 = seq(min(set[, 2]) - 1, max(set[, 2]) + 1, by = 0.01)

grid_set = expand.grid(X1, X2)
colnames(grid_set) = c('PC2', 'PC1')
y_grid = predict(classifier, newdata = grid_set)

plot(set[, -3],
     main = 'SVM (Training set)',
     xlab = 'PC2', ylab = 'PC1',
     xlim = range(X1), ylim = range(X2))
legend("Normal", "Test")

contour(X1, X2, matrix(as.numeric(y_grid), length(X1), length(X2)), add = TRUE)
points(grid_set, pch = '.', col = ifelse(y_grid == 'test', 'coral1', 'aquamarine'))
points(set, pch = 21, bg = ifelse(set[, 3] == 'test', 'green4', 'red3'))



set = test_set
X1 = seq(min(set[, 1]) - 1, max(set[, 1]) + 1, by = 0.01)
X2 = seq(min(set[, 2]) - 1, max(set[, 2]) + 1, by = 0.01)

grid_set = expand.grid(X1, X2)
colnames(grid_set) = c('PC1', 'PC2')
y_grid = predict(classifier, newdata = grid_set)

plot(set[, -3], main = 'SVM (Test set)',
     xlab = 'PC1', ylab = 'PC1',
     xlim = range(X1), ylim = range(X2))

contour(X1, X2, matrix(as.numeric(y_grid), length(X1), length(X2)), add = TRUE)
points(grid_set, pch = '.', col = ifelse(y_grid == 'test', 'coral1', 'aquamarine'))
points(set, pch = 21, bg = ifelse(set[, 3] == 'test', 'green4', 'red3'))




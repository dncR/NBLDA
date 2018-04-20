source("R/all_classes.R")
source("R/functions.R")

tmp <- dir("R/")

for (i in tmp){
  source(paste("R/", i, sep = ""))
}

set.seed(110418)
dat <- generateCountData(n = 100, p = 200, K = 2, param = 1, sdsignal = 0.5, DE = 0.2, allZero.rm = FALSE, tag.samples = TRUE)
x <- t(dat$x + 1)
xte <- t(dat$xte + 1)
y <- as.factor(dat$y) #factor(c("T", "N")[dat$y])
yte <- as.factor(dat$yte) #factor(c("T", "N")[dat$yte])

rpt = 2
nf = 5
allFolds <- NULL

set.seed(110418)
# Each element of the returned list is the indices of test samples in this fold.
allFolds <- lapply(1:rpt, function(x){
  tmp <- balanced.folds(y, nfolds = nf)
  names(tmp) <- paste("Fold.", 1:nf, sep = "")
  tmp
})
names(allFolds) <- paste("Repeat.", 1:rpt, sep = "")


# PLDA vs NBLDA
ctrl <- ctrl2 <- nbldaControl(folds = 5, repeats = 2, foldIdx = allFolds, phi.epsilon = 0.15, transform = FALSE,
                         truephi = NULL, target = 0, normalize.target = TRUE, delta = NULL, return.selected.gene.names = TRUE)

ctrl2$truephi <- 0

# pldaRes <- trainPLDA(x = x, y = y, type = "deseq", transform = FALSE, nfolds = 5, repeats = 2, foldIdx = allFolds, tuneLength = 20)
pldaRes <- trainNBLDA(x = x, y = y, type = "deseq", tuneLength = 20, metric = "accuracy", train.control = ctrl2)
pred1 <- predict.nblda(pldaRes, xte, "predictions")

pldaRes$tuning.results[which(pldaRes$tuning.results[ ,1] == pldaRes$best.rho), ]
confusionMatrix(table(Predicted = pred1, Actual = yte))$overall[1]  # PLDA

nbldaRes <- trainNBLDA(x = x, y = y, type = "deseq", tuneLength = 10, metric = "accuracy", train.control = ctrl)
pred2 <- predict.nblda(nbldaRes, xte, "every")

nbldaRes$tuning.results[which(nbldaRes$tuning.results[ ,1] == nbldaRes$best.rho), ]
confusionMatrix(table(Predicted = pred2, Actual = yte))$overall[1]   # NBLDA

source("R/all_classes.R")
source("R/functions.R")

tmp <- dir("R/")

for (i in tmp){
  source(paste("R/", i, sep = ""))
}

set.seed(110418)
dat <- generateCountData(n = 30, p = 80, K = 2, param = 1, sdsignal = 0.5, DE = 0.6, allZero.rm = FALSE, tag.samples = TRUE)
x <- t(dat$x + 1)
xte <- t(dat$xte + 1)
y <- as.factor(dat$y) #factor(c("T", "N")[dat$y])
yte <- as.factor(dat$yte) #factor(c("T", "N")[dat$yte])

rpt = 3
nf = 5
allFolds <- NULL

set.seed(110418)
# Each element of the returned list is the indices of test samples in this fold.
allFolds <- lapply(1:rpt, function(x){
  tmp <- NBLDA:::balanced.folds(y, nfolds = nf)
  names(tmp) <- paste("Fold.", 1:nf, sep = "")
  tmp
})
names(allFolds) <- paste("Repeat.", 1:rpt, sep = "")


# PLDA vs NBLDA
ctrl <- ctrl2 <- nbldaControl(folds = 5, repeats = 3, foldIdx = allFolds, phi.epsilon = 0.1, transform = FALSE,
                         truephi = NULL, target = 0, normalize.target = TRUE, delta = NULL, return.selected.gene.names = TRUE)

ctrl2$truephi <- 0

# pldaRes <- trainPLDA(x = x, y = y, type = "deseq", transform = FALSE, nfolds = 5, repeats = 2, foldIdx = allFolds, tuneLength = 20)
pldaRes <- trainNBLDA(x = x, y = y, type = "deseq", tuneLength = 10, metric = "accuracy", train.control = ctrl2)
pred1 <- predict(pldaRes, xte, "predictions")

pldaRes@result@crossValidated$tuning.results[which(pldaRes@result@crossValidated$tuning.results[ ,1] == pldaRes@result@crossValidated$best.rho), ]
confusionMatrix(table(Predicted = pred1, Actual = yte))$overall[1]  # PLDA

nbldaRes <- trainNBLDA(x = x, y = y, type = "deseq", tuneLength = 10, metric = "accuracy", train.control = ctrl)
pred2 <- predict(nbldaRes, xte, "predictions")

nbldaRes@result@crossValidated$tuning.results[which(nbldaRes@result@crossValidated$tuning.results[ ,1] == nbldaRes@result@crossValidated$best.rho), ]
confusionMatrix(table(Predicted = pred2, Actual = yte))$overall[1]   # NBLDA


#### Cervical example
library(NBLDA)

data(cervical)
head(cervical[ ,1:10])

set.seed(2128)
idx <- sample(1:ncol(cervical), ceiling(ncol(cervical)*0.70), FALSE)

cond <- factor(rep(c("N", "T"), c(29, 29)))
tr <- t(cervical[ ,idx] + 1)
ts <- t(cervical[ ,-idx] + 1)

ytr <- cond[idx]
yts <- cond[-idx]

fit.nblda <- trainNBLDA(x = tr, y = ytr, type = "deseq", tuneLength = 10, metric = "accuracy",
                        train.control = nbldaControl(folds = 5, repeats = 10, target = 0, transform = FALSE,
                                                     return.selected.features = TRUE, phi.epsilon = 0.10))

fit.plda <- trainNBLDA(x = tr, y = ytr, type = "deseq", tuneLength = 10, metric = "accuracy",
                       train.control = nbldaControl(folds = 5, repeats = 10, target = 0, transform = FALSE, truephi = 0,
                                                    return.selected.features = TRUE, phi.epsilon = 0.10,
                                                    foldIdx = fit.nblda@result@control$foldIdx))

fit.plda
fit.nblda

preds.plda <- predict(fit.plda, ts)
preds.nblda <- predict(fit.nblda, ts)

table(Actual = yts, Predicted = preds.plda)
table(Actual = yts, Predicted = preds.nblda)

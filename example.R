
set.seed(110418)
dat <- generateCountData(n = 100, p = 200, K = 2, param = 1, sdsignal = 0.5, DE = 0.2, allZero.rm = FALSE, tag.samples = TRUE)
x <- t(dat$x + 1)
xte <- t(dat$xte + 1)
y <- as.factor(dat$y) #factor(c("T", "N")[dat$y])
yte <- as.factor(dat$yte) #factor(c("T", "N")[dat$yte])

rpt = 2
nf = 5

set.seed(110418)
# Each element of the returned list is the indices of test samples in this fold.
allFolds <- lapply(1:rpt, function(x){
  tmp <- balanced.folds(y, nfolds = nf)
  names(tmp) <- paste("Fold.", 1:nf, sep = "")
  tmp
})
names(allFolds) <- paste("Repeat.", 1:rpt, sep = "")

#
# system.time({poicla <- PoiClaClu::Classify.cv(x, y, folds = allFolds$Repeat.1, transform = FALSE, type = "deseq")})[3]
# system.time({poicla2 <- trainPLDA(x, y, type = "deseq", transform = FALSE, foldIdx = allFolds)})[3]

# PLDA vs NBLDA
ctrl = control(folds = 5, repeats = 2, foldIdx = allFolds, phi.epsilon = 0, transform = FALSE,
               truephi = NULL, target = 0, normalize.target = TRUE, delta = NULL)

pldaRes <- trainPLDA(x = x, y = y, type = "deseq", transform = FALSE, nfolds = 5, repeats = 2, foldIdx = allFolds, tuneLength = 20)
pred1 <- predict.plda(pldaRes, xte, "predictions")

pldaRes$tuning.results[which(pldaRes$tuning.results[ ,1] == pldaRes$best.rho), ]
confusionMatrix(table(Predicted = pred1, Actual = yte))$overall[1]  # PLDA

nbldaRes <- trainNBLDA(x = x, y = y, type = "deseq", tuneLength = 20, metric = "accuracy", train.control = ctrl)
pred2 <- predict.nblda(nbldaRes, xte, "predictions")

nbldaRes$tuning.results[which(nbldaRes$tuning.results[ ,1] == nbldaRes$best.rho), ]
confusionMatrix(table(Predicted = pred2, Actual = yte))$overall[1]   # NBLDA

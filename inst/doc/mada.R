### R code from vignette source 'mada.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: mada.Rnw:92-93
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mada.Rnw:105-106 (eval = FALSE)
###################################################
## install.packages("mada")


###################################################
### code chunk number 3: mada.Rnw:109-110
###################################################
library("mada")


###################################################
### code chunk number 4: mada.Rnw:142-144
###################################################
y <- 142 * .944 
y


###################################################
### code chunk number 5: mada.Rnw:147-148
###################################################
round(y)


###################################################
### code chunk number 6: mada.Rnw:153-158
###################################################
AuditC6 <- data.frame(TP = c(47, 126, 19, 36, 130, 84),
                      FN = c(9, 51, 10, 3, 19, 2),
                      FP = c(101, 272, 12, 78, 211, 68),
                      TN = c(738, 1543, 192, 276, 959, 89))
AuditC6


###################################################
### code chunk number 7: mada.Rnw:161-163
###################################################
AuditC6$names <- c("Study 1", "Study 2", "Study 4",
                   "Study 4", "Study 5", "Study 6")


###################################################
### code chunk number 8: mada.Rnw:167-169
###################################################
data("AuditC")
tail(AuditC)


###################################################
### code chunk number 9: mada.Rnw:190-191 (eval = FALSE)
###################################################
## madad(AuditC)


###################################################
### code chunk number 10: mada.Rnw:221-222 (eval = FALSE)
###################################################
## madad(AuditC, level = 0.80)


###################################################
### code chunk number 11: mada.Rnw:225-227
###################################################
AuditC.d <- madad(AuditC)
AuditC.d$fpr


###################################################
### code chunk number 12: mada.Rnw:243-245 (eval = FALSE)
###################################################
## forest(madad(AuditC), type = "sens")
## forest(madad(AuditC), type = "spec")


###################################################
### code chunk number 13: mada.Rnw:248-256
###################################################
pdf(file = "pairedforest.pdf", width = 12, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.3,0.3,0.3))
plot.new()
par(fig = c(0, 0.5, 0, 1), pty = "s", new = TRUE)
forest(AuditC.d, type = "sens", xlab = "Sensitivity")
par(fig = c(0.5, 1, 0, 1), pty = "s", new = TRUE)
forest(AuditC.d, type = "spec", xlab = "Specificity")
dev.off()


###################################################
### code chunk number 14: mada.Rnw:272-276
###################################################
rs <- rowSums(AuditC)
weights <- 4 * rs / max(rs)
crosshair(AuditC, xlim = c(0,0.6), ylim = c(0.4,1), 
          col = 1:14, lwd = weights)


###################################################
### code chunk number 15: mada.Rnw:280-282
###################################################
ROCellipse(AuditC, pch = "")
points(fpr(AuditC), sens(AuditC))


###################################################
### code chunk number 16: mada.Rnw:285-294
###################################################
pdf(file = "diagplots.pdf", width = 12, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot.new()
par(fig = c(0, 0.5, 0, 1), pty = "s", new = TRUE)
crosshair(AuditC, xlim = c(0,0.6), ylim = c(0.4,1), col = 1:14, lwd = weights)
par(fig = c(0.5, 1, 0, 1), pty = "s", new = TRUE)
ROCellipse(AuditC, pch = "")
points(fpr(AuditC), sens(AuditC))
dev.off()


###################################################
### code chunk number 17: mada.Rnw:341-343
###################################################
(fit.DOR.DSL <- madauni(AuditC))
(fit.DOR.MH <- madauni(AuditC, method = "MH"))


###################################################
### code chunk number 18: mada.Rnw:346-347
###################################################
summary(fit.DOR.DSL)


###################################################
### code chunk number 19: mada.Rnw:350-351
###################################################
forest(fit.DOR.DSL)


###################################################
### code chunk number 20: mada.Rnw:354-358
###################################################
pdf(file = "DORforest.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
forest(fit.DOR.DSL)
dev.off()


###################################################
### code chunk number 21: mada.Rnw:375-377
###################################################
(fit.phm.homo <- phm(AuditC, hetero = FALSE))
(fit.phm.het <- phm(AuditC))


###################################################
### code chunk number 22: mada.Rnw:380-381
###################################################
summary(fit.phm.homo)


###################################################
### code chunk number 23: mada.Rnw:384-385
###################################################
summary(fit.phm.het)


###################################################
### code chunk number 24: mada.Rnw:388-390
###################################################
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(AuditC, add = TRUE)


###################################################
### code chunk number 25: mada.Rnw:393-398
###################################################
pdf(file = "phmplot.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(AuditC, add = TRUE)
dev.off()


###################################################
### code chunk number 26: mada.Rnw:444-445
###################################################
(fit.reitsma <- reitsma(AuditC))


###################################################
### code chunk number 27: mada.Rnw:448-449
###################################################
summary(fit.reitsma)


###################################################
### code chunk number 28: mada.Rnw:452-457
###################################################
plot(fit.reitsma, sroclwd = 2,
     main = "SROC curve (bivariate model) for AUDIT-C data")
points(fpr(AuditC), sens(AuditC), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))


###################################################
### code chunk number 29: mada.Rnw:460-468
###################################################
pdf(file = "SROCAuditC.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.reitsma, sroclwd = 2,
     main = "SROC curve (bivariate model) for AUDIT-C data")
points(fpr(AuditC), sens(AuditC), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))
dev.off()


###################################################
### code chunk number 30: mada.Rnw:482-486
###################################################
data("IAQ")
data("SAQ")
fit.IAQ <- reitsma(IAQ)
fit.SAQ <- reitsma(SAQ)


###################################################
### code chunk number 31: mada.Rnw:489-496
###################################################
plot(fit.IAQ, xlim = c(0,.5), ylim = c(.5,1),
     main = "Comparison of IAQ and SAQ")
lines(sroc(fit.SAQ), lty = 2)
ROCellipse(fit.SAQ, lty = 2, pch = 2, add = TRUE)
points(fpr(IAQ), sens(IAQ), cex = .5)
points(fpr(SAQ), sens(SAQ), pch = 2, cex = 0.5)
legend("bottomright", c("IAQ", "SAQ"), pch = 1:2, lty = 1:2)


###################################################
### code chunk number 32: mada.Rnw:498-508
###################################################
pdf(file = "SAQIAQ.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.IAQ, xlim = c(0,.5), ylim = c(.5,1),
     main = "Comparison of IAQ and SAQ")
lines(sroc(fit.SAQ), lty = 2)
ROCellipse(fit.SAQ, lty = 2, pch = 2, add = TRUE)
points(fpr(IAQ), sens(IAQ), cex = .5)
points(fpr(SAQ), sens(SAQ), pch = 2, cex = 0.5)
legend("bottomright", c("IAQ", "SAQ"), pch = 1:2, lty = 1:2)
dev.off()


###################################################
### code chunk number 33: mada.Rnw:521-522
###################################################
data("smoking")


###################################################
### code chunk number 34: mada.Rnw:525-526
###################################################
summary(smoking$type)


###################################################
### code chunk number 35: mada.Rnw:529-530
###################################################
fit.smoking <- reitsma(smoking, formula = cbind(tsens, tfpr) ~ type)


###################################################
### code chunk number 36: mada.Rnw:534-535
###################################################
summary(fit.smoking)


###################################################
### code chunk number 37: mada.Rnw:547-551
###################################################
fit.smoking1 <- reitsma(smoking)
fit.smoking2 <- reitsma(smoking, alphasens = 0, alphafpr = 2)
AIC(fit.smoking1)
AIC(fit.smoking2)



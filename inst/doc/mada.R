### R code from vignette source 'mada.Rnw'

###################################################
### code chunk number 1: mada.Rnw:97-98
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mada.Rnw:110-111 (eval = FALSE)
###################################################
## install.packages("mada")


###################################################
### code chunk number 3: mada.Rnw:114-115
###################################################
library("mada")


###################################################
### code chunk number 4: mada.Rnw:147-149
###################################################
y <- 142 * .944 
y


###################################################
### code chunk number 5: mada.Rnw:152-153
###################################################
round(y)


###################################################
### code chunk number 6: mada.Rnw:158-163
###################################################
AuditC6 <- data.frame(TP = c(47, 126, 19, 36, 130, 84),
                      FN = c(9, 51, 10, 3, 19, 2),
                      FP = c(101, 272, 12, 78, 211, 68),
                      TN = c(738, 1543, 192, 276, 959, 89))
AuditC6


###################################################
### code chunk number 7: mada.Rnw:166-168
###################################################
AuditC6$names <- c("Study 1", "Study 2", "Study 4",
                   "Study 4", "Study 5", "Study 6")


###################################################
### code chunk number 8: mada.Rnw:172-174
###################################################
data("AuditC")
tail(AuditC)


###################################################
### code chunk number 9: mada.Rnw:195-196 (eval = FALSE)
###################################################
## madad(AuditC)


###################################################
### code chunk number 10: mada.Rnw:226-227 (eval = FALSE)
###################################################
## madad(AuditC, level = 0.80)


###################################################
### code chunk number 11: mada.Rnw:230-232
###################################################
AuditC.d <- madad(AuditC)
AuditC.d$fpr


###################################################
### code chunk number 12: mada.Rnw:248-250 (eval = FALSE)
###################################################
## forest(madad(AuditC), type = "sens")
## forest(madad(AuditC), type = "spec")


###################################################
### code chunk number 13: mada.Rnw:253-261
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
### code chunk number 14: mada.Rnw:277-281
###################################################
rs <- rowSums(AuditC)
weights <- 4 * rs / max(rs)
crosshair(AuditC, xlim = c(0,0.6), ylim = c(0.4,1), 
          col = 1:14, lwd = weights)


###################################################
### code chunk number 15: mada.Rnw:285-287
###################################################
ROCellipse(AuditC, pch = "")
points(fpr(AuditC), sens(AuditC))


###################################################
### code chunk number 16: mada.Rnw:290-299
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
### code chunk number 17: mada.Rnw:346-348
###################################################
(fit.DOR.DSL <- madauni(AuditC))
(fit.DOR.MH <- madauni(AuditC, method = "MH"))


###################################################
### code chunk number 18: mada.Rnw:351-352
###################################################
summary(fit.DOR.DSL)


###################################################
### code chunk number 19: mada.Rnw:355-356
###################################################
forest(fit.DOR.DSL)


###################################################
### code chunk number 20: mada.Rnw:359-363
###################################################
pdf(file = "DORforest.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
forest(fit.DOR.DSL)
dev.off()


###################################################
### code chunk number 21: mada.Rnw:380-382
###################################################
(fit.phm.homo <- phm(AuditC, hetero = FALSE))
(fit.phm.het <- phm(AuditC))


###################################################
### code chunk number 22: mada.Rnw:385-386
###################################################
summary(fit.phm.homo)


###################################################
### code chunk number 23: mada.Rnw:389-390
###################################################
summary(fit.phm.het)


###################################################
### code chunk number 24: mada.Rnw:393-395
###################################################
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(AuditC, add = TRUE)


###################################################
### code chunk number 25: mada.Rnw:398-403
###################################################
pdf(file = "phmplot.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(AuditC, add = TRUE)
dev.off()


###################################################
### code chunk number 26: mada.Rnw:449-450
###################################################
(fit.reitsma <- reitsma(AuditC))


###################################################
### code chunk number 27: mada.Rnw:453-454
###################################################
summary(fit.reitsma)


###################################################
### code chunk number 28: mada.Rnw:459-464
###################################################
plot(fit.reitsma, sroclwd = 2,
     main = "SROC curve (bivariate model) for AUDIT-C data")
points(fpr(AuditC), sens(AuditC), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))


###################################################
### code chunk number 29: mada.Rnw:467-475
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
### code chunk number 30: mada.Rnw:489-498
###################################################
data("IAQ")
data("SAQ")
# both datasets contain more than one 2x2-table per study
# reduce (somewhat arbitrarily) to one row per study by
# using the first coded table only:
IAQ1 <- subset(IAQ, IAQ$result_id == 1)
SAQ1 <- subset(SAQ, SAQ$result_id == 1)
fit.IAQ <- reitsma(IAQ1)
fit.SAQ <- reitsma(SAQ1)


###################################################
### code chunk number 31: mada.Rnw:501-508
###################################################
plot(fit.IAQ, xlim = c(0,.5), ylim = c(.5,1),
     main = "Comparison of IAQ and SAQ")
lines(sroc(fit.SAQ), lty = 2)
ROCellipse(fit.SAQ, lty = 2, pch = 2, add = TRUE)
points(fpr(IAQ1), sens(IAQ1), cex = .5)
points(fpr(SAQ1), sens(SAQ1), pch = 2, cex = 0.5)
legend("bottomright", c("IAQ", "SAQ"), pch = 1:2, lty = 1:2)


###################################################
### code chunk number 32: mada.Rnw:510-520
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
### code chunk number 33: mada.Rnw:533-536
###################################################
data("smoking")
# again reduce to one result per study:
smoking1 <- subset(smoking, smoking$result_id == 1)


###################################################
### code chunk number 34: mada.Rnw:539-540
###################################################
summary(smoking1$type)


###################################################
### code chunk number 35: mada.Rnw:543-545
###################################################
fit.smoking.type <- reitsma(smoking1, 
                            formula = cbind(tsens, tfpr) ~ type)


###################################################
### code chunk number 36: mada.Rnw:549-550
###################################################
summary(fit.smoking.type)


###################################################
### code chunk number 37: mada.Rnw:558-565
###################################################
fit.smoking.ml.type <- reitsma(smoking1, 
                          formula = cbind(tsens, tfpr) ~ type, 
                          method = "ml")
fit.smoking.ml.intercept <- reitsma(smoking1, 
                                    formula = cbind(tsens, tfpr) ~ 1,
                                    method = "ml")
anova(fit.smoking.ml.type, fit.smoking.ml.intercept)


###################################################
### code chunk number 38: mada.Rnw:575-581
###################################################
fit.smoking1 <- reitsma(smoking1, method = "ml")
fit.smoking2 <- reitsma(smoking1, 
                        alphasens = 0, alphafpr = 2, 
                        method = "ml")
AIC(fit.smoking1)
AIC(fit.smoking2)


###################################################
### code chunk number 39: mada.Rnw:590-592
###################################################
summary_pts_audit <- SummaryPts(reitsma(AuditC))
summary(summary_pts_audit)


###################################################
### code chunk number 40: mada.Rnw:600-602
###################################################
pred_audit1 <- predv_r(AuditC, prop_min=0.05, prop_max=0.15)
summary(pred_audit1)


###################################################
### code chunk number 41: mada.Rnw:609-611
###################################################
pred_audit2 <- predv_d(AuditC, prop_m=0.10, prop_sd=0.05)
summary(pred_audit2)




geneES <- sample(c(rnorm(5)-5, rnorm(40), rnorm(5)+5))
hist(geneES)


b <- boxplot(geneES)
negES <- match(b$out[b$out<0], geneES)
posES <- match(b$out[b$out>0], geneES)
neutralES <- match(setdiff(geneES, b$out), geneES)

resample_idx <- c(rep(negES, round(length(neutralES) / length(negES))),
                  rep(posES, round(length(neutralES) / length(posES))),
                  neutralES)

plot(density(geneES))
plot(density(geneES[resample_idx]))

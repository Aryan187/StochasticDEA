library(VineCopula)
library(copula)

tau <- matrix(0,12,12)

for (i in 1:12){
  if (i == 10)
    next
  for (j in 1:12){
    if (j == 10)
      next
    x <- pobs(dmu20)[,i]
    y <- pobs(dmu20)[,j]
    selectedCopula <- BiCopSelect(x,y,familyset=NA)
    tau[i,j] = selectedCopula$tau
  }
}

for (i in 1:12)
  tau[i,i] = 1

write.csv(tau,"D://Sem8//BTP//Data//tau20.csv",row.names=FALSE)
## binom
size=423
prob=1e-3
n <- 1e7

plot(dbinom(0:423,size=12,prob=.2)) # density
(0:10)[which(pbinom(0:10,size=size,prob=prob,lower.tail=FALSE) < 1e-3)]
qt <- pbinom(0:10,size=size,prob=prob,lower.tail=FALSE)
names(qt) <- 0:10
plot(0:10,qt,log="y")

pbinom(c(2,5),size=size,prob=prob,lower.tail=FALSE) # 9.17e-3,5.37e-6

binom.test(5,423,p=1e-3,alternative="greater")
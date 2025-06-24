#This is just an example of how the D-statist results obtained with angsd were plotted.
#All plots presented in the paper were produced in a similar manner.
#Below you will find the example for Figure 5a

####### Tagus #######

#Plot Dstats TAGUS all together

#TREE A - various pops S. pyrenaicus (P1), Pego (P2), Ardila (P3)
pyr_pego_ardila<-read.table(file="tagus_pego_ardila.txt", header = T)
pops<-pyr_pego_ardila[,1]

Dstat_pyr_pego_ardila_outg_torgal<- pyr_pego_ardila[,2]
Dstat_pyr_pego_ardila_outg_arade<- pyr_pego_ardila[,4]

sd_pyr_pego_ardila_outg_torgal<-pyr_pego_ardila[,3]
sd_pyr_pego_ardila_outg_arade<-pyr_pego_ardila[,5]


#TREE B - various pops S. pyrenaicus (P1), Ardila (P2), Pego (P3)
pyr_ardila_pego<-read.table(file="tagus_ardila_pego.txt", header = T)
#pops<-pyr_ardila_pego[,1]

Dstat_pyr_ardila_pego_outg_torgal<- pyr_ardila_pego[,2]
Dstat_pyr_ardila_pego_outg_arade<- pyr_ardila_pego[,4]

sd_pyr_ardila_pego_outg_torgal<-pyr_ardila_pego[,3]
sd_pyr_ardila_pego_outg_arade<-pyr_ardila_pego[,5]


###plot for both trees using S. torgalensis as outgroup
pdf("Dstat_spyrTagus_ardila_pego_both_trees.pdf")
plot(x=c(1,2,3,4,5,6,7),y=Dstat_pyr_pego_ardila_outg_torgal, ylim=c(-0.6, 0.2), xaxt='n', xlab='', ylab="D", pch=1, col="deepskyblue", xlim=c(0.9,7.2), cex=3, lwd=3, cex.axis=1.5, cex.lab=1.5)
points(x=c(1,2,3,4,5,6,7),y=Dstat_pyr_ardila_pego_outg_torgal, pch=16, col="deepskyblue", cex=3, lwd=3, cex.axis=1.5, cex.lab=1.5)
axis(1, at=1:7, labels = pops)
abline(h=0, lty=2)
dev.off()

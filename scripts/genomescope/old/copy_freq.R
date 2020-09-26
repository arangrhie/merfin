path<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

histfile <- "mBalMus1.k21.hist"
kmer_hist <- read.csv(file=histfile,sep="\t", header=FALSE)

x = kmer_hist[[1]]
y = kmer_hist[[2]]

model_sum=summary(model_4peaks[[1]])

kcovfloor = floor(min_max(model_sum$coefficients['kmercov',])[[1]])
het  = min_max(model_sum$coefficients['r',])
dups = min_max(model_sum$coefficients['bias',])
kcov = min_max(model_sum$coefficients['kmercov',])
mlen = min_max(model_sum$coefficients['length',])
md   = min_max(model_sum$coefficients['d',])

amlen = (mlen[1] + mlen[2]) / 2
ahet  = (het[1]  + het[2])  / 2
amd   = (md[1]   + md[2])   / 2
akcov = (kcov[1] + kcov[2]) / 2
adups = (dups[1] + dups[2]) / 2

unique_hist <- (2 * (1 - amd) * (1 - (1 - ahet)^k))                         * dnbinom(x, size = akcov     / adups, mu = akcov)     * amlen +
  ((amd * (1 - (1 - ahet)^k)^2) + (1 - amd) * ((1 - ahet)^k))  * dnbinom(x, size = akcov * 2 / adups, mu = akcov * 2) * amlen 

one_hist <- ((2*(1-amd)*(1-(1-ahet)^k)) + (2*amd*(1-(1-ahet)^k)^2) + (2*amd*((1-ahet)^k)*(1-(1-ahet)^k))) * dnbinom(x, size = akcov   / adups, mu = akcov)   
two_hist <- (((1-amd)*((1-ahet)^k)) + (amd*(1-(1-ahet)^k)^2))                                             * dnbinom(x, size = akcov*2 / adups, mu = akcov * 2)
thr_hist <- (2*amd*((1-ahet)^k)*(1-(1-ahet)^k))                                                           * dnbinom(x, size = akcov*3 / adups, mu = akcov * 3)
fou_hist <- (amd*(1-ahet)^(2*k))                                                                          * dnbinom(x, size = akcov*4 / adups, mu = akcov * 4)

one_hist[(akcov*5):length(one_hist)] <- 0
two_hist[(akcov*5):length(two_hist)] <- 0
thr_hist[(akcov*5):length(thr_hist)] <- 0
fou_hist[(akcov*5):length(fou_hist)] <- 0

pred=predict(model_4peaks[[1]], newdata=data.frame(x))

## Compute error rate, by counting kmers unexplained by model through first peak
## truncate errors as soon as it goes to zero, dont allow it to go back up
error_xcutoff = kcovfloor
error_xcutoff_ind = which(x==error_xcutoff)

error_kmers = y[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind]
error_kmers = pmax(error_kmers, 1e-10)

repeat_kmers = y[1:length(y)] - pred[1:length(y)]
repeat_kmers <- c(rep(0,akcov*5-1),repeat_kmers[ceiling(akcov*5)[[1]]:length(repeat_kmers)])

total_kmers = sum(as.numeric(x*y))
unique_kmers = sum(as.numeric(x*unique_hist))
total_error_kmers = sum(as.numeric(error_kmers * x[1:error_xcutoff_ind]))
total_repeat_kmers = total_kmers - unique_kmers - total_error_kmers
repeat_len=repeat_kmers/(2*kcov)

repeats <- 0

# for (i in 1:30){
# 
#   frc <- repeat_kmers[akcov*(5+i)]/total_repeat_kmers*5000                   * dnbinom(x, size = akcov * (5+i) / adups, mu = akcov * (5+i))
# 
#   repeats <- repeats+frc
#   
# }
# 
# full <- (one_hist + two_hist + thr_hist + fou_hist + repeats) * amlen
# 
# dev.off()
# plot(kmer_hist,type="n", main="GenomeScope Profile\n", xlab="Coverage", ylab="Frequency", ylim=c(1,50000000),log='xy')
# 
# points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
# lines(error_kmers, col="black")
# lines(full, col="gray")

lookup_table <- NULL
plot_table <- NULL

for (i in 1:(akcov*5-1)){
  
  totalP<-sum(na.zero(error_kmers[i]/total_kmers),
              na.zero(one_hist[i]),
              na.zero(two_hist[i]),
              na.zero(thr_hist[i]),
              na.zero(fou_hist[i]))
  
  plot_table<-rbind(plot_table,c(na.zero(error_kmers[i]/total_kmers/totalP),
       na.zero(one_hist[i]/totalP),
       na.zero(two_hist[i]/totalP),
       na.zero(thr_hist[i]/totalP),
       na.zero(fou_hist[i]/totalP)
       ))

  p.values<-c(na.zero(error_kmers[i]/total_kmers/totalP),
      na.zero(one_hist[i]/totalP),
      na.zero(two_hist[i]/totalP),
      na.zero(thr_hist[i]/totalP),
      na.zero(fou_hist[i]/totalP))
  
  max.p<-max(p.values)  
  readK<-which.max(p.values)-1
  
  lookup_table<-rbind(lookup_table,c(readK,max.p))
  
}

plot_table<-data.frame(plot_table)
lookup_table<-data.frame(lookup_table)

write.table(lookup_table, "lookup_table.txt", sep=",", row.names = FALSE, col.names = FALSE)

colors <- c("black","red","green","purple","blue","orange","gray","yellow","violet")

dev.off()

plot(lookup_table$X1,type="n", xlab="Coverage", ylab="Probability",xlim=c(0,200), ylim=c(0,1))

peaks<-akcov * 1:60000

abline(v=peaks, col="black", lty=2, lwd=0.1)

for (i in 1:5) {
  lines(rownames(plot_table), plot_table[,i], type="l", lwd=1.5,
        col=colors[i])
}

#homozygous diploid case d=1, r=0
(d*(1-r)^(2*k))  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length

#haploid case d=0, r=0
((1-d)*((1-r)^k)) * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length


#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = T) 
a = read.csv(args[1],sep='\t')

pdf(args[2])

a$edfill = a$edfill/1000000.0
a$diagfill = a$diagfill/1000000.0
a$difffill = a$difffill/1000000.0
a$gabafill = a$gabafill/1000000.0

a$edtrace = a$edtrace/1000000.0
a$diagtrace = a$diagtrace/1000000.0
a$difftrace = a$difftrace/1000000.0
a$gabatrace = a$gabatrace/1000000.0

a$edconv = a$edconv/1000000.0
a$diagconv = a$diagconv/1000000.0
a$diffconv = a$diffconv/1000000.0
a$gabaconv = a$gabaconv/1000000.0

a$edtotal = a$edtotal/1000000.0
a$diagtotal = a$diagtotal/1000000.0
a$difftotal = a$difftotal/1000000.0
a$gabatotal = a$gabatotal/1000000.0

regfill = a$len*0.015
regtrace = a$len*0.004
regothers = a$len*0.005

par(pin = c(3,3))

plot(a$len,regfill,col='lightgrey',type='n',pch=20,log='xy',xlim=c(0.9,25000),ylim=c(0.5,500),xlab='Query sequence length (bp)',ylab='Calculation time (us)')
lines(a$len,regfill,col='lightgrey',lty=1)

points(a$len,a$edfill,col='black',type='n',pch=20)
points(a$len,a$diagfill,col='black',type='n',pch=20)
points(a$len,a$difffill,col='black',type='n',pch=20)
points(a$len,a$gabafill,col='black',type='n',pch=20)
lines(a$len,a$edfill,col='black',lty=3)
lines(a$len,a$diagfill,col='black',lty=4)
lines(a$len,a$difffill,col='black',lty=5)
lines(a$len,a$gabafill,col='black',lty=1)

legend('bottomright',legend=c('ED', 'non-diff', 'diff-raw', 'libgaba'),
	col=c('black','black','black','black'),lty=c(3,4,5,1))


plot(a$len,regtrace,col='lightgrey',type='n',pch=20,log='xy',xlim=c(0.9,25000),ylim=c(0.5,500),xlab='Query sequence length (bp)',ylab='Calculation time (us)')
lines(a$len,regtrace,col='lightgrey',lty=1)

points(a$len,a$edtrace,col='black',type='n',pch=20)
points(a$len,a$diagtrace,col='black',type='n',pch=20)
points(a$len,a$difftrace,col='black',type='n',pch=20)
points(a$len,a$gabatrace,col='black',type='n',pch=20)
lines(a$len,a$edtrace,col='black',lty=3)
lines(a$len,a$diagtrace,col='black',lty=4)
lines(a$len,a$difftrace,col='black',lty=5)
lines(a$len,a$gabatrace,col='black',lty=1)

legend('bottomright',legend=c('ED', 'non-diff', 'diff-raw', 'libgaba'),
	col=c('black','black','black','black'),lty=c(3,4,5,1))


plot(a$len,regothers,col='lightgrey',type='n',pch=20,log='xy',xlim=c(0.9,25000),ylim=c(0.5,500),xlab='Query sequence length (bp)',ylab='Calculation time (us)')
lines(a$len,regothers,col='lightgrey',lty=1)

points(a$len,a$edtrace+a$edconv,col='black',type='n',pch=20)
points(a$len,a$diagtrace+a$diagconv,col='black',type='n',pch=20)
points(a$len,a$difftrace+a$diffconv,col='black',type='n',pch=20)
points(a$len,a$gabatrace+a$gabaconv,col='black',type='n',pch=20)
lines(a$len,a$edtrace+a$edconv,col='black',lty=3)
lines(a$len,a$diagtrace+a$diagconv,col='black',lty=4)
lines(a$len,a$difftrace+a$diffconv,col='black',lty=5)
lines(a$len,a$gabatrace+a$gabaconv,col='black',lty=1)

legend('bottomright',legend=c('ED', 'non-diff', 'diff-raw', 'libgaba'),
	col=c('black','black','black','black'),lty=c(3,4,5,1))
dev.off()


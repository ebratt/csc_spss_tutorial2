#######################################################
# PROBLEM 1                                           #
#######################################################
# setup
# clear the environment
rm(list=ls())

DATA_DIR <- './data'
IMAGES_DIR <- './images'
OUTPUT_DIR <- './output'

make_dir <- function(d) {
    if (file.exists(d)) unlink(d, recursive=TRUE, force=TRUE)
    dir.create(d)
}
lapply(c(IMAGES_DIR, OUTPUT_DIR),make_dir)


## function that concatenates strings (useful for directory paths)
concat <- function(x1,x2) {
    result <- paste(x1,x2,sep="")
    return(result)
}

## function that checks to see if a package is installed and,if not,installs it
## portions of this code came from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
load_package <- function(x) {
    if (x %in% rownames(installed.packages())) { 
        print(concat("package already installed: ", x))
    }
    else { 
        install.packages(x, repos="http://mirror.las.iastate.edu/CRAN/") 
    }
    library(x, character.only=TRUE)
}

####################################################
# problem1.sav

load_package("foreign")
data <- read.spss(concat(DATA_DIR,'/problem1.sav'), to.data.frame=TRUE)
dim(data)
str(data)
head(data)
tail(data)
#****
# convert factors to numbers
summary(data$facsex)
summary(data$race)

# perform PCA
fit <- prcomp(data, center=TRUE, scale=TRUE)
summary(fit)
# 5 principal components are required to explain 90% of the total variation

# variance charts
# from http://rstudio-pubs-static.s3.amazonaws.com/27823_dbc155ba66444eae9eb0a6bacb36824f.html
pcaCharts <- function(x) {
    x.var <- x$sdev ^ 2
    x.pvar <- x.var/sum(x.var)
    print("proportions of variance:")
    print(x.pvar)
    
    par(mfrow=c(2,2))
    plot(x.pvar,xlab="Principal component", ylab="Proportion of variance explained", ylim=c(0,1), type='b')
    plot(cumsum(x.pvar),xlab="Principal component", ylab="Cumulative Proportion of variance explained", ylim=c(0,1), type='b')
    screeplot(x)
    screeplot(x,type="l")
    par(mfrow=c(1,1))
}
png(concat(IMAGES_DIR,'/employment - variance graphs.png'), 
    width = 512, height = 512)
pcaCharts(fit)
dev.off()

# eigenvalues
(eigenvalues <- fit$sdev^2)
write.table(round(eigenvalues[1:5], 2), file=concat(OUTPUT_DIR,'/employment - eigenvalues.csv'), sep=",")
png(concat(IMAGES_DIR,'/employment - eigenvalues.png'), 
    width = 512, height = 512)
plot(eigenvalues, type='b', xlab="Principal Components", ylab="Eigenvalues")
dev.off()

# eigenvectors
cov_matrix <- cov(data[2:10])
eigs <- eigen(cov_matrix)
(eigenvalues <- round(eigs$values, 4))
(eigenvectors <- round(eigs$vectors, 4))
write.table(eigenvectors, file=concat(OUTPUT_DIR,'/employment - eigenvectors.csv'), sep=",")

# rotations
(rotations <- round(fit$rotation[,1:5], 4))
write.table(round(fit$rotation[,1:5], 4), file=concat(OUTPUT_DIR,'/employment - rotations.csv'), sep=",")

# bi-directional plot
load_package('ggplot2')
countries <- data[,1]
PCbiplot <- function(PC, rownames, x="PC1", y="PC2") {
    # code is modified but mostly borrowed from:
    #    http://stackoverflow.com/questions/6578355/plotting-pca-biplot-with-ggplot2
    #    posted by http://stackoverflow.com/users/577462/crayola
    # PC being a prcomp object
    data <- data.frame(obsnames=countries, PC$x)
    plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
    plot <- plot + geom_hline(alpha=0.4, size=.2, yintercept=0) + geom_vline(alpha=0.4, size=.2, xintercept=0)
    datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
    mult <- min(
        (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
        (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
    )
    datapc <- transform(datapc,
                        v1 = .7 * mult * (get(x)),
                        v2 = .7 * mult * (get(y))
    )
    plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
    plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
    plot
}
png(concat(IMAGES_DIR,'/employment - biplot.png'), 
    width = 1024, height = 1024)
PCbiplot(fit, countries)
dev.off()

# individual coordinates
write.table(fit$x, file=concat(OUTPUT_DIR,'/employment - coordinates.csv'), sep=",")

load_package("MVA")
png(concat(IMAGES_DIR,'/employment - biplot matrix.png'), width = 2048, height = 2048)
pairs(fit$x[,1:4], ylim=c(-6,6),xlim=c(-6,6),panel=function(x,y,...) {
    text(x,y,abbreviate(countries), cex=0.5)
    bvbox(cbind(x,y), add=TRUE, cex=0.6)
})
dev.off()

# New Data; from each variable to the PC's
# Agriculture
par(mfrow=c(2,2))
out <- sapply(1:4, function(i) {
    plot(data$Agr, fit$x[,i],
         xlab=paste("PC",i,sep=""),
         ylab="Percentage employed in agriculture")
})

# Mining
par(mfrow=c(2,2))
out <- sapply(1:4, function(i) {
    plot(data$Min, fit$x[,i],
         xlab=paste("PC",i,sep=""),
         ylab="Percentage employed in mining")
})

# Manufacturing
par(mfrow=c(2,2))
out <- sapply(1:4, function(i) {
    plot(data$Man, fit$x[,i],
         xlab=paste("PC",i,sep=""),
         ylab="Percentage employed in manufacturing")
})

# Power supply
par(mfrow=c(2,2))
out <- sapply(1:4, function(i) {
    plot(data$PS, fit$x[,i],
         xlab=paste("PC",i,sep=""),
         ylab="Percentage employed in power supply")
})
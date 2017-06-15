confidenceplots <- function(x,y,group,varnames=c("x","y"),groupnames=sort(unique(group)), groupcols=rainbow(length(unique(group))), groupshadecols=rep("gray80",length(unique(group))), shownames=T, xlim=c(NA,NA), ylim=c(NA,NA), lty=1, lwd=1, add=F, alpha=0.95, ellipse=0, shade=F, frac=0.01, cex=1) 
{
  #
  # function to plot x-y scatterplot of groups means and 95% 
  # confidence ellipses  
  # where group classifies the data into groups
  # 
  # package "ellipse" required
  #
  # x is the x-variable
  # y is the y-variable
  # group is the grouping variable
  # varnames is a vector of two labels for the axes (default x and y)
  # groupnames is a vector of labels for the groups (default is 1, 2, etc...)
  # xlim and ylim are possible new limits for the plot
  # lwd is the line width for the ellipse (default is 1)
  # lty is the line type for the ellipse (default is 1)
  # set ellipse<0 for regular data-covering ellipses
  # set ellipse=0 (default) for normal-theory confidence ellipses
  # set ellipse=1 for bootstrap confidence ellipses
  # set ellipse=2 for confidence error bars along ellipse axes [not implemented yet!] 
  # set ellipse=3 for normal-theory confidence error bars lined up with axes
  # set ellipse=4 for bootstrap confidence error bars along axes

  # first find min and max of all ellipses

  groups <- sort(unique(group))
  xmin <- mean(x)
  xmax <- mean(x)
  ymin <- mean(y)
  ymax <- mean(y)
  for(k in groups) {
    points <- cbind(x,y)[group==k,]
    npts <- nrow(points)
    if(!is.matrix(points)) npts <- 1
    if(npts>=2) {
      covpoints <- var(points)
      meanpoints <- as.numeric(apply(points, 2, mean))
      rconf  <- sqrt(2 * (npts-1) * qf(alpha, 2, npts-2)/(npts*(npts-2)))
      conf.elip <- ellipse(covpoints/npts, centre=meanpoints, level=alpha)
      if(ellipse<0) conf.elip <- ellipse(covpoints, center=meanpoints, level=alpha)
      xmin <- min(xmin, conf.elip[,1])
      xmax <- max(xmax, conf.elip[,1])
      ymin <- min(ymin, conf.elip[,2])
      ymax <- max(ymax, conf.elip[,2])
    }
  }
  par(mar=c(4.2,4.2,1,1), cex.axis=0.8)
  if(is.na(xlim[1])) xlim=c(xmin,xmax)
  if(is.na(ylim[1])) ylim=c(ymin,ymax)
  if(!add) plot(x, y, type="n", xlab=varnames[1], ylab=varnames[2], main="", xlim=xlim, ylim=ylim)
  groups.col <- groupcols
  groups.shade.col <- groupshadecols
  index <- 0

# -----------------------------------------------
# ellipse negative: regular data covering regions

  if(ellipse<0) {
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        covpoints <- var(points)
        meanpoints <- c(mean(points[,1]), mean(points[,2]))
        rconf  <- sqrt(2 * (npts-1) * qf(alpha, 2, npts-2)/(npts*(npts-2)))
        conf.elip <- ellipse(covpoints, centre=meanpoints, level=alpha)
        if(!shade) lines(conf.elip, col=groups.col[index], lty=lty, lwd=lwd)
        if(shade) polygon(conf.elip, col=groups.shade.col[index], border=NA)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
      }
    }
  }
# --------------------------------------------
# ellipse = 0 normal theory confidence regions

  if(ellipse==0) {	
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        covpoints <- var(points)
        meanpoints <- c(mean(points[,1]), mean(points[,2]))
        rconf  <- sqrt(2 * (npts-1) * qf(alpha, 2, npts-2)/(npts*(npts-2)))
        conf.elip <- ellipse(covpoints/npts, centre=meanpoints, level=alpha)
        if(!shade) lines(conf.elip, col=groups.col[index], lty=lty, lwd=lwd)
        if(shade) polygon(conf.elip, col=adjustcolor(groups.col[index], alpha.f=0.2), border=NA)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
      }
    }
  }

# -----------------------------------------
# ellipse = 1: bootstrap confidence regions

  if(ellipse==1) {
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        boots <- matrix(0, nrow=1000, ncol=2)
        for(iboot in 1:1000) boots[iboot,] <- apply(points[sample(1:npts, replace=T),],2,mean)
        covpoints <- var(boots)
        meanpoints <- c(mean(points[,1]), mean(points[,2]))
        print(meanpoints)
        conf.elip <- ellipse(covpoints, centre=meanpoints, level=alpha)
        if(!shade) lines(conf.elip, col=groups.col[index], lty=lty, lwd=lwd)
        if(shade) polygon(conf.elip, col=adjustcolor(groups.shade.col[index], alpha.f=0.2), border=NA)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
      }
    }
  }

# --------------------------------------------------
# ellipse = 2: delta method, non-operative at moment

# -----------------------------------------------------------------------
# ellipse = 3: confidence intervals for separate variables, normal theory

  if(ellipse==3) {
    xrange <- xlim[2]-xlim[1]
    yrange <- ylim[2]-ylim[1]
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        sdpoints <- apply(points, 2, sd)
        meanpoints <- apply(points, 2, mean)
        lines(c(meanpoints[1]-qt(alpha, npts-1)*sdpoints[1]/sqrt(npts),meanpoints[1]+qt(alpha, npts-1)*sdpoints[1]/sqrt(npts)), c(meanpoints[2],meanpoints[2]), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1],meanpoints[1]), c(meanpoints[2]-qt(alpha, npts-1)*sdpoints[2]/sqrt(npts),meanpoints[2]+qt(alpha, npts-1)*sdpoints[2]/sqrt(npts)), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-qt(alpha, npts-1)*sdpoints[1]/sqrt(npts),meanpoints[1]-qt(alpha, npts-1)*sdpoints[1]/sqrt(npts)), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]+qt(alpha, npts-1)*sdpoints[1]/sqrt(npts),meanpoints[1]+qt(alpha, npts-1)*sdpoints[1]/sqrt(npts)), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(meanpoints[2]-qt(alpha, npts-1)*sdpoints[2]/sqrt(npts),meanpoints[2]-qt(alpha, npts-1)*sdpoints[2]/sqrt(npts)), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(meanpoints[2]+qt(alpha, npts-1)*sdpoints[2]/sqrt(npts),meanpoints[2]+qt(alpha, npts-1)*sdpoints[2]/sqrt(npts)), col="gray", lwd=lwd, lty=lty)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
      }
    }
  }

# --------------------------------------------------------------------------------
# ellipse = 4: bootstrap confidence intervals for separate variables, by bootstrap

  if(ellipse==4) {
    xrange <- xlim[2]-xlim[1]
    yrange <- ylim[2]-ylim[1]
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        boots <- matrix(0, nrow=1000, ncol=2)
        for(iboot in 1:1000) boots[iboot,] <- apply(points[sample(1:npts, replace=T),],2,mean)
        xquant <- quantile(boots[,1], c((1-alpha)/2,(1+alpha)/2))
        yquant <- quantile(boots[,2], c((1-alpha)/2,(1+alpha)/2))
        meanpoints <- apply(points, 2, mean)
        lines(c(xquant[1],xquant[2]), c(meanpoints[2],meanpoints[2]), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1],meanpoints[1]), c(yquant[1], yquant[2]), col="gray", lwd=lwd, lty=lty)
        lines(c(xquant[1],xquant[1]), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(xquant[2],xquant[2]), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(yquant[1],yquant[1]), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(yquant[2],yquant[2]), col="gray", lwd=lwd, lty=lty)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
      }
    }
  }
}

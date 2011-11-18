#!/usr/bin/r

# this function is calling from BatchAnalaysis.sh
# EX)
# bash R --slave --args [list of argument here ]< $[this R file name]

#
# list file includes:

myArg                <- commandArgs();
dtFilename           <- as.character( myArg[4] ); # first argument excluding R slave
outputFilenamePrefix <- as.character( myArg[5] ); # second argument excluding R slave 

# ---------------------------------------------------------------------------- #
plotLine <- function( dt, lineColor)
{
  points( dt$iteration, dt$SSE, 
          col=lineColor,
          type="l", lty=2, lwd=1);
}

# ---------------------------------------------------------------------------- #
myColorList <- 1:50;
dt<-read.csv( dtFilename, header=F );
colnames( dt )<- c("structure","dummy", "subjectID", "dummy2", "iteration", "dummy3", "SSE");

# get range of graph
range.min <- min( dt$SSE)
range.max <- max( dt$SSE)
range.length <- range.max - range.min;

subjectList   <- levels(factor( dt$subjectID))
structureList <- levels( factor( dt$structure ) );


for( currentStructure in structureList)
{
  colorIndex<-1;
  currentDT <- subset( dt, dt$structure == currentStructure );
  # filename
  currentPlotFilename <-  paste( outputFilenamePrefix, currentStructure, ".pdf",sep="");

  # get mean of error for each iteration
  currentMeanSSE<-aggregate( currentDT$SSE, list(iteration=currentDT$iteration), mean);


  pdf( currentPlotFilename );

  # expand right side of clipping rect to make room for the legend
  par( xpd=T, mar=par()$mar+c(0,0,0,4) );

  plot( currentMeanSSE$iteration, currentMeanSSE$x,
        type="l",lty=1, lwd=2, col="darkblue",
        ylim=c( range.min , range.max ),
        xlab="iteration",
        ylab="SSE"
      );

  # min mean location
  minMeanSSE <- min( currentMeanSSE$x );
  minMeanSSEIteration <-currentMeanSSE$iteration[which.min( currentMeanSSE$x )];

  abline(v=minMeanSSEIteration, col="black", lwd=2, lty=2 );
  mtext( paste(" Min SSE=", round(minMeanSSE,4)," at ", minMeanSSEIteration ), 
         side=3, # at the top
         line=-1
      );

    
  # plot dotted line of individual subjects
  for( currentSubject in subjectList )
  {
    currentSet<-subset(dt, dt$subjectID==currentSubject & dt$structure==currentStructure);
    plotLine( currentSet , myColorList[colorIndex] );
    colorIndex<- colorIndex+1;
  }
  title( currentStructure );
  legend(  max(dt$iteration)*1.05,range.max*1.1, c("mean",subjectList), 
           col=c("darkblue",myColorList), pch=20, lty=1, bty="n",
           cex=0.8);
  dev.off();
}


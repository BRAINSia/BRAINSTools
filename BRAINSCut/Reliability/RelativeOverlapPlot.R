#!/usr/bin/r

# this function is calling from BatchAnalaysis.sh
# EX)
# bash R --slave --args ListFIle < $SIPlotRScript 

#
# list file includes:
# subjectID, CSVThresholdFilename, ROIName

myArg             <- commandArgs();
listFIlename      <- as.character( myArg[4] ); # first argument excluding R slave 
plotName          <- as.character( myArg[5] ); # first argument excluding R slave 

dtList <- read.csv( listFIlename );

# ---------------------------------------------------------------------------- #
# plot similarity index graph for one subject
# ---------------------------------------------------------------------------- #
OneSubjectPlotOfSI <- function( subjectID, csvFilename, ROIName, myColor, cumulativeVolumes )
{
  csvFilename=as.character(csvFilename);
  print( csvFilename );
  dt<-read.csv( csvFilename );

  points( dt$threshold, dt$SI,            # x,y
          type="l", lty=2, col=myColor    #line specification
        )

  newColnames <- c( colnames( cumulativeVolumes), subjectID );

  cumulativeVolumes$subjectID <- dt$ann;
  colnames( cumulativeVolumes ) <- newColnames;

  cumulativeVolumes;
}

# ---------------------------------------------------------------------------- #
# plot Mean SI
# ---------------------------------------------------------------------------- #
MeanSIPlot <- function( data )
{
  meanSI <- rowMeans( data[ 2: ncol(data) ] );
  points( data$threshold, meanSI, 
          type="l" ,lwd="3", lty=4,col="darkblue" );

  maxMeanSI <-  max( meanSI );
  maxMeanThreshold <- which.max( meanSI );

  require( fields );
  xline( data$threshold[maxMeanThreshold], col="darkblue", lty=2, lwd=2 );
  mtext( paste( "Maximum Mean SI : ", round( maxMeanSI,2) , 
                " at ", maxMeanThreshold ) ,
         side=3, line= -1 );
  
}

# ---------------------------------------------------------------------------- #
# ICC 
# ---------------------------------------------------------------------------- #
ICCPlot <- function( man, dt)
{
  print( "dt")
  print( dt )
  print( "man")
  print( man)
  require(psy);

  print( ncol(dt)-1 );
  myICC      <- array( 0, c(2,ncol(dt)-1 ) );

  for( currentCol in 2:ncol( dt ) )
  {
    pairedData <-matrix(c( man, dt[,currentCol] ) , #data 
                        c(2, ncol(dt) ), byrow=T );                # row * col

    print( pairedData);
    iccResult <- icc( pairedData );

    myICC[1,currentCol-1] <- iccResult$icc.agreement;
    myICC[2,currentCol-1] <- iccResult$icc.consistency;
  }

  points( dt$threshold, myICC[1,] ,
          lwd=2, lty=1, type="l", col="darkblue" );
  points( dt$threshold, myICC[2,] ,
          lwd=2, lty=1, type="l", col="red" );
  print ( myICC );

  maxICC <- max.col( myICC );
  mtext( paste( "Maximum ICC  : ", 
                round(myICC[maxICC],3) , ",", 
                round(myICC[maxICC],3) ),
         side=3, line=-2 );
  
}
# ---------------------------------------------------------------------------- #
# plot Subjects
# ---------------------------------------------------------------------------- #

# plot first one
numberOfSubject = nrow( dtList );


currentDTFilename <- as.character(dtList$csvFIlename[1]);
currentSubjectID  <- as.character(dtList$subjectID[1]);
currentROIName    <- as.character(dtList$ROIName[1]);

currentDT <- read.csv( currentDTFilename, header=T );

# data for ICC
cumulativeVolumes <- subset( currentDT, select=c("threshold","ann") );
colnames( cumulativeVolumes )<- c("threshold", currentSubjectID);

cumulativeSI      <- subset( currentDT, select=c("threshold","SI") );

manualData     <- array(0, c(1, numberOfSubject));
manualData[1]  <- currentDT$manual[1];

pdf( plotName );
plot( currentDT$threshold, currentDT$SI,
      type="l", lty=2, col=1,    #line specification
      xlab="threshold", ylab="similarity index, ICC" ,
      xlim=c(0.0,1.0), ylim=c(0.0,1.0) );
title( currentROIName );
    
# plot rest of subjects

for( i in 2:numberOfSubject )
{
  currentDTFilename <- as.character(dtList$csvFIlename[i]);
  currentSubjectID  <- as.character(dtList$subjectID[i]);
  currentROIName    <- as.character(dtList$ROIName[i]);
  currentDT         <- read.csv( currentDTFilename, header=T );

  cumulativeVolumes=OneSubjectPlotOfSI( currentSubjectID, currentDTFilename, currentROIName, i, cumulativeVolumes);

  # gether SI
  tempColNames <- c( colnames( cumulativeSI ), currentSubjectID );
  cumulativeSI$TEMP <- currentDT$SI;
  colnames( cumulativeSI) <- tempColNames;

  # gether manuals
  manualData[i]  <- currentDT$manual[i];
}

# plot mean
MeanSIPlot( cumulativeSI );

# ICC accros threshold
ICCPlot( manualData, cumulativeVolumes );

dev.off();

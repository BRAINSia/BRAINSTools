#!/usr/bin/r

# this function is calling from BatchAnalaysis.sh
# EX)
# bash R --slave --args ListFIle < $SIPlotRScript 

#
# list file includes:
# subjectID, CSVThresholdFilename, ROIName

myArg             <- commandArgs();
listFilename      <- as.character( myArg[4] ); # first argument excluding R slave 
plotName          <- as.character( myArg[5] ); # first argument excluding R slave 

dtList <- read.csv( listFilename );

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
MeanSIPlot <- function( data  )
{
  meanSI <- rowMeans( data[ 2: ncol(data) ] );
  points( data$threshold, meanSI, 
          type="l" ,lwd="3", lty=4,col="darkblue" );

  SI <-  max( meanSI );
  SIThreshold <- data$threshold[ which.max( meanSI ) ];

  MaximumMeanText <- paste( "Maximum Mean SI : ", round( SI,2) ,
                            " at ", SIThreshold  );

  mtext( MaximumMeanText,
         side=1, line= -1 );

  # return summary
  c( SI, SIThreshold);
}

# ---------------------------------------------------------------------------- #
# ICC 
# ---------------------------------------------------------------------------- #
ICCPlot <- function( man, dt )
{

  NoSubject <- ncol(dt)-1;

  manData <- man[1];
  for( i in 2:NoSubject)
  {
    manData<-rbind( manData, man[i] );
  }

  require("psy");

  ICCTable      <- array( 0, c(2,nrow(dt) ) );

  for( currentThresholdRow in 2:nrow( dt ) )
  {
    annData<-dt[ currentThresholdRow,2];
    for( i in 2:NoSubject) 
    {
      annData<-rbind(annData,dt[ currentThresholdRow,i+1]) ; # first element is threshold
    }
    pairedData <- cbind( manData, annData);

    iccResult <- icc( pairedData );

    ICCTable[1,currentThresholdRow-1] <- iccResult$icc.agreement;
    ICCTable[2,currentThresholdRow-1] <- iccResult$icc.consistency;
  }

  points( dt$threshold, ICCTable[1,] ,
          lwd=2, lty=1, type="l", col="darkblue" );
  points( dt$threshold, ICCTable[2,] ,
          lwd=2, lty=1, type="l", col="red" );

  ICCmax <- max.col( ICCTable );

  ICCagreementColumn<-1; ICCconsistencyColumn<-2;

  ICCmaximumAgreementIndex<-ICCmax[1];
  ICCmaximumConsistencyIndex<-ICCmax[2];

  # rounding points
  roundPrecision<-3;

  # summary stats
  iccA <- ICCTable[ ICCagreementColumn, ICCmaximumAgreementIndex];
  iccAThreshold <- dt$threshold[ICCmaximumAgreementIndex] ;

  iccC <- ICCTable[ ICCconsistencyColumn, ICCmaximumConsistencyIndex];
  iccCThreshold <- dt$threshold[ICCmaximumConsistencyIndex];

  # Vertical line
  # require( fields );
  abline( v=iccAThreshold, col="red", lty=3, lwd=2);

  # description on the plot
  MaximumICCAText <- paste( "Maximum ICC(A)  : ", 
                             round( iccA, roundPrecision ) ,
                            " at ", iccAThreshold )
  mtext( MaximumICCAText, side=1, line=-2 );

  MaximumICCCText <- paste( "Maximum ICC(C)  : ", round( iccC , roundPrecision ),
                            " at ", iccCThreshold  );
  mtext( MaximumICCCText,side=1, line=-3 );

  legend( "topright", 
          c( "ICC(A)","ICC(C)" ),
          col=c("darkblue", "red"),
          lty=1, lwd=2,
          bty="n" )

# return threshold value at best agreement
  c(ICCmaximumAgreementIndex, iccA, iccAThreshold, iccC, iccCThreshold);
  
}
# ---------------------------------------------------------------------------- #
# final plot of volumetric measurement at the maxumum agreement 
# ---------------------------------------------------------------------------- #
VolumetricComparisonPlot <- function( manual, ann, threshold)
{
  NoSubject<-ncol(ann)-1;

  # first colume is subject number
  firstSubjectRow<-2;
  annVolume<-ann[ threshold, firstSubjectRow ];
  manualVolume<-manual[1];
  for( i in 2:NoSubject)
  {
    annVolume <- rbind( annVolume, 
                        ann[ threshold, i+1 ] );
    manualVolume<-rbind( manualVolume, manual[i] );

  }

  vol.Maximum <- max( c(annVolume, manualVolume) );
  vol.Minimum <- min( c(annVolume, manualVolume) );

  vol.range <- vol.Maximum-vol.Minimum;

  range.max<-vol.Maximum + 0.1*vol.range;
  range.min<-vol.Minimum - 0.1*vol.range;

  # expand right side of clipping rect to make room for the legend
  par( xpd=T, mar=par()$mar+c(0,0,0,4) );

  plot( manualVolume, annVolume,
        xlim=c(range.min,range.max), ylim=c(range.min,range.max),
        col=seq(1,NoSubject,1),
        pch=c(21,22,23,24,25),
        xlab="Manual Volume", ylab="BRAINSCut Volume");

  # Simple linear fit
  fit<-lm(annVolume~manualVolume);
  fit.intercept <- round( fit$coefficients[[1]], 2 );
  fit.slope <- round( fit$coefficients[[2]], 2 );
  mtext( paste( "y = ",fit.slope , " x + ", fit.intercept ) ,
         side=1, line=-1 );

  lines( c( range.min, range.max ), 
         y= c(( range.min*fit.slope+fit.intercept), range.max*fit.slope+fit.intercept),
         col="red") 
  lines( c( range.min, range.max), y=c( range.min, range.max),
         lty=2);
            
  Legend.x <- range.max + 0.1*vol.range ;
  Legend.y <- range.max + 0.1*vol.range ;
  
  print( ann );
  legend( Legend.x, Legend.y, 
          colnames(ann)[1:NoSubject+1], 
          cex=0.8,
          pch=c(21,22,23,24,25),col=seq(1,NoSubject,1),
          bty="n");
  print( colnames(ann)[1:NoSubject+1] );
  print( seq(1,NoSubject,1) );

  # return slope and intercept
  c( fit.slope, fit.intercept );
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

manualData     <- table(rep(0,numberOfSubject));
manualData[1]  <- currentDT$manual[1];

pdf( plotName, width=14, height=7 );
par(mfrow=c(1,2)); # two plots in one row

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
meanSIPlotResults<-MeanSIPlot( cumulativeSI );

# ICC accros threshold
iccResults<-ICCPlot( manualData, cumulativeVolumes );

# plot volumetric match
#linearFitResults<-VolumetricComparisonPlot( manualData, cumulativeVolumes, iccResults[1]);
linearFitResults<-VolumetricComparisonPlot( manualData, cumulativeVolumes, 50);
dev.off();

# summary stats
summaryStatOutputFilename <- paste(plotName, "txt", sep=".");
outputData<-c( currentROIName, meanSIPlotResults, iccResults[2:5], linearFitResults ); 
outputData<-t(outputData);

write.table( outputData , file=summaryStatOutputFilename , append=FALSE ,
             col.names=c( "roi","SI","SIThreshold", "iccA", "iccAThreshold","iccC","iccCThreshold","slope","intercept") ,
             row.names=FALSE,
             sep=",", quote=FALSE);  


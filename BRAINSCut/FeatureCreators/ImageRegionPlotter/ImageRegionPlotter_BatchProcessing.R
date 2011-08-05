#
# Author : Eun  Young( Regina ) Kim
# 2010 Feb
#
# This program produces joint histogram data from the series of histogram. 
# 
# 1. read in subject data file
# 2. add to the total_histogram
# 3. continue until converges
# 
# 
#  USAGE::
#  $ R --slave --args  \
#                       [listFilename] \
#                       [histogram output filename] \
# 
# 
# 

myArg <- commandArgs();
ListFilename <- as.character( myArg[4] ); # Directory List File
HistogramOutputFilename <- as.character( myArg[5] ); # Directory List File

skiplineForHistogramFile=3;


# 
# read in directory list file to be processed
#
# file ex)
#   /paulsen/IPIG/predict_3T_MR/site-024/0706/43258/histogram.txt
#   /paulsen/IPIG/predict_3T_MR/site-024/0869/62075/histogram.txt
#   /paulsen/IPIG/predict_3T_MR/site-024/1046/83973/histogram.txt
#   ...
# 
# 

listOfFiles <- read.table( ListFilename , stringsAsFactors=F);
colnames( listOfFiles ) <- c("filename");
print( head( listOfFiles ));

#
# iterative procedure by line
#

for ( line in  seq(1:nrow( listOfFiles) ) ){


  filename <- as.character( listOfFiles$filename[line] );
  print( filename )
  # 
  # read in individual file
  # 
  # data file ex) --> Comma separated file! 
  #
  #       [Image1]: /BSITKBRAINSABC/CORRECTED_INDEX_0_LEVEL_4.nii.gz
  #       [Image2]: /BSITKBRAINSABC/CORRECTED_INDEX_2_LEVEL_4.nii.gz
  #       [LabelMap]: /2334_43971_3DSPGR_0_COR_ACPC_DirtyLabels_BRAINSABC.nii.gz
  #       label, INDEX_0_LEVEL_4.nii.gz, CORRECTED_INDEX_2_LEVEL_4.nii.gz, \
#         FrequenceOfQuantile
  #       1,0,0,0
  #       1,0,1,0
  #       1,0,2,0
  #       1,0,3,0
  #       1,0,4,0
  #       ... ...
  #       
  #       

  dt <- read.csv( filename , header=T, skip= skiplineForHistogramFile );
  colnames(dt)<-c("label","bin1", "bin2","freq");

  #
  # if this is the first file, histogramOutput has to be generated
  # from dt file
  #
  #
  # create output file
  # - label : label of tissue
  # - normalizedFreqSum = sum ( freq / totalfreq )
  # - freqSum = sum ( freq )
  #
  if( line == 1 )
  {
    outputHistogram <- subset(dt, select=c("label","bin1","bin2") );
    outputHistogram$normalizedFreqSum <-0;
    outputHistogram$freqSum <-0;

    print ( head( outputHistogram ));
  }


  #
  # initialize percentage colume
  #
  dt$percent <- 0;

  # 
  # For each label in the file
  # ( This has to be consistent across each file in the list )
  # 
  listOfLabels <- levels(factor(dt$label)) ;

  subjectTotalFrequencyByLabel <- matrix( 0:0, ncol= length( listOfLabels ) );
  subjectTotalFrequencyByLabel <- as.data.frame( subjectTotalFrequencyByLabel );

  #
  # give each name as label2, label5, and so on.
  #
  colnames( subjectTotalFrequencyByLabel ) <- paste("label",listOfLabels,sep="");

  # print
  print( head( dt ));

  for( currentLabel in listOfLabels ){

    # print
    print( paste( "Processing label number ", currentLabel ) );

    # 
    # get only data for current label
    # 

    labelDt<- subset( dt,  dt$label == currentLabel  );

    # print
    print( head( labelDt) )


    # 
    # get sum for each label
    # 

    assign( paste( "subjectTotalFrequencyByLabel$label",currentLabel,sep=""),
            sum( labelDt$freq )  );
    print( get(paste( "subjectTotalFrequencyByLabel$label",currentLabel,sep="") ) );

  } # end of label iteration

  #
  # calculate percentage of each label
  #
  dt$percent <- 
    dt$freq / get ( paste( "subjectTotalFrequencyByLabel$label", dt$label,sep="") );

  # 
  # add to the outputHistogram
  # 
  outputHistogram <- merge( outputHistogram, dt, 
                            by=c("label","bin1","bin2") );
  # 
  #  print 
  # 
  print( "-------------------outputHistogram after merge -------------------");
  print( head ( outputHistogram) );
  outputHistogram$normalizedFreqSum <- outputHistogram$normalizedFreqSum + 
                                       outputHistogram$percent;

  outputHistogram$freqSum <- outputHistogram$freqSum + outputHistogram$freq;
                                       

  print( head ( outputHistogram) );
  outputHistogram <- subset( outputHistogram, 
                             select=c("label", "bin1", "bin2",
                                      "normalizedFreqSum","freqSum" ) );

  print( head ( outputHistogram) );

  #
  # write out the table
  #
  write.csv( outputHistogram , HistogramOutputFilename , 
               col.names=T, row.names=F);

}# end of file list iteration

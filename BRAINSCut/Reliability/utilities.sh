#!/bin/bash

function echoOutTo( ){
  echo $1 >> $2
}

function qsubHeader( ){

   if [ $# != 1 ]; then
     echo "Incorrect Number of Argument:: $#";
     echo "Usage:::::" 
     echo "::::::::::"
     echo "$0 outputFilename "
     echo "::::::::::"
     exit 1;
   fi

   outputFilename=$1;
   rm -f $outputFilename;

   macBRANISBuild="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin-20111028/lib"
   heliumBRAINSBuild="/scratch/PREDICT/regina/BRAINS/buildICC/lib"

   echoOutTo "#!/bin/bash"                             $outputFilename 
   echoOutTo "#$ -N BCut${testIteration}"              $outputFilename
   echoOutTo "#$ -j yes"                               $outputFilename
   echoOutTo "#$ -o $outputFilename.log"               $outputFilename
   echoOutTo "#$ -l mf=2G "                            $outputFilename

   echoOutTo "arch=\`uname\`;"                         $outputFilename
   echoOutTo "if [ \"\$arch\" == \"Darwin\" ]; then "  $outputFilename
   echoOutTo "    BRAINSBuild=\"$macBRANISBuild\"   "  $outputFilename
   echoOutTo "else "                                   $outputFilename
   echoOutTo "    BRAINSBuild=\"$heliumBRAINSBuild\""  $outputFilename
   echoOutTo "fi "                                     $outputFilename
}

function printCommandAndRun( ) {
  command=$1;
  echo $command;
  $command;
}

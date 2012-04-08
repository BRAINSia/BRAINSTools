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

   macBRANISBuild="/ipldev/scratch/eunyokim/src/BRAINS20111028/build-Darwin/lib"
   heliumBRAINSBuild="/scratch/PREDICT/regina/BRAINS/buildICC/lib"

   echoOutTo "#!/bin/bash"                             $outputFilename 
   echoOutTo "#$ -N BCut${testIteration}"              $outputFilename
   echoOutTo "#$ -j yes"                               $outputFilename
   echoOutTo "#$ -o $outputFilename.log"               $outputFilename
   echoOutTo "#$ -l mf=2G "                            $outputFilename
   echoOutTo "#$ -pe smp1 1-12"                        $outputFilename
   echoOutTo "PLATFORM=\$(uname)"                      $outputFilename
   echoOutTo "hostname"                                $outputFilename
   echoOutTo "uname -a"                                $outputFilename
   echoOutTo "## Set global number of threads to 4 in order to minimize CPU wasting." $outputFilename
   echoOutTo "export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NSLOTS;"                       $outputFilename
   echoOutTo "echo \"USING NUM THREADS \${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}\" "       $outputFilename

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

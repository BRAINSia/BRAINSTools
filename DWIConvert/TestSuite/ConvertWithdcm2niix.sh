#!/bin/bash -e -x
# \author Hans J. Johnson
test_name=test99

cp $0 ~/Dropbox/Downloads
cp tensor.mrml ~/Dropbox/Downloads
cp /Users/johnsonhj/src/BT-11/BRAINSTools-build/DWIConvert/SlicerDWIProcessing.sh ~/Dropbox/Downloads

DWICONVERTBIN="/Users/johnsonhj/src/BT-11/bin/DWIConvert"
SLICERPROC=/Users/johnsonhj/src/BT-11/BRAINSTools-build/DWIConvert/SlicerDWIProcessing.sh

## INFO: READ THIS: http://neurohut.blogspot.com/2015/11/how-to-extract-bval-bvec-from-dicom.html

make -j 4
if [[ $? -ne 0 ]]; then
  echo "FAILURE TO BUILD"
  exit -1
fi

DICOM_LIST="DTI_004 GeSignaHDx GeSignaHDxBigEndian GeSignaHDxt PhilipsAchieva1 PhilipsAchieva2 PhilipsAchieva3 PhilipsAchieva4 PhilipsAchieva6 PhilipsAchieva7 PhilipsAchievaBigEndian1 SiemensTrio-Syngo2004A-1 SiemensTrio-Syngo2004A-2 SiemensTrioTim1 SiemensTrioTim2 SiemensTrioTim3 SiemensTrioTimBigEndian1 SiemensVerio"
# DICOM_LIST="SiemensTrioTim1 SiemensTrioTim2 SiemensTrioTim3 SiemensTrioTimBigEndian1 SiemensVerio"
#DICOM_LIST="GeSignaHDx"
# DICOM_LIST="PhilipsAchieva1"
# DICOM_LIST="GeSignaHDxBigEndian GeSignaHDxt"
# DICOM_LIST="DTI_004 GeSignaHDxt PhilipsAchieva1"
#DICOM_LIST="SiemensTrioTim1"
#DICOM_LIST="SiemensVerio DTI_004 PhilipsAchieva1"

GUI_VIEW_COMMANDS=${test_name}/${test_name}_gui.list
mkdir -p ${test_name}
echo "# Commands to view with GUI" > ${GUI_VIEW_COMMANDS}
echo "### ========================"  >> ${GUI_VIEW_COMMANDS}

FSLVIEW=fslview
#FSLVIEW=echo

bdir=/Users/johnsonhj/src/BT-11/BRAINSTools-build/DWIConvert
make -C ${bdir}

source  ${FSLDIR}/etc/fslconf/fsl.sh


for dicomDir in $(echo ${DICOM_LIST}); do
  echo "#===  ${dicomDir} ===================="  >> ${GUI_VIEW_COMMANDS}
  outdir=$(pwd)/${test_name}/${dicomDir}
  if [[ -d ${outdir} ]]; then
     echo "SKIPPING: ${outdir}"
     continue;
  fi
  mkdir -p ${outdir}

  ~/src/BT-11/DCMTK-build/bin/dcm2xml $(find $(pwd)/${dicomDir} |head -n 2 |tail -n 1) > ${outdir}/SecondFile.xml

cp ~/Desktop/tt.nii.gz ${outdir}/T1.nii.gz

echo "#================="
echo "#================="
echo "#================="
echo "#================="

if [[ 1 -eq 1 ]] ; then
  PREVPFX=""
  PFX=NIIXFSL
  /opt/dcm2niix/bin/dcm2niix -m -b y -t y -v y -z i  -o ${outdir} -f "${PFX}${dicomDir}" ${dicomDir}
  ### dcm2niix sometimes adds undesirable prefixes to names, and it is not always clear why?
  for ff in $(find ${outdir} -name "_c2${PFX}*"); do
      mv ${ff} ${ff//_c2${PFX}/${PFX}}
  done
  ~/src/BT-11/BRAINSTools-build/BRAINSCommonLib/TestSuite/testbin/DumpImageInfo 4 ${outdir}/${PFX}${dicomDir}.nii.gz > ${outdir}/${PFX}${dicomDir}.nii_itk
   /opt/Clibs/bin/nifti_tool -disp_hdr -in ${outdir}/${PFX}${dicomDir}.nii.gz > ${outdir}/${PFX}${dicomDir}.nii_nifti
  mkdir -p ${outdir}/${PFX}_dti
  /opt/fsl/bin/bet2 ${outdir}/${PFX}${dicomDir}.nii.gz ${outdir}/${PFX}_dti/${PFX}${dicomDir}mask.nii.gz
  /opt/fsl/bin/dtifit \
    --data=${outdir}/${PFX}${dicomDir}.nii.gz \
    --out=${outdir}/${PFX}_dti/${PFX}${dicomDir} \
    --mask=${outdir}/${PFX}_dti/${PFX}${dicomDir}mask.nii.gz \
    --bvecs=${outdir}/${PFX}${dicomDir}.bvec \
    --bvals=${outdir}/${PFX}${dicomDir}.bval \
    --wls \
    --sse \
    --save_tensor
   echo "### dcm2niix " >> ${GUI_VIEW_COMMANDS}
   echo "${FSLVIEW} ${outdir}/${PFX}_dti/${PFX}${dicomDir}mask.nii.gz ${outdir}/${PFX}_dti/${PFX}${dicomDir}_V1.nii* " >> ${GUI_VIEW_COMMANDS}
   echo "" >> ${GUI_VIEW_COMMANDS}

fi

if [[ 1 -eq 1 ]]; then
  PREVPFX=${PFX}
  PFX=FSLToNrrd
  ${DWICONVERTBIN} \
        --conversionMode FSLToNrrd \
        --allowLossyConversion \
        --useIdentityMeaseurementFrame \
        --inputVolume ${outdir}/${PREVPFX}${dicomDir}.nii.gz \
        --outputDirectory ${outdir} \
        --outputVolume ${PFX}${dicomDir}.nhdr |tee ${outdir}/${PFX}.log
  ${SLICERPROC} ${outdir}/${PFX}${dicomDir}.nhdr
  echo "### DWIConvert FSLToNrrd " >> ${GUI_VIEW_COMMANDS}
  sed  "s#DTIIMAGEHOLDER#${PFX}${dicomDir}_SlicerDTI.nrrd#g" tensor.mrml \
  |sed "s#MODELNAME#${PFX}${dicomDir}#g" \
  |sed "s#DWIIMAGEHOLDER#${PFX}${dicomDir}_Baseline.nrrd#g" > ${outdir}/${PFX}${dicomDir}.mrml
  echo "cd  ${outdir};  /Applications/Slicer.app/Contents/MacOS/Slicer  ${outdir}/${PFX}${dicomDir}.mrml" >> ${GUI_VIEW_COMMANDS}
  echo "" >> ${GUI_VIEW_COMMANDS}
fi

if [[ 1 -eq 1 ]] ; then
  PREVPFX=${PFX}
  PFX=DicomToNrrd
  ${DWICONVERTBIN} \
        --conversionMode DicomToNrrd \
        --allowLossyConversion \
        --useIdentityMeaseurementFrame \
        --inputDicomDirectory ${dicomDir} \
        --outputDirectory ${outdir} \
        --outputVolume ${PFX}${dicomDir}.nhdr  |tee ${outdir}/${PFX}.log
  ${SLICERPROC} ${outdir}/${PFX}${dicomDir}.nhdr
  echo "### DWIConvert DicomToNrrd " >> ${GUI_VIEW_COMMANDS}
  sed  "s#DTIIMAGEHOLDER#${PFX}${dicomDir}_SlicerDTI.nrrd#g" tensor.mrml \
  |sed "s#MODELNAME#${PFX}${dicomDir}#g" \
  |sed "s#DWIIMAGEHOLDER#${PFX}${dicomDir}_Baseline.nrrd#g" > ${outdir}/${PFX}${dicomDir}.mrml
  echo "cd  ${outdir};  /Applications/Slicer.app/Contents/MacOS/Slicer  ${outdir}/${PFX}${dicomDir}.mrml" >> ${GUI_VIEW_COMMANDS}
   echo "" >> ${GUI_VIEW_COMMANDS}
fi

if [[ 1 -eq 1 ]] ; then
  PREVPFX=${PFX}
  PFX=DicomToFSL
  ${DWICONVERTBIN} \
        --conversionMode DicomToFSL \
        --allowLossyConversion \
        --useIdentityMeaseurementFrame \
        --inputDicomDirectory ${dicomDir} \
        --outputDirectory ${outdir} \
        --outputVolume ${PFX}${dicomDir}.nii  |tee ${outdir}/${PFX}.log
  ~/src/BT-11/BRAINSTools-build/BRAINSCommonLib/TestSuite/testbin/DumpImageInfo 4 ${outdir}/${PFX}${dicomDir}.nii > ${outdir}/${PFX}${dicomDir}.nii_itk
   /opt/Clibs/bin/nifti_tool -disp_hdr -in ${outdir}/${PFX}${dicomDir}.nii > ${outdir}/${PFX}${dicomDir}.nii_nifti
   /Users/johnsonhj/src/BT-11/bin/BFileCompareTool ${outdir}/${PFX}${dicomDir}.bvec ${outdir}/${PFX}${dicomDir}.bvec 3 1e-2 > ${outdir}/bvec_compare.output

  mkdir -p ${outdir}/${PFX}_dti
  /opt/fsl/bin/bet2 ${outdir}/${PFX}${dicomDir}.nii.gz ${outdir}/${PFX}_dti/${PFX}${dicomDir}mask.nii.gz
  /opt/fsl/bin/dtifit \
    --data=${outdir}/${PFX}${dicomDir}.nii.gz \
    --out=${outdir}/${PFX}_dti/${PFX}${dicomDir} \
    --mask=${outdir}/${PFX}_dti/${PFX}${dicomDir}mask.nii.gz \
    --bvecs=${outdir}/${PFX}${dicomDir}.bvec \
    --bvals=${outdir}/${PFX}${dicomDir}.bval \
    --wls \
    --sse \
    --save_tensor
   echo "### DWIConvert DicomToFSL " >> ${GUI_VIEW_COMMANDS}
   echo "${FSLVIEW} ${outdir}/${PFX}_dti/${PFX}${dicomDir}mask.nii.gz ${outdir}/${PFX}_dti/${PFX}${dicomDir}_V1.nii* " >> ${GUI_VIEW_COMMANDS}
   echo "" >> ${GUI_VIEW_COMMANDS}
fi

if [[ 1 -eq 1 ]] ; then
  PREVPFX=${PFX}
  PFX=NIIXFSLToNrrd
  ${DWICONVERTBIN} \
        --conversionMode FSLToNrrd \
        --allowLossyConversion \
        --useIdentityMeaseurementFrame \
        --inputVolume ${outdir}/${PREVPFX}${dicomDir}.nii \
        --outputDirectory ${outdir} \
        --outputVolume ${PFX}${dicomDir}.nhdr  |tee ${outdir}/${PFX}.log
  ~/src/BT-11/BRAINSTools-build/BRAINSCommonLib/TestSuite/testbin/DumpImageInfo 4 ${outdir}/${PFX}${dicomDir}.nhdr > ${outdir}/${PFX}${dicomDir}.nhdr_itk
  ${SLICERPROC} ${outdir}/${PFX}${dicomDir}.nhdr
  echo "### DWIConvert FSLToNrrd " >> ${GUI_VIEW_COMMANDS}
  sed  "s#DTIIMAGEHOLDER#${PFX}${dicomDir}_SlicerDTI.nrrd#g" tensor.mrml \
  |sed "s#MODELNAME#${PFX}${dicomDir}#g" \
  |sed "s#DWIIMAGEHOLDER#${PFX}${dicomDir}_Baseline.nrrd#g" > ${outdir}/${PFX}${dicomDir}.mrml
  echo "cd  ${outdir};  /Applications/Slicer.app/Contents/MacOS/Slicer  ${outdir}/${PFX}${dicomDir}.mrml" >> ${GUI_VIEW_COMMANDS}
  echo "" >> ${GUI_VIEW_COMMANDS}
fi

echo "" >> ${GUI_VIEW_COMMANDS}
echo "### -!-!-!-!" >> ${GUI_VIEW_COMMANDS}
echo "" >> ${GUI_VIEW_COMMANDS}

done

for i in $(find ${test_name}/ -name bvec_compare.output ); do
   echo "$i $(tail -n 1 $i)"
done


NOTES
================

These are notes to assist with knowing how to manually intervien when BCD
encounters one of it's rare failure cses.

step1_forceFailure.sh
step2_runEMSP.sh
step3_finalRun.sh


HINTS OF WHERE TO LOOK
==============================
JIRA
----
https://www.icts.uiowa.edu/jira/browse/PREDICTIMG-3181
https://www.icts.uiowa.edu/jira/browse/PREDICTIMG-3811


1)
Create an input list from the "denoised" images:

find /Shared/johnsonhj/HDNI/Neuromorphometrics/2012Subscription/Data -name "*_UNM_denoised.nii.gz" >> inputFiles.txt

2)
I ran the "1_forceReportFailure.sh" file, so BCD creates EMSP.nrrd/ EMSP.fcsv files. The bash script is attached here.

3)
For those cases that BCD could estimate a good MSP image:
A- I Loaded EMSP.nrrd, EMSP.fcsv files by Slicer and define eyes.
B- Changed EMSP.fcsv to have proper format.
C- Renamed the EMSP.nrrd to "${SUBJECT}_maneyes_denoised.nrrd".
D- Created the following file adjacent to the nrrd file to indicates that denoising should be skipped:
touch ${SUBJECT}_maneyes_denoised.nrrd_noDenoise

4)
For many cases the "EMSP.nrrd" image is just a crap because of bad MSP estimation.
For such cases:
A- I removed the created EMSP.fcsv file
B- I re-ran the 1_forceReportFailure.sh script on the previously created "EMSP.nrrd" file!
C- Now, new EMSP.nrrd/ EMSP.fcsv files are created with good qualities.
D- I ran all steps of the stage 3 on the above EMSP.nrrd/ EMSP.fcsv files.

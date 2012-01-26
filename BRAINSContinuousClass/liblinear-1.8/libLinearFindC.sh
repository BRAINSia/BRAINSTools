#time ./train -c 1/6 ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC1sixth.model > /dev/null
#time ./predict ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC1sixth.model /dev/null > ../../Training/libLinearResultsC1sixth.txt
#time ./train -c 1/4 ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC0.25.model > /dev/null
#time ./predict ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC0.25.model /dev/null > ../../Training/libLinearResultsC0.25.txt
#time ./train -c 1/2 ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC0.50.model > /dev/null
#time ./predict ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC0.50.model /dev/null > ../../Training/libLinearResultsC0.5.txt
#time ./train -c 1 ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC1.model > /dev/null
time ./predict ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC1.model /dev/null > ../../Training/libLinearResultsC1.txt
#time ./train -c 2 ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC2.model > /dev/null
time ./predict ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC2.model /dev/null > ../../Training/libLinearResultsC2.txt
time ./train -c 4 ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC4.model > /dev/null
time ./predict ../../Training/allDataScaled.dat ../../Training/libLinearAllDataScaledC4.model /dev/null > ../../Training/libLinearResultsC4.txt

//
// Created by Hui Xie on 2/28/17.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>       /* fabs */
#include <cstdlib>      /* atof */
#include <string.h>     /* strcmp */

void printUsage(char* argv0){
  std::cout<<"Function: check whether a specific tag-value pair in header file match with a input value."<< std::endl;
  std::cout<<"Usage:"<<std::endl
           << argv0 << " --f HeaderFileName --tag TagName --section sectionIndex --subsection subsectionIndex --v verifyingValue --numtype <1|0>" << std::endl
           << "Notes: index is from 0 to count after the tag" <<std::endl;
  std::cout<<"Example:" <<std::endl
          <<"For a line text: space directions: (0.9343162500000001,0.07719215625,0) (-0.07719215625,0.9343162500000001,0) (0,0,5.999994770777341) none ,"<<std::endl
          <<"--tag \"space direction\" --section 2 --subsection 2  --v 6  --numtype 1 will check whether 5.999994770777341 is consistent with 6." << std::endl;
  return;
}

int main(int argc, char * argv[])
{
  int result = 1; // 1 means unmatching;
  //Get input parameters
  std::string inputFilename;
  std::string tagStr;
  int sectionIndex = 0;
  int subsectionIndex = 0;
  std::string verifyingValue;
  bool numType = true;

  if(13 != argc){
    std::cout<<"Parameter Error: it should have 12 parameters."<<std::endl;
    printUsage(argv[0]);
    return -1;
  }
  for (int i=1; i<12; i=i+2){
    if(0 == strcmp("--f", argv[i])){
       inputFilename = argv[i+1];
    }
    else if(0 == strcmp("--tag", argv[i])){
       tagStr = argv[i+1];
    }
    else if(0 == strcmp("--section", argv[i])){
      sectionIndex = atoi(argv[i+1]);
    }
    else if(0 == strcmp("--subsection", argv[i])){
      subsectionIndex = atoi(argv[i+1]);
    }
    else if(0 == strcmp("--v", argv[i])){
      verifyingValue = argv[i+1];
    }
    else if(0 == strcmp("--numtype", argv[i])){
      numType = atoi(argv[i+1]) > 0;
    }
    else{
      std::cout<<"Parameter Error: unmatched option."<<std::endl;
      std::cout<<"All options are small case characters."<<std::endl;
      printUsage(argv[0]);
      return -1;
    }
  }

  //check whether tag has blank space inside
  int nSpace = 0;
  std::size_t  pos = 0;
  pos = tagStr.find(' ',pos);
  while(std::string::npos != pos){
    ++nSpace;
    pos = tagStr.find(' ',pos+1);
  }

  //read file
  std::ifstream file(inputFilename.c_str());
  std::string lineText;
  while (std::getline(file, lineText))
  {
    std::istringstream lineStream(lineText);
    std::string fileTag;
    fileTag.clear();
    for(int i=0; i<=nSpace;++i){
      std::string tag;
      std::getline(lineStream,tag,' ');
      if (0 != tag.length()){
         fileTag += tag +" ";
      }
    }
    fileTag.erase(fileTag.length()-1,1); // trim tail whitespace
    if (fileTag.length() > 2){
      fileTag.erase(fileTag.length()-1,1); // trim tail ':'
    }

    if (0 == tagStr.compare(fileTag)){
      std::string sectionStr;
      for(int i=0; i<=sectionIndex;++i){
        sectionStr.clear();
        std::getline(lineStream, sectionStr,' ');
        if (0 == sectionStr.length()) --i;
      }
      std::cout<<"finding section = " << sectionStr<< std::endl;

      //erease bracket at both ends
      if ('('==sectionStr.at(0)){
         sectionStr.erase(0, 1);
      }
      if (')'==sectionStr.at(sectionStr.length()-1)){
        sectionStr.erase(sectionStr.length()-1, 1);
      }


      std::istringstream sectionStream(sectionStr);
      std::string valueStr;
      for(int i=0; i<=subsectionIndex;++i){
        valueStr.clear();
        std::getline(sectionStream,valueStr,',');
      }

      //compare the value
      if (numType){
        double verifyingNum = atof(verifyingValue.c_str());
        double tagNum = atof(valueStr.c_str());
        if (fabs(verifyingNum-tagNum) <= 1e-3)
        {
          result = 0;
          break;
        }
        else {
          std::cout<<"input verifying value: "<< verifyingValue <<std::endl;
          std::cout<<"value read from file: "<< valueStr << std::endl;
          result = 1;
          break;
        }
      }
      else{
        if (0 == verifyingValue.compare(valueStr))
        {
          result = 0;
          break;
        }
        else {
          std::cout<<"input verifying value: "<< verifyingValue <<std::endl;
          std::cout<<"value read from file: "<< valueStr << std::endl;
          result = 1;
          break;
        }
      }
    }
  }
  file.close();
  if (0 != result){
    std::cout<<"Can not find matching key-value with input"<< std::endl;
  }
  else{
    std::cout<<"The input value is consistent with key-value in the header file."<< std::endl;
  }

  return result;
}

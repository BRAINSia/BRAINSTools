// Author: Ali Ghayoor

// INCLUDES ////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>

#include <glob.h>
#include <wordexp.h>

#include "LandmarksCompareCLP.h"

// DEFINES /////////////////////////////////////////////////////////////////////

// Check to see if a field is a comment
static bool is_comment(const std::string& field)
{
  for( unsigned int i = 0; i < field.length(); i++ )
    {
    if( (field[i] != ' ') && (field[i] != '\t') )
      {
      if( i + 1 >= field.length() )
        {
        return false;
        }
      if( field[i] == '#' )
        {
        return true;
        }
      else
        {
        return false;
        }
      }
    }
  return false;
}

// Read in a csv file
static std::vector<std::vector<std::string> > read_csv(const std::string& file)
{
  std::vector<std::vector<std::string> > ret;
  std::ifstream                          infile(file.c_str() );
  std::string                            line;

  while( getline(infile, line) )
    {
    if( is_comment(line) )
      {
      continue;
      }
    std::istringstream       linestream(line);
    std::string              item;
    std::vector<std::string> line_vec;
    while( getline(linestream, item, ',') )
      {
      line_vec.push_back(item);
      }

    ret.push_back(line_vec);
    }

  return ret;
}

// MAIN ////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const std::vector<std::vector<std::string> > file1 = read_csv(inputLandmarkFile1);
  const std::vector<std::vector<std::string> > file2 = read_csv(inputLandmarkFile2);

  if( file1 != file2 )
    {
    std::cout << "Two landmark files are not identical!" << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "The landmark files are identical!" << std::endl;
  return EXIT_SUCCESS;
}

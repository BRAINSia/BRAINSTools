//
// Created by Hans Johnson on 10/8/16.
//
#include "DWIConvertUtils.h"


void PrintVec(const vnl_vector_fixed<double,3> & vec)
{
  std::cerr << "[";

  for( unsigned i = 0; i < vec.size(); ++i )
  {
    PrintVec<double>(vec[i]);
    if( i < vec.size() - 1 )
    {
      std::cerr << " ";
    }
  }
  std::cerr << "]" << std::endl;

}

void PrintVec(const DWIMetaDataDictionaryValidator::GradientTableType & vec)
{
  std::cerr << "[";

  for( unsigned i = 0; i < vec.size(); ++i )
  {
    PrintVec(vec[i]);
    if( i < vec.size() - 1 )
    {
      std::cerr << " ";
    }
  }
  std::cerr << "]" << std::endl;
}

bool CloseEnough(const vnl_vector_fixed<double,3> & a, const vnl_vector_fixed<double,3> & b, double magdiv)
{
  if( a.size() != b.size() )
  {
    std::cerr << "Vector size mismatch: "
              << a.size() << " "
              << b.size() << std::endl;
    return false;
  }
  for( unsigned i = 0; i < a.size(); ++i )
  {
    if( !CloseEnough(a[i], b[i], magdiv) )
    {
      std::cerr << "Value mismatch" << std::endl;
      return false;
    }
  }
  return true;
}

void normalize(const DWIMetaDataDictionaryValidator::GradientDirectionType &vec,double *normedVec)
{
  double norm = 0.0;
  for(unsigned j = 0; j < 3; ++j)
  {
    normedVec[j] = vec[j];
    norm += vec[j] * vec[j];
  }
  norm = std::sqrt(norm);
  if( norm < 0.00001) //Only norm if not equal to zero
  {
    for(unsigned j = 0; j < 3; ++j)
    {
      normedVec[j] = 0.0;
    }
  }
  else if( std::abs(1.0 - norm) > 1e-4 ) //Only normalize if not very close to 1
  {
    for(unsigned j = 0; j < 3; ++j)
    {
      normedVec[j] /= norm;
    }
  }
}


int
WriteBVectors(const DWIMetaDataDictionaryValidator::GradientTableType & bVectors,
              const std::string & filename)
{
  itk::NumberToString<double> DoubleConvert;
  std::ofstream  bVecFile;

  bVecFile.open(filename.c_str(), std::ios::out | std::ios::binary);
  bVecFile.precision(17);
  if( !bVecFile.is_open() || !bVecFile.good() )
  {
    return EXIT_FAILURE;
  }
  //Matching formatting from dcm2niix that seems to be gaining acceptance as
  //the format for these files
  //
  // The lines [1,2,3] is the [x,y,z]component of the gradient vector for each gradient image
  for( unsigned int index=0; index < 3; ++index)
  {
    for( unsigned int k = 0; k < bVectors.size(); ++k )
    {
      const char * const spacer = (  k==bVectors.size()-1 ) ? "" : " ";
      double normedVec[3];
      normalize(bVectors[k],normedVec);
      bVecFile << DoubleConvert(normedVec[index]) << spacer;
    }
    bVecFile << std::endl;
  }

  bVecFile.close();
  return EXIT_SUCCESS;
}

int
ReadBVals(std::vector<double> & bVals, unsigned int & bValCount, const std::string & bValFilename, double & maxBValue)
{
  std::ifstream bValFile(bValFilename.c_str(), std::ifstream::in);

  if( !bValFile.good() )
  {
    std::cerr << "Failed to open " << bValFilename
              << std::endl;
    return EXIT_FAILURE;
  }
  bVals.clear();
  bValCount = 0;
  while( !bValFile.eof() )
  {
    double x;
    bValFile >> x;
    if( bValFile.fail() )
    {
      break;
    }
    if( x > maxBValue )
    {
      maxBValue = x;
    }
    bValCount++;
    bVals.push_back(x);
  }

  return EXIT_SUCCESS;
}

int
ReadBVecs(DWIMetaDataDictionaryValidator::GradientTableType & bVecs, unsigned int & bVecCount, const std::string & bVecFilename , bool horizontalBy3Rows )
{
  std::ifstream bVecFile(bVecFilename.c_str(), std::ifstream::in);

  if( !bVecFile.good() )
  {
    std::cerr << "Failed to open " << bVecFilename
              << std::endl;
    return EXIT_FAILURE;
  }
  bVecs.clear();
  bVecCount = 0;
  if( horizontalBy3Rows )
  {
    std::vector<std::vector<double> > bVecst( 3 ) ;
    for( unsigned i = 0 ; i < 3 ; i++ )
    {
      std::string bvect;
      std::getline(bVecFile, bvect);
      bool error = false ;
      if( bVecFile.fail() )
      {
        return EXIT_FAILURE ;
      }
      std::istringstream issLineToString(bvect);
      for( std::istream_iterator<std::string> it = std::istream_iterator<std::string>(issLineToString); it != std::istream_iterator<std::string>(); it++ )
      {
        std::istringstream iss( *it ) ;
        double val ;
        iss >> val ;
        if( iss.fail() )
        {
          error = true ;
          break;
        }
        bVecst[ i ].push_back( val ) ;
      }
      if( error )
      {
        break ;
      }
    }
    for( unsigned int i = 1 ; i < 3 ; i++ )
    {
      if( bVecst[ i ].size() !=  bVecst[ 0 ].size() )
      {
        return EXIT_FAILURE ;
      }
    }
    bVecCount = bVecst[ 0 ].size() ;
    for( unsigned int i = 0 ; i < bVecCount ; i++ )
    {
      double list[] = {bVecst[0][i],bVecst[1][i],bVecst[2][i]} ;
      DWIMetaDataDictionaryValidator::GradientDirectionType x;
      x[0]=list[0];
      x[1]=list[1];
      x[2]=list[2];
      bVecs.push_back(x);
    }
  }
  else
  {
    while( !bVecFile.eof() )
    {
      DWIMetaDataDictionaryValidator::GradientDirectionType x;
      for( unsigned i = 0; i < 3; ++i )
      {
        double val;
        bVecFile >> val;
        if( bVecFile.fail() )
        {
          break;
        }
        x[i]=(val);
      }
      if( bVecFile.fail() )
      {
        break;
      }
      bVecCount++;
      bVecs.push_back(x);
    }
  }
  return EXIT_SUCCESS;
}

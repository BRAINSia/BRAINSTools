
#ifndef __itkTimeSeriesDatabase_h
#define __itkTimeSeriesDatabase_h

#include <itkImage.h>
#include <itkArray.h>
#include <itkImageSource.h>
#include <iostream>
#include <fstream>
#include <itkTimeSeriesDatabaseHelper.h>

#define TimeSeriesBlockSize 16
#define TimeSeriesBlockSizeP2 TimeSeriesBlockSize * TimeSeriesBlockSize
#define TimeSeriesBlockSizeP3 TimeSeriesBlockSize * TimeSeriesBlockSize * TimeSeriesBlockSize
#define TimeSeriesVolumeBlockSize TimeSeriesBlockSize * TimeSeriesBlockSize * TimeSeriesBlockSize
#define TimeSeriesVolumeBlockSizeP3 TimeSeriesVolumeBlockSize * TimeSeriesVolumeBlockSize * TimeSeriesVolumeBlockSize

namespace itk
{

/*
 * TimeSeriesDatabase transforms a series of images stored on disk into a high performance database
 *
 * The main idea behind TimeSeriesDatabase is to have a representation of a 4 dimensional dataset that
 * is larger than main memory, but may still be accessed in a rapid manner.  Though not strictly
 * ITK conforming, this initial pass is strictly 4 dimensional datasets.
 */
template < typename TPixel >
class TimeSeriesDatabase : public ImageSource< Image< TPixel, 3 > >
{
public:
  using Self = TimeSeriesDatabase;
  using Superclass = ImageSource< Image< TPixel, 3 > >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;
  using ConstWeakPointer = WeakPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TimeSeriesDatabase, ImageSource );

  using OutputImageType = Image< TPixel, 3 >;
  using OutputImageTypePointer = typename OutputImageType::Pointer;
  using OutputSliceType = Image< TPixel, 2 >;
  using OutputSliceTypePointer = typename OutputSliceType::Pointer;
  using ArrayType = Array< TPixel >;

  /** Connect to an existing TimeSeriesDatabase file on disk
   * The idea behind the Connect method is to associate this
   * class with a pre-existing self-describing file containing
   * a 4-dimensional dataset that is indexed for rapid retrieval.
   */
  void
  Connect( const char * filename );

  /** Disconnect from a TimeSeriesDatabase file
   * Essentially closes the file and returns the object to its
   * original state.
   */
  void
  Disconnect();

  /** Create a new TimeSeriesDatabase from an Archetype filename
   * Find all the volumes matching the archetype pattern, loading
   * and checking that they are all the same size.  Write the data
   * into a series of files.  The default filesize is 1 GiB, but may
   * be changed using the overloaded method.
   * A call to Connect in required to open the newly created TimeSeriesDatabase.
   */
  static void
  CreateFromFileArchetype( const char * filename, const char * archetype );
  static void
  CreateFromFileArchetype( const char * filename, const char * archetype, unsigned long BlocksPerFile );

  /** Set the image to be read when GenerateData is called.
   * This method selects the image to be returned by an Update
   * call.  By changing the CurrentImage, a pipeline can process
   * each image in the series one after another.
   */
  itkSetMacro( CurrentImage, unsigned int );
  itkGetMacro( CurrentImage, unsigned int );

  /** Return information about the TimeSeriesDatabase file */
  int
  GetNumberOfVolumes()
  {
    return this->m_Dimensions[3];
  };
  itkGetMacro( OutputSpacing, typename OutputImageType::SpacingType );
  itkGetMacro( OutputRegion, typename OutputImageType::RegionType );
  itkGetMacro( OutputOrigin, typename OutputImageType::PointType );
  itkGetMacro( OutputDirection, typename OutputImageType::DirectionType );

  /** Standard method for a ImageSource object */
  virtual void
  GenerateOutputInformation( void ) override;
  virtual void
  GenerateData( void ) override;

  /** A convience method for reading a voxel's time course
   * Subsequent calls to voxels in the immediate region of this will be
   * cached for quick access
   */
  void
  GetVoxelTimeSeries( typename OutputImageType::IndexType idx, ArrayType & array );

  /** Set the size of the cache in MiB (1 MiB = 2^20 bytes)
   */
  void
  SetCacheSizeInMiB( float sz );
  /** Get the size of the cache in MiB (1 MiB = 2^20 bytes)
   */
  float
  GetCacheSizeInMiB();


protected:
  TimeSeriesDatabase();
  ~TimeSeriesDatabase();
  virtual void
                        PrintSelf( std::ostream & os, Indent indent ) const override;
  Array< unsigned int > m_Dimensions;
  Array< unsigned int > m_BlocksPerImage;

  typename OutputImageType::SpacingType   m_OutputSpacing;
  typename OutputImageType::RegionType    m_OutputRegion;
  typename OutputImageType::PointType     m_OutputOrigin;
  typename OutputImageType::DirectionType m_OutputDirection;
  using StreamPtr = itk::TimeSeriesDatabaseHelper::counted_ptr< std::fstream >;

  static std::streampos
  CalculatePosition( unsigned long index, unsigned long BlocksPerFile );

  unsigned int
  CalculateFileIndex( unsigned long Index );
  static unsigned int
  CalculateFileIndex( unsigned long Index, unsigned long BlocksPerFile );

  unsigned long        CalculateIndex( Size< 3 > Position, int ImageCount );
  static unsigned long CalculateIndex( Size< 3 > Position, int ImageCount, unsigned int BlocksPerImage[3] );
  /// Return true if this is a full block, false otherwise.  Assumes there is overlap!
  bool CalculateIntersection( Size< 3 > BlockIndex, typename OutputImageType::RegionType RequestedRegion,
                              typename OutputImageType::RegionType & BlockRegion,
                              typename OutputImageType::RegionType & ImageRegion );
  bool
  IsOpen() const;

  /// How many pixels are in the last block?
  Array< unsigned int > m_PixelRemainder;
  std::string           m_Filename;
  unsigned int          m_CurrentImage;

  std::vector< StreamPtr >   m_DatabaseFiles;
  std::vector< std::string > m_DatabaseFileNames;
  unsigned long              m_BlocksPerFile;

  /// our cache
  struct CacheBlock
  {
    TPixel data[TimeSeriesBlockSize * TimeSeriesBlockSize * TimeSeriesBlockSize];
  };
  TimeSeriesDatabaseHelper::LRUCache< unsigned long, CacheBlock > m_Cache;
  CacheBlock *
  GetCacheBlock( unsigned long index );
};

} // end namespace itk
#include "itkTimeSeriesDatabase.txx"

#endif

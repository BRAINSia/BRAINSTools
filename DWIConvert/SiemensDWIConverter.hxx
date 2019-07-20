//
// Created by Hui Xie on 12/19/16.
//

#ifndef BRAINSTOOLS_SIEMENSDWICONVERTER_HXX
#define BRAINSTOOLS_SIEMENSDWICONVERTER_HXX

template < typename T >
T
SiemensDWIConverter::CSAExtractFromString( const char * ptr )
{
  T rval = *( reinterpret_cast< const T * >( ptr ) );
  itk::ByteSwapper< T >::SwapFromSystemToLittleEndian( &rval );
  return rval;
}


#endif // BRAINSTOOLS_SIEMENSDWICONVERTER_HXX

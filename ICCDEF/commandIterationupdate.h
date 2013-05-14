/*
 *  commandIterationupdate.h
 *  ThirionCLP
 *
 *  Created by Yong Qiang Zhao on 8/21/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include "itkCommand.h"
#include "itkPDEDeformableRegistrationFilter.h"

#include "itkMultiResolutionPDEDeformableRegistration.h"
#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkWarpImageFilter.h"
#include "itkICCDeformableRegistrationFilter.h"

#include "itkDisplacementFieldJacobianDeterminantFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkWarpHarmonicEnergyCalculator.h"
#include "itkGridForwardWarpImageFilter.h"
#include "itkVectorCentralDifferenceImageFunction.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>

#if defined( USE_DEBUG_IMAGE_VIEWER )
#include <DebugImageViewerClient.h>
extern DebugImageViewerClient DebugImageDisplaySender;
#include <itkWarpImageFilter.h>
#endif

/*
class CommandIterationUpdate : public itk::Command
{
public:
typedef  CommandIterationUpdate   Self;
typedef  itk::Command             Superclass;
typedef  itk::SmartPointer<CommandIterationUpdate>  Pointer;
itkNewMacro( CommandIterationUpdate );
protected:
CommandIterationUpdate() {};

typedef itk::Image< float, 3 > InternalImageType;
typedef itk::Vector< float, 3 >    VectorPixelType;
typedef itk::Image<  VectorPixelType, 3 > DisplacementFieldType;

typedef itk::PDEDeformableRegistrationFilter<
InternalImageType,
InternalImageType,
DisplacementFieldType>   RegistrationFilterType;

public:

void Execute(itk::Object *caller, const itk::EventObject & event)
{
Execute( (const itk::Object *)caller, event);
}

void Execute(const itk::Object * object, const itk::EventObject & event)
{
const RegistrationFilterType * filter =
dynamic_cast< const RegistrationFilterType * >( object );

if( !(itk::IterationEvent().CheckEvent( &event )) )
{
return;
}
std::cout <<   filter->GetMetric() << std::endl;
}
};
*/

template <class TPixel = float, unsigned int VImageDimension = 3>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate  Self;
  typedef  itk::Command            Superclass;
  typedef  itk::SmartPointer<Self> Pointer;

  typedef itk::Image<TPixel, VImageDimension>          InternalImageType;
  typedef itk::Vector<TPixel, VImageDimension>         VectorPixelType;
  typedef itk::Image<VectorPixelType, VImageDimension> DisplacementFieldType;

  typedef itk::DemonsRegistrationFilter<
      InternalImageType,
      InternalImageType,
      DisplacementFieldType>   DemonsRegistrationFilterType;

  typedef itk::DiffeomorphicDemonsRegistrationFilter<
      InternalImageType,
      InternalImageType,
      DisplacementFieldType>   DiffeomorphicDemonsRegistrationFilterType;

  typedef itk::ICCDeformableRegistrationFilter<
      InternalImageType,
      InternalImageType,
      DisplacementFieldType>   ICCDeformableRegistrationFilterType;

  typedef itk::FastSymmetricForcesDemonsRegistrationFilter<
      InternalImageType,
      InternalImageType,
      DisplacementFieldType>   FastSymmetricForcesDemonsRegistrationFilterType;

  typedef itk::MultiResolutionPDEDeformableRegistration<
      InternalImageType, InternalImageType,
      DisplacementFieldType, TPixel>   MultiResRegistrationFilterType;

  typedef itk::DisplacementFieldJacobianDeterminantFilter<
      DisplacementFieldType, TPixel, InternalImageType> JacobianFilterType;

  typedef itk::MinimumMaximumImageCalculator<InternalImageType>
    MinMaxFilterType;

  typedef itk::WarpHarmonicEnergyCalculator<DisplacementFieldType>
    HarmonicEnergyCalculatorType;

  typedef itk::VectorCentralDifferenceImageFunction<DisplacementFieldType>
    WarpGradientCalculatorType;

  typedef typename WarpGradientCalculatorType::OutputType WarpGradientType;

  itkNewMacro( Self );
private:
  std::ofstream m_Fid;
  std::ofstream m_Bid;
  bool          m_headerwritten;
  typename JacobianFilterType::Pointer m_JacobianFilter;
  typename MinMaxFilterType::Pointer m_Minmaxfilter;
  typename HarmonicEnergyCalculatorType::Pointer m_HarmonicEnergyCalculator;
  typename DisplacementFieldType::ConstPointer m_TrueField;
  typename WarpGradientCalculatorType::Pointer m_TrueWarpGradientCalculator;
  typename WarpGradientCalculatorType::Pointer m_CompWarpGradientCalculator;
  typename InternalImageType::Pointer m_MovingImage;
  typename InternalImageType::Pointer m_FixedImage;
public:
  void SetMovingImage(typename InternalImageType::Pointer & img)
  {
    m_MovingImage = img;
  }

  void SetFixedImage(typename InternalImageType::Pointer & img)
  {
    m_FixedImage = img;
  }

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event );
  }

  void Execute(const itk::Object *object, const itk::EventObject & event)
  {
    if( !( itk::IterationEvent().CheckEvent( &event ) ) )
      {
      return;
      }

    typename DisplacementFieldType::Pointer deffield = 0;
    typename DisplacementFieldType::Pointer backdeffield = 0;
    unsigned int iter = vcl_numeric_limits<unsigned int>::max();
    double       metricbefore = -1.0;
    double       backmetricbefore = -1.0;

    if( const ICCDeformableRegistrationFilterType * filter
          = dynamic_cast<const ICCDeformableRegistrationFilterType *>(
              object ) )
      {
      iter = filter->GetElapsedIterations() - 1;
      metricbefore = filter->GetForwardMetric();
      backmetricbefore = filter->GetBackwardMetric();
      deffield = const_cast<ICCDeformableRegistrationFilterType *>
        ( filter )->GetOutput(0);
      backdeffield = const_cast<ICCDeformableRegistrationFilterType *>
        ( filter )->GetOutput(1);
      }
    else if( const MultiResRegistrationFilterType * multiresfilter
               = dynamic_cast<const MultiResRegistrationFilterType *>( object ) )
      {
      std::cout << "Finished Multi-resolution iteration :"
                << multiresfilter->GetCurrentLevel() - 1 << std::endl;
      std::cout << "==============================" << std::endl << std::endl;
      }
    else
      {
      return;
      }

    if( deffield.IsNotNull() )
      {
      std::cout << iter << ": MSE " << metricbefore << " - ";

      double fieldDist = -1.0;
      double fieldGradDist = -1.0;
      double tmp;
      if( m_TrueField )
        {
        typedef itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType>
          FieldIteratorType;
        FieldIteratorType currIter(
          deffield, deffield->GetLargestPossibleRegion() );
        FieldIteratorType trueIter(
          m_TrueField, deffield->GetLargestPossibleRegion() );

        m_CompWarpGradientCalculator->SetInputImage( deffield );

        fieldDist = 0.0;
        fieldGradDist = 0.0;
        for( currIter.GoToBegin(), trueIter.GoToBegin();
             not currIter.IsAtEnd(); ++currIter, ++trueIter )
          {
          fieldDist += ( currIter.Value() - trueIter.Value() ).GetSquaredNorm();

          // No need to add Id matrix here as we do a substraction
          tmp = (
              ( m_CompWarpGradientCalculator->EvaluateAtIndex( currIter.GetIndex() )
                - m_TrueWarpGradientCalculator->EvaluateAtIndex( trueIter.
                                                                 GetIndex() )
              ).GetVnlMatrix() ).frobenius_norm();
          fieldGradDist += tmp * tmp;
          }
        fieldDist = sqrt( fieldDist / (double)(
                            deffield->GetLargestPossibleRegion().
                            GetNumberOfPixels() ) );
        fieldGradDist = sqrt( fieldGradDist / (double)(
                                deffield->GetLargestPossibleRegion().
                                GetNumberOfPixels() ) );

        std::cout << "d(.,true) " << fieldDist << " - ";
        std::cout << "d(.,Jac(true)) " << fieldGradDist << " - ";
        }

#if defined( USE_DEBUG_IMAGE_VIEWER )
      if( DebugImageDisplaySender.Enabled() )
        {
        DebugImageDisplaySender.SendImage<DisplacementFieldType>( deffield, 0, 0);
        DebugImageDisplaySender.SendImage<DisplacementFieldType>( deffield, 1, 1);
        DebugImageDisplaySender.SendImage<DisplacementFieldType>( deffield, 2, 2);
        typedef typename itk::WarpImageFilter<InternalImageType,
                                              InternalImageType, DisplacementFieldType> WarpFilterType;
        typename WarpFilterType::Pointer warper = WarpFilterType::New();
        warper->SetInput(m_MovingImage);
        warper->SetOutputSpacing( deffield->GetSpacing() );
        warper->SetOutputOrigin( deffield->GetOrigin() );
        warper->SetOutputDirection( deffield->GetDirection() );
        warper->SetDisplacementField( deffield);
        warper->Update();
        typename InternalImageType::Pointer
        DeformedMovingImagePtr = warper->GetOutput();
        DebugImageDisplaySender.SendImage<InternalImageType>(DeformedMovingImagePtr, 3);
        //       std::cerr << std::endl << "************IMAGES SENT*************" << std::endl;
        }
#endif // defined(USE_DEBUG_IMAGE_VIEWER)

      m_HarmonicEnergyCalculator->SetImage( deffield );
      m_HarmonicEnergyCalculator->Compute();
      const double harmonicEnergy
        = m_HarmonicEnergyCalculator->GetHarmonicEnergy();
      std::cout << "harmo. " << harmonicEnergy << " - ";

      m_JacobianFilter->SetInput( deffield );
      m_JacobianFilter->UpdateLargestPossibleRegion();

      const unsigned int numPix = m_JacobianFilter->
        GetOutput()->GetLargestPossibleRegion().
        GetNumberOfPixels();

      TPixel *pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
      TPixel *pix_end = pix_start + numPix;

      TPixel *jac_ptr;

      // Get percentage of det(Jac) below 0
      unsigned int jacBelowZero(0u);
      for( jac_ptr = pix_start; jac_ptr != pix_end; ++jac_ptr )
        {
        if( *jac_ptr <= 0.0 )
          {
          ++jacBelowZero;
          }
        }
      const double jacBelowZeroPrc = static_cast<double>( jacBelowZero )
        / static_cast<double>( numPix );

      // Get min an max jac
      const double minJac = *( std::min_element(pix_start, pix_end) );
      const double maxJac = *( std::max_element(pix_start, pix_end) );

      // Get some quantiles
      // We don't need the jacobian image
      // we can modify/sort it in place
      std::cout << "max|Jac| " << maxJac << " - "
                << "min|Jac| " << minJac << " - "
                << "ratio(|Jac|<=0) " << jacBelowZeroPrc << std::endl;

      if( this->m_Fid.is_open() )
        {
        if( not m_headerwritten )
          {
          this->m_Fid << "Iteration"
                      << ", MSE before"
                      << ", Harmonic energy"
                      << ", min|Jac|"
                      << ", max|Jac|"
                      << ", ratio(|Jac|<=0)";
          if( backdeffield.IsNotNull() )
            {
            this->m_Bid << "Iteration"
                        << ", MSE before"
                        << ", Harmonic energy"
                        << ", min|Jac|"
                        << ", max|Jac|"
                        << ", ratio(|Jac|<=0)";
            }

          if( m_TrueField )
            {
            this->m_Fid << ", dist(warp,true warp)"
                        << ", dist(Jac,true Jac)";
            }

          this->m_Fid << std::endl;
          if( backdeffield.IsNotNull() )
            {
            this->m_Bid << std::endl;
            }

          m_headerwritten = true;
          }

        this->m_Fid << iter
                    << ", " << metricbefore
                    << ", " << harmonicEnergy
                    << ", " << minJac
                    << ", " << maxJac
                    << ", " << jacBelowZeroPrc;

        if( m_TrueField )
          {
          this->m_Fid << ", " << fieldDist
                      << ", " << fieldGradDist;
          }

        this->m_Fid << std::endl;
        }
      }
    if( backdeffield.IsNotNull() )
      {
      std::cout << "  : (B)MSE " << backmetricbefore << " - ";

      double fieldDist = -1.0;
      double fieldGradDist = -1.0;
      double tmp;
//    m_headerwritten = false;

      if( m_TrueField )
        {
        typedef itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType>
          FieldIteratorType;
        FieldIteratorType currIter(
          backdeffield, backdeffield->GetLargestPossibleRegion() );
        FieldIteratorType trueIter(
          m_TrueField, backdeffield->GetLargestPossibleRegion() );

        m_CompWarpGradientCalculator->SetInputImage( backdeffield );

        fieldDist = 0.0;
        fieldGradDist = 0.0;
        for( currIter.GoToBegin(), trueIter.GoToBegin();
             not currIter.IsAtEnd(); ++currIter, ++trueIter )
          {
          fieldDist += ( currIter.Value() - trueIter.Value() ).GetSquaredNorm();

          // No need to add Id matrix here as we do a substraction
          tmp = (
              ( m_CompWarpGradientCalculator->EvaluateAtIndex( currIter.GetIndex() )
                - m_TrueWarpGradientCalculator->EvaluateAtIndex( trueIter.
                                                                 GetIndex() )
              ).GetVnlMatrix() ).frobenius_norm();
          fieldGradDist += tmp * tmp;
          }
        fieldDist = sqrt( fieldDist / (double)(
                            backdeffield->GetLargestPossibleRegion().
                            GetNumberOfPixels() ) );
        fieldGradDist = sqrt( fieldGradDist / (double)(
                                backdeffield->GetLargestPossibleRegion().
                                GetNumberOfPixels() ) );

        std::cout << "d(.,true) " << fieldDist << " - ";
        std::cout << "d(.,Jac(true)) " << fieldGradDist << " - ";
        }
      m_HarmonicEnergyCalculator->SetImage( backdeffield );
      m_HarmonicEnergyCalculator->Compute();
      const double backharmonicEnergy
        = m_HarmonicEnergyCalculator->GetHarmonicEnergy();
      std::cout << "harmo. " << backharmonicEnergy << " - ";

      m_JacobianFilter->SetInput( backdeffield );
      m_JacobianFilter->UpdateLargestPossibleRegion();

      const unsigned int numPix = m_JacobianFilter->
        GetOutput()->GetLargestPossibleRegion().
        GetNumberOfPixels();

      TPixel *pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
      TPixel *pix_end = pix_start + numPix;

      TPixel *jac_ptr;

      // Get percentage of det(Jac) below 0
      unsigned int jacBelowZero(0u);
      for( jac_ptr = pix_start; jac_ptr != pix_end; ++jac_ptr )
        {
        if( *jac_ptr <= 0.0 )
          {
          ++jacBelowZero;
          }
        }
      const double jacBelowZeroPrcb = static_cast<double>( jacBelowZero )
        / static_cast<double>( numPix );

      // Get min an max jac
      const double minJacb = *( std::min_element(pix_start, pix_end) );
      const double maxJacb = *( std::max_element(pix_start, pix_end) );

      std::cout << "max|Jac| " << maxJacb << " - "
                << "min|Jac| " << minJacb << " - "
                << "ratio(|Jac|<=0) " << jacBelowZeroPrcb << std::endl;

      if( this->m_Bid.is_open() )
        {
        if( not m_headerwritten )
          {
          this->m_Bid << "Iteration"
                      << ", MSE before"
                      << ", Harmonic energy"
                      << ", min|Jac|"
                      << ", max|Jac|"
                      << ", ratio(|Jac|<=0)";

          if( m_TrueField )
            {
            this->m_Bid << ", dist(warp,true warp)"
                        << ", dist(Jac,true Jac)";
            }

          this->m_Bid << std::endl;

          m_headerwritten = true;
          }

        this->m_Bid << iter
                    << ", " << backmetricbefore
                    << ", " << backharmonicEnergy
                    << ", " << minJacb
                    << ", " << maxJacb
                    << ", " << jacBelowZeroPrcb;

        if( m_TrueField )
          {
          this->m_Bid << ", " << fieldDist
                      << ", " << fieldGradDist;
          }

        this->m_Bid << std::endl;
        }
      }
  }

protected:
  CommandIterationUpdate() :
    m_Fid( "f_metricvalues.csv" ),
    m_Bid( "b_metricvalues.csv" ),
    m_headerwritten(false)
  {
    m_JacobianFilter = JacobianFilterType::New();
    m_JacobianFilter->SetUseImageSpacing( true );
    m_JacobianFilter->ReleaseDataFlagOn();

    m_Minmaxfilter = MinMaxFilterType::New();

    m_HarmonicEnergyCalculator = HarmonicEnergyCalculatorType::New();

    m_TrueField = 0;
    m_TrueWarpGradientCalculator = 0;
    m_CompWarpGradientCalculator = 0;
  }

  ~CommandIterationUpdate()
  {
    this->m_Fid.close();
    this->m_Bid.close();
  }
};

/*=========================================================================
 *
 *  Copyright SINAPSE: Scalable Informatics for Neuroscience, Processing and Software Engineering
 *            The University of Iowa
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include <sstream>

#include <itkMath.h>
#include <itkPoint.h>
#include <itkVersor.h>

itk::Matrix<double, 3, 3>
TestCreateRotationMatrixFromAngles(const double alpha, const double beta, const double gamma)
{
  // alpha is rotate the X axis -- Attitude
  // beta is rotate the Y axis  -- Bank
  // gamma is rotate the Z axis -- Heading
  const double ca = std::cos(alpha);
  const double sa = std::sin(alpha);
  const double cb = std::cos(beta);
  const double sb = std::sin(beta);
  const double cg = std::cos(gamma);
  const double sg = std::sin(gamma);

  itk::Matrix<double, 3, 3> R;

  R(0, 0) = cb * cg;
  R(0, 1) = -ca * sg + sa * sb * cg;
  R(0, 2) = sa * sg + ca * sb * cg;
  R(1, 0) = cb * sg;
  R(1, 1) = ca * cg + sa * sb * sg;
  R(1, 2) = -sa * cg + ca * sb * sg;
  R(2, 0) = -sb;
  R(2, 1) = sa * cb;
  R(2, 2) = ca * cb;
  itk::Matrix<double, 3, 3>::InternalMatrixType test = R.GetVnlMatrix() * R.GetTranspose();
  if (!test.is_identity(1.0e-10))
  {
    std::cout << "Computed matrix is not orthogonal!!!" << std::endl;
    std::cout << R << std::endl;
  }
  return R;
}

itk::Versor<double>
TestCreateRotationVersorFromAngles(const double alpha, const double beta, const double gamma)
{
  // http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
  // psi = alpha is rotate the X axis -- Attitude
  // theta= beta is rotate the Y axis  -- Bank
  // phi=  gamma is rotate the Z axis -- Heading
  const double cha = std::cos(alpha * 0.5);
  const double chb = std::cos(beta * 0.5);
  const double chg = std::cos(gamma * 0.5);
  const double sha = std::sin(alpha * 0.5);
  const double shb = std::sin(beta * 0.5);
  const double shg = std::sin(gamma * 0.5);

  vnl_vector_fixed<double, 4> q;
  q[0] = cha * chb * chg + sha * shb * shg;
  q[1] = sha * chb * chg - cha * shb * shg;
  q[2] = cha * shb * chg + sha * chb * shg;
  q[3] = cha * chb * shg - sha * shb * chg;

  itk::Versor<double> v;
  v.Set(q[1], q[2], q[3], q[0]);
  return v;
}

/**
 * This test that the conversion to and from Rotaion Matrix and
 * Versor produces consistent results.
 */
int
RotationMatrixToVersorTest(void)
{
  int errorCount = 0;
  // const double onedegree=1e-10*itk::Math::pi/180.0;
  const double onedegree = 1e-2 * itk::Math::pi / 180.0;
  // const double td=180.0/itk::Math::pi;
  double centers[6];

  centers[0] = 0;
  centers[1] = itk::Math::pi * 0.25;
  centers[2] = itk::Math::pi * 0.5;
  centers[3] = itk::Math::pi;
  centers[4] = itk::Math::pi * 1.5;
  centers[5] = itk::Math::pi * 2.0;

  constexpr double steps = 5;
  const double     small_degree_steps = onedegree / 1000.0; // 1/1000 of a degree
  for (int j = 0; j < 6; j++)
  {
    for (double alpha = centers[j] - steps * small_degree_steps; alpha <= centers[j] + steps * small_degree_steps;
         alpha += small_degree_steps)
    {
      for (double beta = centers[j] - steps * small_degree_steps; beta <= centers[j] + steps * small_degree_steps;
           beta += small_degree_steps)
      // double beta=0.0;
      {
        for (double gamma = centers[j] - steps * small_degree_steps; gamma <= centers[j] + steps * small_degree_steps;
             gamma += small_degree_steps)
        // double gamma=itk::Math::pi;
        {
          itk::Matrix<double, 3, 3> MR = TestCreateRotationMatrixFromAngles(alpha, beta, gamma);
          itk::Versor<double>       VR = TestCreateRotationVersorFromAngles(alpha, beta, gamma);

          itk::Point<double, 3> testPoint;
          testPoint[0] = -1020.27;
          testPoint[1] = 3.21;
          testPoint[2] = 1000.786432;

          itk::Versor<double> VFROMMR;
          VFROMMR.Set(MR);
          const itk::Point<double, 3> newMRtestPoint = (MR)*testPoint;
          const itk::Point<double, 3> newVRtestPoint = (VR.GetMatrix()) * testPoint;

          const itk::Point<double, 3> newVRFROMMRPoint = (VFROMMR.GetMatrix()) * testPoint;
          const itk::Point<double, 3> newVRFROMMRTransformPoint = VFROMMR.Transform(testPoint);

          const double error_newMRtestPoint_newVRtestPoint = (newMRtestPoint - newVRtestPoint).GetNorm();
          const double error_newMRtestPoint_newVRFROMMRPoint = (newMRtestPoint - newVRFROMMRPoint).GetNorm();
          const double error_newVRFROMMRPoint_newVRFROMMRTransformPoint =
            (newVRFROMMRPoint - newVRFROMMRTransformPoint).GetNorm();

          const double maxAllowedPointError = 1e-5;
          if ((error_newMRtestPoint_newVRtestPoint + error_newMRtestPoint_newVRFROMMRPoint +
               error_newVRFROMMRPoint_newVRFROMMRTransformPoint) > maxAllowedPointError)
          {
            std::cout << "(alpha,beta,gamma)= (" << alpha << "," << beta << "," << gamma << ")" << std::endl;
            std::cout << newMRtestPoint << " " << newVRtestPoint << " " << newVRFROMMRPoint << " "
                      << newVRFROMMRTransformPoint << std::endl;
            std::cout << "ERRORS: " << error_newMRtestPoint_newVRtestPoint << " "
                      << error_newMRtestPoint_newVRFROMMRPoint << " "
                      << error_newVRFROMMRPoint_newVRFROMMRTransformPoint << std::endl;
            std::cout << "MR=\n"
                      << MR << "\nVR=\n"
                      << VR.GetMatrix() << "\nVFROMMR=\n"
                      << VFROMMR.GetMatrix() << std::endl;
            ++errorCount;
          }
        }
      }
    }
  }
  return errorCount;
}

int
main(int, char *[])
{
  return RotationMatrixToVersorTest();
}

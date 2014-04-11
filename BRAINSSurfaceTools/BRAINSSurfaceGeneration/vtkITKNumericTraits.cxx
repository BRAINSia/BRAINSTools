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

namespace itk
{
#if defined(VTK_TYPE_USE___INT64)
const __int64 NumericTraits<__int64 >::Zero = 0;
const __int64 NumericTraits<__int64 >::One = 1;

const unsigned __int64 NumericTraits<unsigned __int64 >::Zero = 0ui64;
const unsigned __int64 NumericTraits<unsigned __int64 >::One = 1ui64;
#endif

#if ITK_VERSION_MAJOR < 4
#if defined(VTK_TYPE_USE_LONG_LONG) && !defined(ITK_TYPE_USE_LONG_LONG)
const long long NumericTraits<long long >::Zero = 0;
const long long NumericTraits<long long >::One = 1;

const unsigned long long NumericTraits<unsigned long long >::Zero = 0;
const unsigned long long NumericTraits<unsigned long long >::One = 1;
#endif
#endif
};

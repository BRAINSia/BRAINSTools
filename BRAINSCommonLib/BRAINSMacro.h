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
#ifndef __BRAINSMacro_h
#define __BRAINSMacro_h

/**
  * This class encapulates all the functionality needed to run the BRAINSFit
  *program without any command line options.
  *
  * It is a class that has all the functionality after the command line options
  *are processed and the images pre-processed, and
  * returns all binary versions of the objects for post-processing.
  *
  * NOTE:  This class is not templated!
  */

#if defined( NDEBUG )
#define VECTORitkDebugMacro(s, type, x)
#else
#define VECTORitkDebugMacro(s, type, x)                                 \
    {                                                                   \
    if( this->GetDebug() && ::itk::Object::GetGlobalWarningDisplay() )  \
      {                                                                 \
      ::std::ostringstream itkmsg;                                      \
      itkmsg << "Debug: In " __FILE__ ", line " << __LINE__ << "\n"     \
             << this->GetNameOfClass() << " (" << this << "):" << s;    \
      for( type::const_iterator it = x.begin(); it != x.end(); ++it )   \
        {                                                               \
        itkmsg << " "                                                   \
               << *( it );                                              \
        }                                                               \
      itkmsg << "\n\n";                                                 \
      ::itk::OutputWindowDisplayDebugText( itkmsg.str().c_str() );      \
      }                                                                 \
    }
#endif
/** Set built-in type.  Creates member Set"name"() (e.g., SetVisibility()); */
#define VECTORitkSetMacro(name, type)                           \
  virtual void Set##name(const type &_arg)                     \
    {                                                           \
    VECTORitkDebugMacro("setting " #name " to ", type, _arg);   \
      {                                                         \
      this->m_##name.resize( _arg.size() );                     \
      this->m_##name = _arg;                                    \
      this->Modified();                                         \
      }                                                         \
    }

/** Get built-in type.  Creates member Get"name"() (e.g., GetVisibility()); */
#define VECTORitkGetConstMacro(name, type)                                      \
  virtual const type &Get##name() const                                        \
    {                                                                           \
    VECTORitkDebugMacro("returning " << #name " of ", type, this->m_##name);    \
    return this->m_##name;                                                      \
    }

#endif // __BRAINSMACRO_h

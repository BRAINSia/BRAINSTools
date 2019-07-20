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
#ifndef _TRON_H
#define _TRON_H

class function
{
public:
  virtual double
  fun( double * w ) = 0;

  virtual void
  grad( double * w, double * g ) = 0;

  virtual void
  Hv( double * s, double * Hs ) = 0;

  virtual int
  get_nr_variable( void ) = 0;

  virtual ~function( void ) {}
};

class TRON
{
public:
  TRON( const function * fun_obj, double eps = 0.1, int max_iter = 1000 );
  ~TRON();

  void
  tron( double * w );

  void
  set_print_string( void ( *i_print )( const char * buf ) );

private:
  int
  trcg( double delta, double * g, double * s, double * r );

  double
  norm_inf( int n, double * x );

  double     eps;
  int        max_iter;
  function * fun_obj;
  void
  info( const char * fmt, ... );

  void ( *tron_print_string )( const char * buf );
};
#endif

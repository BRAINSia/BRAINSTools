/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    ImageCalculatorUtils.h
  Language:  C++

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#if !defined(__ImageCalculator_h__)
#define __ImageCalculator_h__
#include <iostream>
#include <string>

// This function prints the valid pixel types.
extern void PrintDataTypeStrings(void);

// This function replaces a substring with another substring
extern void ReplaceSubWithSub(std::string& s, const char *o, const char  *n);

// This function compares strings.
extern int CompareNoCase( const std::string & s, const std::string& s2 );

extern int PrimaryImageCalculatorRoutine(int argc, char *argv[]);

#endif // __ImageCalculator_h__

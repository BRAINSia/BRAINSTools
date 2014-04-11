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

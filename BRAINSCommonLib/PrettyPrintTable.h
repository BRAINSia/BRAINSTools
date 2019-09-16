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
#ifndef __PrettyPrintTable_h
#define __PrettyPrintTable_h
#include <cstdio>
#include <vector>
#include <string>
#include <ostream>
#include "itkNumberToString.h"
#if defined(_WIN32) || defined(_WIN64)
// Windows uses a different function name for this behavior.
#  define SNPRINTF_FUNC _snprintf
#else
#  define SNPRINTF_FUNC snprintf
#endif
/**
 * \class PrettyPrintTable
 * \author Kent Williams
 * Simple class to print out column-aligned tables
 */
class PrettyPrintTable
{
public:
private:
  using rowType = std::vector<std::string>;
  using tableType = std::vector<rowType>;

  tableType    m_Table;
  unsigned int m_Pad;
  bool         m_rightJustify;

public:
  PrettyPrintTable()
    : m_Pad(1)
    , m_rightJustify(false)
  {}

  void
  setTablePad(unsigned int pad)
  {
    this->m_Pad = pad;
  }

  void
  leftJustify(void)
  {
    m_rightJustify = false;
  }

  void
  rightJustify(void)
  {
    m_rightJustify = true;
  }

  void
  add(const unsigned int row, const unsigned int column, const char * const s)
  {
    // Make sure the table has enough rows.
    if (m_Table.size() <= row)
    { // Add empty rows
      m_Table.resize(row + 1);
    }
    // For each row, make sure that it now has enough columns.
    for (unsigned int q = 0; q < m_Table.size(); ++q)
    {
      if (m_Table[q].size() <= column)
      {
        m_Table[q].resize(column + 1, std::string(""));
      }
    }
    m_Table[row][column] = s;
  }

  void
  add(const unsigned int row, const unsigned int column, const std::string & s)
  {
    add(row, column, s.c_str());
  }

  void
  add(const unsigned int row, const unsigned int column, const int x, const char * printf_format = nullptr)
  {
    const char * format(printf_format == nullptr ? "%d" : printf_format);
    char         buf[4096];

    SNPRINTF_FUNC(buf, 4096, format, x);
    this->add(row, column, buf);
  }

  void
  add(const unsigned int row, const unsigned int column, const unsigned int x, const char * printf_format = nullptr)
  {
    const char * format(printf_format == nullptr ? "%d" : printf_format);
    char         buf[4096];

    SNPRINTF_FUNC(buf, 4096, format, x);
    this->add(row, column, buf);
  }

  void
  add(const unsigned int row, const unsigned int column, const double x, const char * printf_format = nullptr)
  {
    if (printf_format != nullptr)
    {
      char buf[4096];
      SNPRINTF_FUNC(buf, 4096, printf_format, x);
      this->add(row, column, buf);
    }
    else
    {
      itk::NumberToString<double> doubleToString;
      std::string                 val = doubleToString(x);
      this->add(row, column, val.c_str());
    }
  }

  void
  Print(std::ostream & output)
  {
    using ColWidthsType = std::vector<unsigned int>;
    ColWidthsType colWidths(m_Table[0].size(), 0);
    // find largest columns
    for (unsigned i = 0; i < m_Table.size(); ++i)
    {
      for (unsigned j = 0; j < m_Table[i].size(); ++j)
      {
        if (colWidths[j] < m_Table[i][j].size())
        {
          colWidths[j] = m_Table[i][j].size();
        }
      }
    }
    for (unsigned i = 0; i < m_Table.size(); ++i)
    {
      for (unsigned j = 0; j < m_Table[i].size(); ++j)
      {
        // if right justify, output leading blanks
        if (m_rightJustify)
        {
          int count = colWidths[j] - m_Table[i][j].size();
          while (count--)
          {
            output << " ";
          }
        }
        unsigned int k(0);
        for (k = 0; k < m_Table[i][j].size(); ++k)
        {
          output << m_Table[i][j][k];
        }
        unsigned int limit;
        // if right justify, just output pad
        if (m_rightJustify)
        {
          limit = this->m_Pad;
          k = 0;
        }
        else
        {
          // print column fill + pad
          limit = colWidths[j] + this->m_Pad;
        }
        for (; k < limit; ++k)
        {
          output << " ";
        }
      }
      output << std::endl;
    }
  }
};

#endif // PrettyPrintTable_h

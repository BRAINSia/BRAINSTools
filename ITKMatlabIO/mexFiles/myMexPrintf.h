/**
 * Goal: refactor itkLoadWithMetadata and itkSaveWithMetadata
 * Changes: move all myMexPrintf to single header file
 * Data: Oct 13th 16
 * Author: Hui Xie
 */

#ifndef BRAINSTOOLS_MYMEXPRINTF_H
#define BRAINSTOOLS_MYMEXPRINTF_H

bool debug = false;

void
myMexPrintf(std::string msg)
{
  if (debug)
    mexPrintf(msg.c_str());
}

void
myMexPrintf(std::string msg, int value)
{
  if (debug)
    mexPrintf(msg.c_str(), value);
}

void
myMexPrintf(std::string msg, std::string value)
{
  if (debug)
    mexPrintf(msg.c_str(), value.c_str());
}

void
myMexPrintf(std::string msg, double value)
{
  if (debug)
    mexPrintf(msg.c_str(), value);
}

void
myMexPrintf(std::string msg, unsigned int value)
{
  if (debug)
    mexPrintf(msg.c_str(), value);
}

void
myMexPrintf(std::string msg, int value, int value2, int value3)
{
  if (debug)
    mexPrintf(msg.c_str(), value, value2, value3);
}

void
myMexPrintf(std::string msg, double value, int value2, int value3)
{
  if (debug)
    mexPrintf(msg.c_str(), value, value2, value3);
}

#endif // BRAINSTOOLS_MYMEXPRINTF_H

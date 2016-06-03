//
// Created by Leinoff, Alexander on 6/3/16.
//

#ifdef WIN32
#define MODULE_IMPORT __declspec(dllimport)
#else
#define MODULE_IMPORT
#endif

extern "C" MODULE_IMPORT int ModuleEntryPoint(int, char * []);

int BRAINSDefaceTest(int argc, char* argv[])
{
  return ModuleEntryPoint(argc, argv);
}
#ifndef HAPRAD_ERROR
#define HAPRAD_ERROR

// Pre-processor macro to report messages

#include <iostream>
#define HAPRAD_OS std::cerr

// giving a location + string

#define HAPRAD_INFO_MSG(loc,str) \
   HAPRAD_OS << "Info in " << loc << ": "  << str \
             << std::endl;
#define HAPRAD_WARN_MSG(loc,str) \
   HAPRAD_OS << "Warning in " << loc << ": " << str \
             << std::endl;
#define HAPRAD_ERROR_MSG(loc,str) \
   HAPRAD_OS << "Error in " << loc << ": " << str \
             << std::endl;

// giving location + a value

# define HAPRAD_INFO_VAL(loc,x) \
   HAPRAD_OS << "Info in " << loc << ": " << #x << " = " << (x) \
             << std::endl;
# define HAPRAD_WARN_VAL(loc,x) \
   HAPRAD_OS << "Warning in " << loc << ": " << #x << " = " << (x) \
             << std::endl;
# define HAPRAD_ERROR_VAL(loc,x) \
   HAPRAD_OS << "Error in " << loc << ": " << #x << " = " << (x) \
             << std::endl;

// giving a location + string + value

# define HAPRAD_INFO_MSGVAL(loc,str,x) \
   HAPRAD_OS << "  Info in " << loc << ": " << str \
             << "; " << #x << " = " << (x) << std::endl;
# define HAPRAD_WARN_MSGVAL(loc,str,x) \
   HAPRAD_OS << "Warning in " << loc << ": " << str \
             << "; " << #x << " = " << (x) << std::endl;
# define HAPRAD_ERROR_MSGVAL(loc,str,x) \
   HAPRAD_OS << "Error in " << loc << ": " << str \
             << "; " << #x << " = " << (x) << std::endl;

#endif

#!/bin/sh
set -e
if test "$CONFIGURATION" = "Debug"; then :
  cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode
  make -f /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "Release"; then :
  cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode
  make -f /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "MinSizeRel"; then :
  cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode
  make -f /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode/CMakeScripts/ReRunCMake.make
fi
if test "$CONFIGURATION" = "RelWithDebInfo"; then :
  cd /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode
  make -f /Users/giuliabenvenuto/Documents/GitHub/SSGD_Problem/build_xcode/CMakeScripts/ReRunCMake.make
fi


% tested under Windows7 64bit, MSVC v10.0 compiler
mex OPTIMFLAGS="/Ox /Oi /Oy /DNDEBUG /fp:fast /arch:SSE2 /DMEX_MODE" mexFGS.cpp
mex OPTIMFLAGS="/Ox /Oi /Oy /DNDEBUG /fp:fast /arch:SSE2 /DMEX_MODE" mexFGS_simple.cpp
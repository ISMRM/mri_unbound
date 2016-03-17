REAME.txt

Author: Nick Zwart, Dallas Turley
Date: 2011 apr 11
Rev: 2011 aug 22

The 3D grid algorithm implemented in this package is 
located in the file: 

    grid_utils.c (main subroutine: grid3())
    
Wrapper (or gateway) functions have been provided for use 
in MATLAB and AVS.  However, the grid_utils.c code
is independent of either software and may be included as
a support routine in other C/C++ code.  The threaded code
requires the POSIX-Threads (pthread.h, not included) library and an 
additional support file threads.c (included).

OPTIONAL NON-THREADED COMPILATION:
    If pthread.h is not available due to compatibility issues, delete
    the functions _addPartition3_thread() and grid3_threaded() 
    and implement grid3() directly.  The MATLAB and AVS wrappers include
    calls to the threaded functions, so those calls will have to be 
    removed as well in order to compile a totally non-threaded binary.  

MATLAB:

    The MATLAB mex-wrapper implements the grid3() and
    grid3_threaded() subroutines in:

        grid3_MAT.c

    This is compiled by navigating to the directory 
    containing this package, and typing (LINUX or OSX):
        
        make -f Makefile.MATLAB grid3_MAT

    from a shell command-line.

    This has been tested with the following configurations:

        LINUX:
            Ubuntu 10.10 x86_64
            gcc version 4.4.5 (Ubuntu/Linaro 4.4.4-14ubuntu5)
            MATLAB 7.9.0.529 (R2009b) 64bit

    The MEX file may be tested using the supplied data by
    running the included testmex script.  This script
    executes the grid3_MAT() command on supplied coordinates,
    simulated data, and density compensation data, and 
    compares the output to a provided solution set.  

AVS:

    The AVS wrapper implements grid3() and 
    grid3_threaded() in:

        grid3_AVS.cpp

    The binary is built using the make-file:

        make -f Makefile.AVS grid3_AVS

    Tested platforms:

        LINUX:
            Ubuntu 10.10 i386
            gcc version 4.4.3 (Ubuntu 4.4.3-4ubuntu5)
            AVS version: 5.6 (50.96 i686 ogl RedHat9.0)



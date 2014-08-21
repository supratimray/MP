Matching Pursuit ported to 32 buit Windows, using Visual Studio C++ 2010 Express


BK (bart@vision.rutgers.edu) - May 2013.

***INSTALL:***

The source/ directory now contains a gabord.exe which the Matlab wrappers will call automatically. (See runGabor.m). 
Nothing else needs to be done. 


***CHANGES:***

1. Linked with Winsock2 (ws2_32.lib) to get gethostbyname 
2. Defined some built-in VC defines to get rid of some of the Windows predefines: VC_EXTRALEAN;WIN32_LEAN_AND_MEAN;NOGDI
3. Defined WINDOWS and added #ifdef's in the code to generate OS specific code (mainly includes). 
	All Windows specific changes are hidden behind #ifdef WINDOWS. 	
4. Added a strndup function.
5. Changes in the base code (i.e. not behind #ifdef WINDOWS)
	a. MSVC++ does not initialized longs to zero. Probably gcc does? I explicitly set this now in update.c for those variables that needed it in the test run.
		Of course there could be other uninitalized variables that my tests did not use... 
	b. Renamed 'WORD' in the original MP code to  'WRD' and 'word' to 'wrd' because Windows uses WORD extensively for a very different type.This is a safe change, just renaming.

Compiled/Linked fine for Win32 in Debug and Release.
Excutables are written to ../source/gabord.exe to be consistent with the gcc output.

Tested with exampleMP.m 
	
Likely issues: 
	Spaces in the file name specification on the command line.
	The presence of ':' in the specified paths in the CTL file. 

Solution:
	I rewrote the Matlab wrapper code slightly to work in the directory that contains the ctl and data.


***COMPILE:***

Install Visual Studio C++ 2010 (The free express edition works fine,and later versions should be backward compatible).
Open MP.sln in VC++
Choose Release version.
Build. 

For debugging, choose the Debug version which will create gabord.debug.exe in the source directory.

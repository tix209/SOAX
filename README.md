 # SOAX
A software for 2D/3D biopolymer network extraction and quantification.

## About
SOAX is an open source software tool to extract the centerlines, junctions and filament lengths of biopolymer networks in 2D and 3D images. It facilitates quantitative, reproducible and objective analysis of the image data. The underlying method of SOAX uses multiple Stretching Open Active Contours (SOACs) that are automatically initialized at image intensity ridges and then stretch along the centerlines of filaments in the network. SOACs can merge, stop at junctions, and reconfigure with others to allow smooth crossing at junctions of filaments.

SOAX provides 3D visualization for exploring image data and visually checking results against the image. Quantitative analysis functions based on extracted networks are also implemented in SOAX, including spatial distribution, orientation, and curvature of filamentous structures. SOAX also provides interactive manual editing to furthure improve the extraction results, which can be saved in a file for archiving or further analysis.

**\*Note\*** [TSOAX](https://github.com/tix209/TSOAX) is an extension of SOAX to support tracking filaments and networks over multiple frames.

## Reference

T. Xu, D. Vavylonis, F. Tsai, G. H. Koenderink, W. Nie, E. Yusuf, I-Ju Lee, J.-Q. Wu, and X. Huang, "[SOAX: A Software for Quantification of 3D Biopolymer Networks](http://www.nature.com/srep/2015/150313/srep09081/full/srep09081.html)", In Scientific Reports, 5, 9081; DOI:10.1038/srep09081, 2015.

T. Xu, D. Vavylonis, X. Huang, "[3D Actin Network Centerline Extraction with Multiple Active Contours](http://www.sciencedirect.com/science/article/pii/S136184151300159X)", In Medical Image Analysis, 18(2):272-84, 2014.


## Acknowledgement

This work has been supported by NIH grants R01GM098430 and R35GM136372.

## Sample Data

A synthetic test image and ground truth are available:
- [Test image](https://www.lehigh.edu/~div206/soax/data/cable1_std_4.mha)
- [Ground truth](https://www.lehigh.edu/~div206/soax/data/cable1.txt). Load the ground truth via *Load JFilament Snakes* in SOAX's *File* menu.



## Building from Source

SOAX 3.8.x has been upgraded to support: Qt 5, VTK 8.2+ or 9.x, ITK 4.x or 5.x, CMake 3.16+, Boost 1.70+, C++17 compatible compiler

### Ubuntu (Tested on Ubuntu 22.04)

##### Install Dependencies

sudo apt-get install build-essential libeigen3-dev qtbase5-dev qtcreator qt5-qmake libqt5x11extras5-dev cmake x11-apps cmake-gui libxt-dev libboost-all-dev

##### Compile VTK from Source

Get VTK 8.2.0 [https://github.com/Kitware/VTK/tree/v8.2.0](https://github.com/Kitware/VTK/tree/v8.2.0) or later
and onfigure with cmake-gui or cmake

Set VTK_Group_Qt: ON, BUILD_SHARED_LIBS: OFF, CMAKE_BUILD_TYPE: Release

Enable VTK_Group_Qt and Module_vtkGUISupportQtOpenGL

Possible compiling errors: Add ```#include <QPainterPath>``` in VTK-8.2.0\Rendering\Qt\vtkQtLabelRenderStrategy.h and VTK-8.2.0\Rendering\Qt\vtkQtStringToImage.h

In VTK build directory:  make -j 4

##### Compile ITK from Source 
Download ITK 5.x from https://github.com/InsightSoftwareConsortium/ITK/releases

Configure with cmake or cmake-gui

Set BUILD_SHARED_LIBS: OFF, CMAKE_BUILD_TYPE: Release, 
BUILD_STATIC_LIBS: ON, Module_ITKFEM: ON, Module_ITKVtkGlue: ON

Build: make -j 4


##### Compile SOAX

git clone https://github.com/tix209/SOAX.git

Configure with cmake or cmake-gui:
Set ITK_DIR to your ITK build path

Set CMAKE_BUILD_TYPE: Release

In SOAX build directory: make -j 4




## Build instructions for versions 3.7 and earlier

Library dependencies: Qt4, VTK 6.x, ITK 4.x

### Windows, Visual Studio 2010, Static linking (Qt4/VTK6/ITK4)
<details>
<summary>Click to expand build instructions</summary>

1. Install fully up-to-date Visual Studio 2010 with SP1
2. Install Perl (32-bit)
3. Install Python 2.7.x
4. Install CMake
5. Download and extract Qt source 4.8.6 to C:\Qt
6. Setup environment variables. In Visual Studio Command Prompt, type
set QMAKESPEC=win32-msvc2010
set QTDIR=C:\Qt
set PATH=C:\Qt\bin;%PATH%
7. Download jom and extract it to C:\Qt\jom
8. Open qmake.conf in %QTDIR%\mkspecs\win32-msvc2010 and change the line
QMAKE_CFLAGS_RELEASE = /O2 /MD
to
QMAKE_CFLAGS_RELEASE = /O2 /MT
9. Do the same for QMAKE_CFLAGS_DEBUG and QMAKE_CFLAGS_RELEASE_WITH_DEBUGINFO;
Start VS2010 command prompt (for 64-bit, use Visual Studio x64 Win64 Command Prompt) and run
configure -debug-and-release -static -opensource -no-phonon -no-phonon-backend -no-qt3support -no-multimedia -no-webkit -no-libtiff
then
.\jom\jom.exe -j N
where N is the number of CPU cores you want to use for Qt compilation.
10. Download boost library and extract it to C:\
11. Go to boost directory and type
bootstrap
then
.\b2 msvc link=static runtime-link=static threading=multi stage
12. Download VTK source and configure using CMake. In CMake,
Change all /MD to /MT in CMAKE_CXX_FLAGS_MINSIZEREL, CMAKE_CXX_FLAGS_RELEASE, CMAKE_CXX_FLAGS_RELWITHDEBINFO, CMAKE_C_FLAGS_MINSIZEREL, CMAKE_C_FLAGS_RELEASE, and CMAKE_C_FLAGS_RELWITHDEBINFO;
Change /MDd to /MTd in CMAKE_CXX_FLAGS_DEBUG and CMAKE_C_FLAGS_DEBUG;
Uncheck BUILD_SHARED_LIBS;
Check Module_vtkGUISupportQt, Module_vtkGUISupportQtOpenGL, Module_vtkRenderingQt, Module_vtkViewsQt;
13. Build VTK in Debug and Release mode using VS2010
14. Download ITK source and configure using CMake. In CMake,
Change all /MD to /MT and /MDd to /MTd as in configuring VTK;
Check Module_ItkVtkGlue;
15. Build ITK in Debug and Release mode using VS2010
16. Download the above SOAX source and configure using CMake. In CMake,
Change all /MD to /MT and /MDd to /MTd
Set Boost_INCLUDE_DIR to the topmost boost directory
Set Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE to "C:/boost_1_55_0/stage/lib/libboost_program_options-vc100-mt-s-1_55.lib"
Set Boost_PROGRAM_OPTIONS_LIBRARY_DEBUG to "C:/boost_1_55_0/stage/lib/libboost_program_options-vc100-mt-sgd-1_55.lib"
Do the same for Boost_FILESYSTEM_LIBRARY_DEBUG, Boost_FILESYSTEM_LIBRARY_RELEASE, Boost_SYSTEM_LIBRARY_DEBUG and Boost_SYSTEM_LIBRARY_RELEASE
17. Build SOAX in Debug and Release using VS2010


</details>

### Building SOAX in Mac OS X
<details>
<summary>Click to expand build instructions</summary>

1. Install Homebrew
2. Install Boost C++ libraries:
brew install boost
3. Install Qt 4 libraries (http://github.com/cartr/homebrew-qt4/)
using:
```
brew tap cartr/qt4
brew tap-pin cartr/qt4
brew install qt@4
brew install qt-webkit@2.3
```
4. Install CMake
5. Configure and build VTK 6.3 from source (Enable VTK_Group_Qt and change CMAKE_BUILD_TYPE to "Release" in CMake; Disable BUILD_SHARED_LIBS; build in Xcode by opening "VTK.xcodeproj" generated by CMake and then Product > Build)
6. Configure and build ITK 4.7.2 from source (Disable BUILD_EXAMPLES, BUILD_TESTING, and enable Module_ITKVtkGlue in CMake; build Release version in Xcode)
7. Download SOAX source code above, configure by CMake and compile a Release version in Xcode
8. Use macdeployqt for making a dmg application
</details>

### Building SOAX in Windows Linux Subsystem (Ubuntu 18.04.1 LTS)
<details>
<summary>Click to expand build instructions</summary>

David Rutkowski, 12/2020

1. Newer version of cmake: https://askubuntu.com/questions/829310/how-to-upgrade-cmake-in-ubuntu
ccmake:
sudo apt-get install cmake-curses-gui libxt-dev
2. OpenGL implementation on Linux (if this is required):
sudo apt-get install libegl1-mesa-dev
https://stackoverflow.com/questions/31170869/cmake-could-not-find-opengl-in-ubuntu
3. Subsystem GUI:
Install XMing on Windows
sudo apt-get install vim-gtk
export DISPLAY=:0
https://www.howtogeek.com/261575/how-to-run-graphical-linux-desktop-applications-from-windows-10s-bash-shell/
4. QT4 package:
sudo apt install qt4-default
https://askubuntu.com/questions/1102660/enable-to-locate-qt-sdk-package
sudo apt install build-essential cmake gdb git libphonon-dev libqt4-dev libqt4-opengl-dev libqtwebkit-dev qt4-designer qt4-doc qt4-doc-html qt4-qmake qtcreator qtcreator-doc subversion
5. Boost libraries:
sudo apt-get install libboost-all-dev
6. VTK 6.3.0:
https://github.com/Kitware/VTK/releases?after=v7.1.0
Use ccmake from a separate build folder to Enable VTK_Group_Qt and change CMAKE_BUILD_TYPE to "Release"
7. ITK 4.7.2:
https://sourceforge.net/projects/itk/files/itk/4.7/
Use ccmake to Disable BUILD_EXAMPLES, BUILD_TESTING; Enable Module_ITKVtkGlue
8. Add the following to vcl_compiler.h (for gcc 7.5.0) (\InsightToolkit-4.7.2\Modules\ThirdParty\VNL\src\vxl\vcl\vcl_compiler.h) at line 129
```
# elif (__GNUC__==7)
# define VCL_GCC_7
# if (__GNUC_MINOR__ > 4 )
# define VCL_GCC_75
# else
# define VCL_GCC_70
# endif
```
9. Change vcl_new.h (\InsightToolkit-4.7.2\Modules\ThirdParty\VNL\src\vxl\vcl\vcl_new.h)
Line 19 to # include < new > from # include < new.h>
10. SOAX:
Ccmake while specifying location of VTK and ITK build folders
(These were already populated this way but just in case they are needed)
CMAKE_INSTALL_PREFIX:/usr/local
QT_IMPORTS_DIR:/usr/lib/x86_64-linux-gnu/qt4/imports
QT_QMAKE_EXECUTABLE:/user/bin/qmake

</details>






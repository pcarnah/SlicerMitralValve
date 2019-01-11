cmake_minimum_required(VERSION 2.8.9)

project(MVModeller)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MVModeller")
set(EXTENSION_CATEGORY "Cardiac")
set(EXTENSION_CONTRIBUTORS "Patrick Carnahan (University of Western Ontario)")
set(EXTENSION_DESCRIPTION "This extension aids the process of modelling a mitral valve from 3d ultrasound data.")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/MVModeller.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/MVModeller/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(BiplaneRegistration)
add_subdirectory(MVSegmenter)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})

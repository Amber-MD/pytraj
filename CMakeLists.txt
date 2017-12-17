# AmberTools CMake buildfile for the pytraj wrapper for cpptraj
# NOTE: This is NOT a standalone build file
# It is used by the Amber build system for when it builds pytraj
# If you are just trying to build pytraj, run setup.py, not this!

project(pytraj)


#figure out arguments

if(OPENMP AND NOT TARGET_OSX) # pytraj will not build with OpenMP on macOS
	set(OPENMP_ARG "")
	set(CPPTRAJ_TARGET "libcpptraj_omp")
else()
	set(OPENMP_ARG --disable-openmp)
	set(CPPTRAJ_TARGET "libcpptraj")
endif()

set(RPATH_ARG "")

if(TARGET_OSX)
	# OS X: set RPATH to path up from pysander's pyd modules to AMBERHOME/lib
	set(RPATH_ARG LDFLAGS="-Wl,-rpath,@loader_path/../../../../..")
endif()

#------------------------------------------------------------------------------------------

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/python-build)
set(STAMP_FILE ${BUILD_DIR}/pytraj-build.stamp)
set(LIBDIR_FILE ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/libcpptraj.location)

set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${BUILD_DIR})

file(MAKE_DIRECTORY ${BUILD_DIR})
#find all python library files as dependencies
file(GLOB_RECURSE PYTHON_SOURCES "*.py" "*.pyx" "*.pxd")

#we use cmake -E env to manipulate the environment of setup.py

add_custom_command(OUTPUT ${STAMP_FILE}
	COMMAND ${CMAKE_COMMAND} -E env --unset=AMBERHOME 
		CPPTRAJ_HEADERDIR=${CMAKE_CURRENT_SOURCE_DIR}/../cpptraj/src CPPTRAJ_LIBDIR=$<TARGET_FILE_DIR:${CPPTRAJ_TARGET}>
		CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ${RPATH_ARG}
		${PYTHON_EXECUTABLE} setup.py build -b ${BUILD_DIR} ${OPENMP_ARG} ${PYTHON_COMPILER_ARG} ${WIN64_DEFINE_ARG}
	COMMAND ${CMAKE_COMMAND} -E touch ${STAMP_FILE}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS ${PYTHON_SOURCES} ${CPPTRAJ_TARGET}
	VERBATIM
	COMMENT "Building pytraj native library")

#We want to build the python library during the build step so as to catch any build errors
add_custom_target(pytraj ALL DEPENDS ${STAMP_FILE} SOURCES ${PYTHON_SOURCES})

#you can't put a generator expression into an install SCRIPT rule, which I need to do.
#waiting on https://gitlab.kitware.com/cmake/cmake/issues/15785
file(GENERATE OUTPUT ${LIBDIR_FILE} CONTENT "
	#this file contains the real location of libcpptraj, which is a Special Thing in CMake.  It is read by the install script.
	set(CPPTRAJ_LIBDIR $<TARGET_FILE_DIR:${CPPTRAJ_TARGET}>)")

install(CODE "
	include(${LIBDIR_FILE})
	${FIX_BACKSLASHES_CMD}
	execute_process(
    COMMAND \"${CMAKE_COMMAND}\" -E env  --unset=AMBERHOME
    \"CPPTRAJ_HEADERDIR=${CMAKE_CURRENT_SOURCE_DIR}/../cpptraj/src\" \"CPPTRAJ_LIBDIR=\${CPPTRAJ_LIBDIR}\"
    \"CXX=${CMAKE_CXX_COMPILER}\"
    \"CC=${CMAKE_C_COMPILER}\"
    ${RPATH_ARG}
    ${PYTHONPATH_SET_CMD}
	${PYTHON_EXECUTABLE} ./setup.py build -b ${BUILD_DIR} ${OPENMP_ARG} ${PYTHON_COMPILER_ARG} ${WIN64_DEFINE_ARG} install ${PYTHON_PREFIX_ARG}
    WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\")")

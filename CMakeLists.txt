# AmberTools CMake buildfile for the pytraj wrapper for cpptraj
# NOTE: This is NOT a standalone build file
# It is used by the Amber build system for when it builds pytraj
# If you are just trying to build pytraj, run setup.py, not this!

project(pytraj)


#figure out arguments

if(OPENMP)
	set(OPENMP_ARG "")
else()
	set(OPENMP_ARG --disable-openmp)
endif()

#------------------------------------------------------------------------------------------

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/python-build)
set(STAMP_FILE ${BUILD_DIR}/pytraj-build.stamp)
set(LIBDIR_FILE ${BUILD_DIR}/libcpptraj.location)

file(MAKE_DIRECTORY ${BUILD_DIR})
#find all python library files as dependencies
file(GLOB_RECURSE PYTHON_SOURCES "*.py" "*.pyx" "*.pxd")

#we use cmake -E env to manipulate the environment of setup.py

add_custom_command(OUTPUT ${STAMP_FILE}
	COMMAND ${CMAKE_COMMAND} -E env --unset=AMBERHOME 
		CPPTRAJ_HEADERDIR=${CMAKE_CURRENT_SOURCE_DIR}/../cpptraj/src CPPTRAJ_LIBDIR=$<TARGET_FILE_DIR:libcpptraj>
		CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
		${PYTHON_EXECUTABLE} setup.py build_ext --no-setuptools -b ${BUILD_DIR} ${OPENMP_ARG} ${PYTHON_COMPILER_ARG} ${WIN64_DEFINE_ARG}
	COMMAND ${CMAKE_COMMAND} -E touch ${STAMP_FILE}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS ${PYTHON_SOURCES} libcpptraj
	VERBATIM
	COMMENT "Building pytraj native library")

#We want to build the python library during the build step so as to catch any build errors
add_custom_target(pytraj ALL DEPENDS ${STAMP_FILE} SOURCES ${PYTHON_SOURCES})

#you can't put a generator expression into an install SCRIPT rule, which I need to do.
#waiting on https://gitlab.kitware.com/cmake/cmake/issues/15785
get_property(LIBCPPTRAJ_INCDIR TARGET libcpptraj PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

file(GENERATE OUTPUT ${LIBDIR_FILE} CONTENT "
	#this file contains the real location of libcpptraj, which is a Special Thing in CMake.  It is read by the install script.
	set(CPPTRAJ_LIBDIR $<TARGET_FILE_DIR:libcpptraj>)")

install(CODE "
	include(${LIBDIR_FILE})
	${FIX_BACKSLASHES_CMD}
	execute_process(
    COMMAND \"${CMAKE_COMMAND}\" -E env  --unset=AMBERHOME
    \"CPPTRAJ_HEADERDIR=${CMAKE_CURRENT_SOURCE_DIR}/../cpptraj/src\" \"CPPTRAJ_LIBDIR=\${CPPTRAJ_LIBDIR}\"
    \"CXX=${CMAKE_CXX_COMPILER}\"
    \"CC=${CMAKE_C_COMPILER}\"
	${PYTHON_EXECUTABLE}
  	./setup.py build_ext --no-setuptools -b \"${BUILD_DIR}\" ${OPENMP_ARG} ${PYTHON_COMPILER_ARG} ${WIN64_DEFINE_ARG} install ${PYTHON_PREFIX_ARG}
    WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\")")
	

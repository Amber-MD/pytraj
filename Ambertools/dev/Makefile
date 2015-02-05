# CPPTRAJ - Main makefile
# Daniel R. Roe
# 2014-12-10

# Create cpptraj binary in ./src/
all: install

# Create cpptraj binary and move to ./bin/
install: config.h
	cd src && $(MAKE) install 

# Create libcpptraj.so
libcpptraj: config.h
	cd src && $(MAKE) libcpptraj

# Run Tests
check: config.h
	cd test && $(MAKE) test

# Clean up
clean: config.h
	cd src && $(MAKE) clean
	cd test && $(MAKE) clean

docs: src/cpptraj.Doxyfile
	cd src && doxygen cpptraj.Doxyfile

# Remove cpptraj binary
uninstall: config.h
	cd src && $(MAKE) uninstall

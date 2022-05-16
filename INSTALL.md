Dependencies
============

* cmake3 : if not already installed, download it at https://cmake.org, or intall it from available packages (rpm, debian, brew, ...)
    
    example (Debian): sudo apt -y install cmake
    
    example (Fedora / Centos): sudo yum install cmake3
    

* C++ compiler

    example (Debian): sudo apt -y install g++
    
    example (Fedora / Centos): sudo yum install gcc-c++
    

* cppunit : A C++ library for unitary tests.

    example (Debian): sudo apt -y install libcppunit-dev
    
    example (Fedora / Centos): sudo yum install cppunit cppunit-devel
    
    
* GIT

    example (Debian): sudo apt -y install git
    
    example (Fedora / Centos): sudo yum install git
    

Download
========

    git clone https://github.com/urgi-anagen/TE_finder.git

Compilation
===========

<TE_finder-root-dir> is the root directory of TE-finder package.

1. go to the <TE_finder-root-dir> directory

2. type 'cmake3 .'

3. type 'make install'

4. binaries are in <TE_finder-root-dir>/bin


blaster:
========

set environment variable in .bashrc:

	BLASTER_BANK_DIRS=/usr/bioinfo/db/blasterdb



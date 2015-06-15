PopAlu
======

PopAlu detects polymorphic Alu elements in high-throughput sequencing data at the population scale.

Dependencies
------------

This program requires:
* SeqAn library to be installed (http://www.seqan.de/)

Addtionally, you might want to install:
* boost library (http://www.boost.org/)
* googletest framework (https://code.google.com/p/googletest/) 

Modify the path of library at file 'Makefile.vars'

    B_LIB = path_of_boost_lib
    G_LIB = path_of_gtest_lib


Compiling the program
---------------------

    make OD=opt3 opt3/alu_usage
    make OD=opt3 opt3/build_dist
    make OD=opt3 opt3/alu_delete
    make OD=opt3 opt3/alu_insert
    make OD=opt3 opt3/unittest_insert

Running the program
-------------------

    opt3/alu_usage build_dist
    opt3/alu_usage alu_delete
    opt3/alu_usage alu_insert

Reference
---------

Please cite:

Qian Y, Kehr B, Halldórsson BV. (2015) PopAlu: population-scale discovery of Alu polymorphisms. PeerJ PrePrints 3:e1430

https://peerj.com/preprints/1174/

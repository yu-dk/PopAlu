## PopAlu
detect polymorphic Alu elements in high-throughput sequencing data at population scale 

## dependency
this program requires:
      SeqAn library to be installed (http://www.seqan.de/)
addtionally, you might want to install:
      boost library (http://www.boost.org/) 
      googletest framework ( https://code.google.com/p/googletest/ ) 

modify the path of library at file 'Makefile.vars'
      B_LIB = path_of_boost_lib
      G_LIB = path_of_gtest_lib

## to compile the program:
make OD=opt3 opt3/alu_usage
make OD=opt3 opt3/build_dist
make OD=opt3 opt3/alu_delete
make OD=opt3 opt3/alu_insert
make OD=opt3 opt3/unittest_insert

## to run the program:
opt3/alu_usage build_dist
opt3/alu_usage alu_delete
opt3/alu_usage alu_insert
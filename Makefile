#-*-makefile-*-
include Makefile.vars

INCLUDE = -I/home/qianyuxx/local/include/
WARN= -W -Wall

CPPFLAGS = $(WARN) $(FLAG) $(INCLUDE) 
CC = g++ -fno-merge-constants -fopenmp
LIB = -lz -lbz2 $(B_LIB)

#######################
## g unit test
GTEST_DIR = ${G_LIB}/gtest-1.7.0
CPPFLAGS += -isystem $(GTEST_DIR)/include # Flags passed to the preprocessor.
CXXFLAGS += -g -Wall -Wextra -pthread # Flags passed to the C++ compiler.

UNITTEST_INSERT_FILES_CPP = unittest_insert.cpp common.cpp utils.cpp insert_utils.cpp 
UNITTEST_INSERT_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(UNITTEST_INSERT_FILES_CPP))

ALU_USAGE_FILES_CPP = alu_usage.cpp common.cpp 
ALU_USAGE_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_USAGE_FILES_CPP))

ALU_DELETE_FILES_CPP = alu_delete.cpp common.cpp utils.cpp delete_utils.cpp
ALU_DELETE_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_DELETE_FILES_CPP))

ALU_INSERT_FILES_CPP = alu_insert.cpp common.cpp utils.cpp insert_utils.cpp 
ALU_INSERT_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_INSERT_FILES_CPP))

ALU_TMP_FILES_CPP = alu_tmp.cpp common.cpp utils.cpp delete_utils.cpp
ALU_TMP_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(ALU_TMP_FILES_CPP))

BUILD_DIST_FILES_CPP = build_dist.cpp common.cpp utils.cpp
BUILD_DIST_FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(BUILD_DIST_FILES_CPP))

FILES_CPP = *cpp *h 
FILES_O = $(patsubst %.cpp,$(OD)/%.o,$(FILES_CPP))

$(OD) :
	mkdir $(OD)

$(OD)/%.o : %.cpp
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

$(OD)/alu_usage: $(OD) $(ALU_USAGE_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_USAGE_FILES_O)

$(OD)/alu_delete: $(OD) $(ALU_DELETE_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_DELETE_FILES_O)

$(OD)/alu_insert: $(OD) $(ALU_INSERT_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(ALU_INSERT_FILES_O)

$(OD)/build_dist: $(OD) $(BUILD_DIST_FILES_O) 
	$(CC) -o $@ $(CPPFLAGS) $(LIB) $(BUILD_DIST_FILES_O)

$(OD)/unittest_insert: $(UNITTEST_INSERT_FILES_O)
	$(CC) -o $@ $(CPPFLAGS) $(CXXFLAGS) ${LIB} -lpthread -L${GTEST_DIR}/lib/.libs -lgtest -lgtest_main $(UNITTEST_INSERT_FILES_O) 

clean:
	/bin/rm -rf $(OD)/* $(addprefix $(OD)/,$(EXE))

TESTS = unittest_insert

# House-keeping build targets.  
all : $(TESTS)
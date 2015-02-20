#include "common.h"

int PrintUsg() {
  cout << "######### software to detect polymorphic Alu elements in paired end sequencing data #########\n"
       << "USAGE: opt3/alu_usage [func name] \n"
       << "## [func name] STRING: one of 3 functions to excute, {build_dist, alu_insert, alu_delete} \n"
       << "Example: opt3/alu_usage build_dist\n"
       << "NB: you need to run build_dist first, before running alu_delete or alu_insert\n"
       << "contact: qianyuxx@gmail.com \n\n";
  return 1;
}

string Parameter1(){
  string temp = "## [configuration file] STRING: a configuration file defining input and output files, e.g. config.txt \n";
  temp += "## ${X}: the value of X is defined in [configuration file] \n";
  temp += "## [sample name] STRING: sample name, defined at ${file_pn}, prefix of *bam files\n";
  return temp;
}

string Parameter2(){
  return "## [chromosome] STRING: one item from {chr1, chr2, ..., chrX}\n";
}

string _Parameter2(){
  string temp = "## [sample index] INT, the index of the sample in ${file_pn} \n";
  temp +=  "   e.g. [sample index] = {0, 1, 2, ..., s-1}; s is the number of samples in ${file_pn} \n";
  return temp;
}

string Parameter0(){
  return "####################################################\n";
}

int UsgBuildDist(){
  cout << "############################################################################################\n"
       << "######### build insert length distribution, before runing alu_delete or alu_insert #########\n"
       << "############################################################################################\n"
       << "USAGE: opt3/build_dist [configuration file] [sample name] \n"
       << Parameter1() << endl;
  return 0;
}

int UsgAluDelete(){
  cout << "################################################################################\n"
       << "######### detect polymorphic Alu that is deleted wrt. reference genome #########\n"
       << "################################################################################\n"
       << "## there are 3 steps in total \n" 
       << Parameter0()
       << "## step 0: create folders for temporary files\n" 
       << "USAGE: opt3/alu_delete [configuration file] 0\n"
       << Parameter0()
       << "## step 1: detect each alu delettions for each individual \n" 
       << "USAGE: opt3/alu_delete [configuration file] 1 [sample name] \n"
       << Parameter1() 
       << Parameter0()
       << "## step 2: create a big vcf file for genotype calls of all samples\n" 
       << "USAGE: opt3/alu_delete [configuration file] 2 \n"
       << Parameter0()
       << Parameter1() << endl; 
  return 0;
}

int UsgAluInsert(){
  cout << "##################################################################################\n"
       << "######### detect polymorphic Alu that is inserted wrt. reference genome ##########\n"
       << "##################################################################################\n"
       << "## there are 9 steps in total \n" 
       << Parameter0()
       << "## step 0: create folders for temporary files\n" 
       << "USAGE: opt3/alu_insert [configuration file] 0\n"
       << Parameter0()
       << "## step 1: collect alu mate for each sample\n" 
       << "USAGE: opt3/alu_insert [configuration file] 1 [sample name] \n"
       << Parameter0()
       << "## step 1+: make file ${file_pn_used} to include all sample names that succeed in previous step\n"
       << "USAGE: cp ${file_pn} ${file_pn_used}, if all samples succeeds\n"   
       << Parameter0()
       << "## step 2: combine insertion regions from multiple individuals, write to ${output_path}/insert_alu1/insert_pos.*\n"
       << "USAGE: opt3/alu_insert [configuration file] 2 \n"  // combine_pos_pns
       << Parameter0()
       << "## step 3: write soft-clipped reads to ${output_path}/insert_alu1/clip/chr*/[sample name] \n"  // clipReads_pn
       << "USAGE: opt3/alu_insert [configuration file] 3 [sample name] \n"  
       << Parameter0()
       << "## step 4: write to insert_alu1/clip/chr*_pos/*, write refined insertion regions at @file_alu_insert1/clip/chr*.clip_region \n" //clipReads_pns
       << "USAGE: opt3/alu_insert [configuration file] 4 [chromosome] \n" 
       << Parameter0()
       << "## step 5: write alu and clip reads to ${output_path}/insert_alu1/cons/chr*/[sample name] \n" // write_tmp0_pn
       << "USAGE: opt3/alu_insert [configuration file] 5 [sample name] \n" 
       << Parameter0()
       << "## step 6: write exact insertion breakpoints at ${output_path}/insert_alu1/cons/chr*_clip_pass\n" // clipPos_pns
       << "USAGE: opt3/alu_insert [configuration file] 6 [chromosome] \n" 
       << Parameter0()
       << "## step 7: genotype calling for each sample, write details at ${output_path}/insert_alu2/*\n" //  write_tmp2_pn 
       << "USAGE: opt3/alu_insert [configuration file] 7 [sample name] \n" 
       << Parameter0()
       << "## step 8: create a big vcf file for genotype calls of all samples \n" //  fixed_vcf_pns
       << "USAGE: opt3/alu_insert [configuration file] 8 \n" 
       << Parameter0()
       << Parameter1() 
       << Parameter2() << endl;
  return 0;
}

int main( int argc, char* argv[] ) {

  if (argc == 1 ) {
    return PrintUsg();
  }
  
  string func = argv[1];
  
  if (func.empty())
    return PrintUsg();

  cout << endl;
  if ( func == "build_dist" ) {
    return UsgBuildDist();
  } else if ( func == "alu_delete" ) {
    return UsgAluDelete();
  } else if ( func == "alu_insert" ) {
    return UsgAluInsert();
  } else {
    cerr << "ERROR: invalid functions!\n";
    return PrintUsg();
  }
 
}

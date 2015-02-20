#define SEQAN_HAS_ZLIB 1
#include "insert_utils.h"
#include "delete_utils.h"
#include "gtest/gtest.h"

int GetGp(string & config_file, string  line, int errCode = -1) {
  ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
  
  string file_pn_used = cf_fh.get_conf( "file_pn_used");
  std::vector <string> pns_used;
  GetUsedSampleName( file_pn_used, pns_used);
  string pn = *pns_used.begin();  // use the first pn

  GetIOPath *io_path = new GetIOPath(config_file);
  string file_dist_prefix = io_path->insertlen_dist_;

  map <int, EmpiricalPdf *> pdf_rg;    
  map<string, int> isPairEnd;
  read_pdf_pn(file_dist_prefix, pn, "100_1000_5", pdf_rg, isPairEnd);
  float *gp = new float[3];
  parseline_ins(line, cout, pdf_rg, 3.5, 300, errCode, true, gp); 
  if ( gp[0] >= max(gp[1], gp[2]) ) return 2;
  else if ( gp[1] >= max(gp[0], gp[2]) ) return 1;
  
  delete io_path;
  return 0;
}

TEST(GenoCallTest, Input1) {
  string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
  EXPECT_EQ(0, GetGp(config_file, "chr21 40876480 40876480,0 1 0 15 0"));
}

TEST(GenoCallTest, Input12) {
  string config_file = "/home/qianyuxx/faststorage/Alu/config.dk";
  EXPECT_EQ(1, GetGp(config_file, "chr21 40876480 40876480,0 1 0 10 0"));
}



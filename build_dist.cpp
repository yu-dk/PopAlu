// build distribution file based on flow cells (reading group)
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include "utils.h"

ParametersAlu p_alu = ParametersAlu();
string rg_default = "default";

void write_counts(map <seqan::CharString, map<int, int> > &rg_lenCounts, string &rg_lenCnt_file){
  map <seqan::CharString, map<int, int> >::iterator sii;
  map<int, int>::iterator si;
  cout << "output to: " << rg_lenCnt_file << "*\n"
       << "in total " << rg_lenCounts.size() << " RG groups\n";
  for (sii = rg_lenCounts.begin(); sii != rg_lenCounts.end() ; sii++) {
    cout << "RG: " << sii->first << endl;
    ofstream fout((rg_lenCnt_file + toCString(sii->first)).c_str());
    for (si = (sii->second).begin(); si != (sii->second).end(); si++) 
      fout << si->first << " "  << si->second << endl;
    fout.close();
  }
}

void read_pn(map <seqan::CharString, map<int, int> > &rg_lenCnt, string &bamInput_file, int jump_first, int max_read_num){
  seqan::BamAlignmentRecord record;  
  seqan::BamStream bamStreamIn(bamInput_file.c_str());
  assert(isGood(bamStreamIn));
  unsigned idx_RG;
  int i = 0;
  while (!atEnd(bamStreamIn)) {
    readRecord(record, bamStreamIn);
    if ( record.rID != record.rNextId ) continue;
    if ( record.rID > 20 ) break;  // only count from chr0 - chr21
    i++;    
    if ( i < jump_first) continue;
    if ( max_read_num and i > max_read_num ) break;    
    
    //if ( i % 1000 == 0) cout << "done " << i << " " << rg_lenCnt.size() << endl;
    if ((not hasFlagQCNoPass(record) ) and hasFlagAllProper(record) and (not hasFlagDuplicate(record)) and hasFlagMultiple(record) ) {
      if (hasFlagFirst(record)) continue; // look at only one end of the pair      
      seqan::BamTagsDict tags(record.tags);
      seqan::CharString rg;
      if (findTagKey(idx_RG, tags, "RG"))
	rg = getTagValue(tags, idx_RG); // idx_RG is the idx at the tag dict.
      else
	rg = ::rg_default;
      addKey(rg_lenCnt[rg], abs(record.tLen));    
    }    
  }  
  seqan::close(bamStreamIn);
}

void write_pdf(string &f_count, string &f_prob, int _len_min, int bin_width){
  int insertlen, count;
  map <int, int> bin_counts;
  fstream fin;
  fin.open(f_count.c_str());
  if (!fin) {
    cerr << "ERROR: " << f_count << " not exist\n";
    exit(0);
  }
  int mLen=0, mCnt=0;
  // get len_min and len_max
  while (fin >> insertlen >> count) {
    if (count > mCnt) {
      mLen = insertlen;
      mCnt = count;
    }
  }
  fin.close();
  int len_min = max(mLen - 500, _len_min);
  int len_max = mLen + 500;
  int counts = 0;
  // reads of insert length in [len_min, len_max]
  fin.open(f_count.c_str());
  while (fin >> insertlen >> count) {
    if (insertlen < len_min) continue;
    if (insertlen >= len_max) break;    
    int id_bin = (insertlen - len_min)/bin_width;
    addKey(bin_counts, id_bin, count);
    counts += count;
  }
  fin.close();
  int half_bin_width = bin_width / 2;
  ofstream fout(f_prob.c_str());
  for (map<int, int>::iterator bc = bin_counts.begin(); bc != bin_counts.end(); bc++) 
    fout << len_min + (bc->first) * bin_width + half_bin_width << " " << (bc->second) / (float)counts << endl;
  fout << "0 " << counts << " " << mLen << endl; // last line for extra info
  fout.close();
  cout << "written into " << f_prob << endl;
}


int main( int argc, char* argv[] )
{
  if (argc == 3) {
    string config_file = argv[1];
    string pn = argv[2];

    ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
    GetIOPath *io_path = new GetIOPath(config_file);

    string path0 = io_path->insertlen_;
    string path_count = io_path->insertlen_count_;
    string path_prob =  io_path->insertlen_dist_;
    string bamInput_file = io_path->GetBamFileName(pn);
    
    string rg_lenCnt_file = path_count + pn + ".count."; /* + RG*/
    map <seqan::CharString, map<int, int> > rg_lenCnt; // strata by reading group 
    read_pn(rg_lenCnt, bamInput_file, 5e3, 2e8);  // use only 20% of 30X reads to estimate 
    write_counts(rg_lenCnt, rg_lenCnt_file); 

    stringstream ss;
    ss.str(::p_alu.pdf_param);
    int len_min, len_max, bin_width; ///int len_min = 100, len_max = 1000, bin_width = 5; 
    string token;
    getline(ss, token, '_');
    seqan::lexicalCast2(len_min, token);
    getline(ss, token, '_');
    seqan::lexicalCast2(len_max, token);
    getline(ss, token, '_');
    seqan::lexicalCast2(bin_width, token);    
    
    ofstream fout( (path_prob + "RG." + pn).c_str()); 
    for (map <seqan::CharString, map<int, int> >::iterator ri = rg_lenCnt.begin(); ri != rg_lenCnt.end(); ri++ ) {
      string rg = seqan::toCString(ri->first); // need to use toCString, otherwise has ^@ in the ending (due to binary file)
      fout << rg << endl;
      string f_count = path_count + pn + ".count." + rg;
      string f_prob = path_prob + pn + ".count." + rg + "." + ::p_alu.pdf_param;
      write_pdf(f_count, f_prob, len_min, bin_width);    
    }
    fout.close(); 
    delete io_path;
    
  } else if ( argc == 2 ) {  
    
    string opt =  argv[1];
    if (opt != "mason_db") // write database for mason simulation
      return 0;
    
    string path_input = "/home/qianyuxx/faststorage/AluDK/inputs/alu_hg19/";
    string path_output = "/home/qianyuxx/faststorage/AluDK/inputs/alu_hg19_mason/";
    AluRefPos *alurefpos;

    vector<string> chrns;
    for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
    chrns.push_back("chrX");
    chrns.push_back("chrY");
    for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++) {
      string chrn = *ci;
      alurefpos = new AluRefPos(path_input + "alu_" + chrn, 200, 600);  // no Alu within 600 bp 
      ofstream fout ( (path_output + "alu_" + chrn).c_str() );
      while ( alurefpos->nextdb() ) {
	string at;
	int beginPos, endPos;
	beginPos = alurefpos->get_beginP();
	endPos = alurefpos->get_endP();
	at = alurefpos->get_type();
	fout << chrn << " " << beginPos << " " << endPos << " " <<  at << endl;
      }
      delete alurefpos;
      fout.close();
    }
    
  }  
  return 0;
}


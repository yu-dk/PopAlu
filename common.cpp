#include "common.h"

ConfigFileHandler::ConfigFileHandler(string config_file) {
  if ( config_file.size() < 2 ) {
    cerr << "config file is empty!\n";
    exit(0);
  }
  ifstream fin(config_file.c_str());
  if (!fin) {
    cerr << "config file: " << config_file << " does not exist!\n";
    exit(0);
  }
  stringstream ss;
  string line, key, value;
  while (getline(fin, line)) {
    if (line[0] == '#')  continue;
    ss.clear(); ss.str( line );  
    ss >> key >> value;
    configs[key] = value;
  }
  fin.close();
}

string ConfigFileHandler::get_conf(string key, bool required) {
  map <string, string>::iterator ci;
  if ( (ci = configs.find(key)) != configs.end() )
    return ci->second;
  if (!required) return "";
  cerr << "#### ERROR #### key: " << key <<" doesn't exist\n";
  exit(1);
}

string int_to_string(int i) {
  stringstream ss;
  ss << i;
  return ss.str();
}

string get_pn(string pn_file, int idx_pn){
  ifstream fin( pn_file.c_str());
  string pn;
  bool pn_exist = false;
  int i = 0;
  while (fin >> pn) {
    if (i == idx_pn ) {
      fin.close();  
      pn_exist = true;
      break;
    }
    i++;
  }
  assert(pn_exist);
  return pn;
}

void get_pn(string pn_file, map<int, string> &ID_pn){
  ID_pn.clear();
  ifstream fin( pn_file.c_str());
  assert(fin);
  string pn;
  int i = 0;
  while (fin >> pn)  ID_pn[i++] = pn;
  fin.close();
}

bool GetUsedSampleName(string fn, vector <string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  if (!fin) 
    return false;
  string pn;
  while (fin >> pn) pns_used.push_back(pn);
  fin.close();
  return pns_used.size() > 1;
}

bool GetUsedSampleName(string fn, std::set <string> & pns_used) {
  ifstream fin(fn.c_str()) ;
  if (!fin) 
    return false;
  string pn;
  while (fin >> pn) pns_used.insert(pn);
  fin.close();
  return pns_used.size() > 1;
}

int check_file_size(string fn){
  FILE * pFile = fopen(fn.c_str(), "r");
  if (pFile == NULL) return -1; // not exist
  else {
    fseek (pFile, 0, SEEK_END);   // non-portable
    long size = ftell(pFile);
    fclose (pFile);
    return size;  // empty OR >0
  }
}

void parse_chrns(string s_chrns, vector<string> &chrns){
  stringstream ss;
  ss.str(s_chrns);
  string m;
  while (getline(ss, m, ','))
    chrns.push_back(m);
}

void split_by_sep(string &str, string &m, string &n, char sep ){
  stringstream ss;
  ss.str(str);
  getline(ss, m, sep);
  getline(ss, n, ' ');
}

int major_key_freq (vector <int> & ps, int & k1, int bin_width, float freq_th, int minInput) {
  k1 = 0;
  if ( (int) ps.size() <= minInput ) return 0;
  map <int, int> pos_cnt;  
  for (vector <int>::iterator pi = ps.begin(); pi != ps.end(); pi ++ )  
    addKey(pos_cnt, round_by_resolution(*pi, bin_width), 1);
  multimap<int, int> cnt_pos = flip_map(pos_cnt);  
  multimap<int, int>::reverse_iterator it = cnt_pos.rbegin();
  pos_cnt.clear();  
  int _k1, _k2;
  _k1 = it->second; 

  if (cnt_pos.size() > 1 ) {
    it ++;
    _k2 =  it->second;
    if ( abs(_k1 - _k2 ) == bin_width ) {
      _k1 = (_k1 + _k2) / 2;
    }
  }   

  // find most common position;
  int match_cnt = 0;
  for (vector <int>::iterator pi = ps.begin(); pi != ps.end(); pi ++ ) 
    if ( abs( *pi - _k1) <= (int) bin_width * 1.1) {
      addKey(pos_cnt, *pi , 1);
      match_cnt++;
    }
  
  assert( !pos_cnt.empty() );
  cnt_pos = flip_map(pos_cnt);
  it = cnt_pos.rbegin();
  if ( match_cnt / (float) ps.size() >= freq_th ) 
    k1 = it->second;
  return it->first;  // return matched count
}

int major_two_keys (vector <int> & ps, int & k1, int & k2, int & kf1, int & kf2, int bin_width, float freq_th, bool debugprint ) {
  k1 = k2 = 0;
  kf1 = kf2 = 0;
  if (ps.empty() ) 
    return 0;
  if (ps.size() == 1) {
    k1 = *(ps.begin());
    return 0;
  }
  map <int, int> pos_cnt, pos_cnt2;  
  float f1, f2;
  for (vector <int>::iterator pi = ps.begin(); pi != ps.end(); pi ++ ) 
    addKey(pos_cnt, round_by_resolution(*pi, bin_width), 1);
  multimap<int, int> cnt_pos;
  multimap<int, int>::reverse_iterator it;
  cnt_pos = flip_map(pos_cnt);
  it = cnt_pos.rbegin();
  pos_cnt.clear();

  int _k1 = 0, _k2 = 0;
  _k1 =  it->second;
  f1 = it->first / float( ps.size() );
  if ( cnt_pos.size() > 1 ) {
    it++;
    _k2 =  it->second;    
    f2 = it->first / float( ps.size() );
    if ( abs(_k1 - _k2 ) == bin_width ) {
      _k1 = (_k1 + _k2) / 2;
      f1 += f2;
      f2 = 0;
    } 
    if ( f2 < freq_th )
      _k2 = 0;
  }

  if (f1 < freq_th)  return 0;

  if ( debugprint ) 
    cout << "round " << _k1 << " " << _k2 << endl;
  
  for (vector <int>::iterator pi = ps.begin(); pi != ps.end(); pi ++ ) {
    if ( abs( *pi - _k1) <= (int) bin_width * 1.3) {
      addKey( pos_cnt, *pi, 1);      
    } else if ( _k2 and abs( *pi - _k2) <= (int) bin_width * 1.3 ) {
      addKey( pos_cnt2, *pi, 1);      
    }
  }
  
  if ( !pos_cnt.empty())  { 
    cnt_pos = flip_map(pos_cnt);
    k1 = (cnt_pos.rbegin())->second;
    kf1 = (cnt_pos.rbegin())->first;
  } else k1 = 0;

  if ( _k2 and !pos_cnt2.empty() ) {
    cnt_pos = flip_map(pos_cnt2);
    k2 = (cnt_pos.rbegin())->second;
    kf2 = (cnt_pos.rbegin())->first;
  } else k2 = 0;
  
  return 0;
}

void intersect_fast0(std::set <pair<int, int> > & query, list <pair<int, int> > & db, std::set <pair<int, int> > & query_no_overlap){
  db.sort();
  query_no_overlap.clear();
  list <pair<int, int> >::iterator d2 = db.begin();
  for ( std::set <pair<int, int> >::iterator d1 = query.begin(); d1 != query.end(); d1++ ) {
    if ( (*d1).second <= 0 )  continue;
    while ( (*d1).first >= (*d2).second and d2 != db.end() )  d2++;
    if ( d2 != db.end() and  min ((*d2).second, (*d1).second ) - max ( (*d2).first, (*d1).first ) < 0 ) 
      query_no_overlap.insert( *d1 ) ;
    else if ( d2 == db.end() ) 
      query_no_overlap.insert( *d1 ) ;
  }
}

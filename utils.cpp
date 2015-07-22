#include "utils.h"

string GetIOPath::CheckFolder(string path)
{
    if ( *path.rbegin() != '/') path += '/';
    check_folder_exists(path);
    cerr << "log: creating " << path << endl;
    return path;
}


GetIOPath::GetIOPath(string config_file)
{
    cf_fh = new ConfigFileHandler(config_file);
    bam_path_ = CheckFolder (cf_fh->get_conf("file_bam_prefix") );
    output_path_ = CheckFolder (cf_fh->get_conf("output_path") );
    insertlen_ = CheckFolder (output_path_ + "insert_len/");
    insertlen_dist_ = CheckFolder (output_path_ + "insert_len/prob/");
    insertlen_count_ = CheckFolder (output_path_ + "insert_len/count/");
}

EmpiricalPdf::EmpiricalPdf(string pdf_file, string rg_name, map <string, int> & isPairEnd)
{ 
    int len;
    int cnt = 0, freq_len;
    float lenp, freq_p;
    bin_width = 0;
    ifstream fin( pdf_file.c_str());
    if (!fin)
    {
	cerr << "ERROR: " << pdf_file << " not exists!\n";
	exit(0);
    }
    fin >> len >> lenp;
    freq_len = len;
    freq_p = lenp;
    prob_vec[len] = lenp;
    this->min_len = len;
    this->max_len = 100;
    this->min_prob = lenp;
    while (fin >> len >> lenp)
    {
	if (len == 0 )
	{
	    cnt = (int) lenp;
	    break;
	}
	prob_vec[len] = lenp;
	if (this->min_prob > lenp)
	{
	    this->min_prob = lenp;
	}
	if (freq_p < lenp)
	{
	    freq_p = lenp;
	    freq_len = len;
	}
	this->max_len = len;
	if (!bin_width) bin_width = len - min_len;
    }
    fin.close();

    this->RG_mate_pair = false;
    if ( !isPairEnd.empty() )
    {
	if ( isPairEnd.find(rg_name) == isPairEnd.end() )
	{
	    if ( freq_len > 800 ) this->RG_mate_pair = true;
	    if (this->RG_mate_pair) cerr << "unknown reading group (consider as mate pair): "<< rg_name;
	    else cerr << "unknown reading group (consider as pair end): "<< rg_name;    
	}
	else if ( !isPairEnd[rg_name] )
	{
	    this->RG_mate_pair = true;
	}
    }
    /*
       if ( this->RG_mate_pair ) cout << "mate pair";
       else cout << "pair end" ;
       cout << ", max insert length " << this->max_len << " " << pdf_file << endl;   
       */
}

float EmpiricalPdf::pdf_obs(int insertlen)
{
    if (insertlen >= max_len || insertlen <= min_len) return min_prob;
    int nearby_pos = min_len + (insertlen-min_len)/bin_width*bin_width;
    map <int, float >::iterator it = prob_vec.find(nearby_pos);
    if ( it != prob_vec.end())  return it->second;
    return min_prob;
}


void EmpiricalPdf::ratio_obs(int y, int z, float log10_ratio_ub, float & py, float & pz)
{
    assert( y > z);
    float min_ratio = pow(10, - abs(log10_ratio_ub));
    if ( z < 0 )
    {
	py = 1; pz = min_ratio;
    }
    else
    {
	float yz_ratio = pdf_obs(y) / pdf_obs(z);
	if ( yz_ratio >= 1 )
	{
	    py = 1; pz = max(1/yz_ratio, min_ratio);
	}
	else
	{
	    py = max(yz_ratio, min_ratio); pz = 1;
	}
    }
}


void EmpiricalPdf::delete_map(map <int, EmpiricalPdf *> & epdf_rg)
{
    for (map <int, EmpiricalPdf *>::iterator ri = epdf_rg.begin(); ri != epdf_rg.end(); ri++) 
	delete ri->second;
    epdf_rg.clear();
}

    BamFileHandler::BamFileHandler(vector<string> &chrns, string bam_input, string bai_input, string bam_output) 
    : fn_bai(bai_input)
    , fn_bam_out(bam_output)
    , nameStoreCache(nameStore)
      , context(nameStore, nameStoreCache)
{
    if (fn_bai != "" and read(baiIndex, fn_bai.c_str()) != 0)
    {
	cerr << "ERROR: Could not read BAI index file " << fn_bai << endl;
	exit(1);
    }
    if (!open(inStream, bam_input.c_str(), "r"))
    {
	std::cerr << "ERROR: Could not open " << bam_input << " for reading.\n";
	exit(1);
    }
    seqan::BamHeader header;
    seqan::BamAlignmentRecord record;
    assert(!readRecord(header, context, inStream, seqan::Bam()) );
    if (!fn_bam_out.empty())
    {
	assert(open(outStream, bam_output.c_str(), "w") );
	assert(write2(outStream, header, context, seqan::Bam()) == 0);
    }
    int rID = -1;
    for ( vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci++)
    {
	assert ( (*ci).substr(0,3) == "chr");
	if ( !getIdByName(nameStore, *ci, rID, nameStoreCache) )
	    if (!getIdByName(nameStore, (*ci).substr(3), rID, nameStoreCache))
	    {
		cerr << "ERROR: Reference sequence named "<< *ci << " not known.\n";
		exit(1);
	    }
	//cout << *ci << " " << rID  << endl;
	chrn_rID[*ci] = rID;
	rID_chrn[rID] = *ci;
    }

}


BamFileHandler::~BamFileHandler(void)
{
    seqan::close(inStream); 
    if (!fn_bam_out.empty())
	seqan::close(outStream); 
}

bool BamFileHandler::jump_to_region(int rid, int region_begin, int region_end)
{
    assert ( !fn_bai.empty()); // need bai file if want to jump 
    bool hasAlignments = false;    
    if (!jumpToRegion(inStream, hasAlignments, context, rid, region_begin, region_end, baiIndex))
	return false;
    return hasAlignments;
}

bool BamFileHandler::jump_to_region(string chrn, int region_begin, int region_end)
{
    assert ( !fn_bai.empty()); // need bai file if want to jump 
    bool hasAlignments = false;    
    if (!jumpToRegion(inStream, hasAlignments, context, chrn_rID[chrn], region_begin, region_end, baiIndex))
	return false;
    return hasAlignments;
}

bool BamFileHandler::fetch_a_read(seqan::BamAlignmentRecord & record)
{
    if (atEnd(inStream)) return false;
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    return true;
}

string BamFileHandler::fetch_a_read(string chrn, int region_begin, int region_end, seqan::BamAlignmentRecord & record)
{
    if (atEnd(inStream)) return "stop";
    assert (!readRecord(record, context, inStream, seqan::Bam())); 
    if (record.rID != chrn_rID[chrn] || record.beginPos >= region_end) return "stop";
    if (record.beginPos + (int)getAlignmentLengthInRef(record) < region_begin) return "skip";
    return "record";
}

bool BamFileHandler::get_chrn(int query_rid, string & pairChrn)
{
    map <int, string>::iterator ri;
    if ( (ri = rID_chrn.find(query_rid)) == rID_chrn.end() )
	return false;
    pairChrn = ri->second;
    return true;
}

void BamFileHandler::print_mapping_rID2chrn()
{
    cout << "mapping from rID to chrn \n";
    for (map<int, string>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++) 
	cout << rc->first << " " << rc->second << endl;
}

bool BamFileHandler::write_a_read(seqan::BamAlignmentRecord & record)
{
    if (fn_bam_out.empty())  return false;
    return write2(outStream, record, context, seqan::Bam()) == 0;
}

BamFileHandler * BamFileHandler::openBam_24chr(string bam_input, string bai_input, string bam_output)
{
    vector<string> chrns;
    for (int i = 1; i < 23; i++)  chrns.push_back("chr" + int_to_string(i) );
    chrns.push_back("chrX");
    chrns.push_back("chrY");
    return new BamFileHandler(chrns, bam_input, bai_input, bam_output);
}

FastaFileHandler::FastaFileHandler(string fn_fa)
{
    if ( read(faiIndex, fn_fa.c_str()) )
    {
	build(faiIndex, fn_fa.c_str() ); 
	write(faiIndex, (fn_fa+".fai").c_str() ); 
    }

    only_one_seq = false;
}

FastaFileHandler::FastaFileHandler(string fn_fa, string seq_name)
{
    if ( read(faiIndex, fn_fa.c_str()) )
    {
	build(faiIndex, fn_fa.c_str() ); 
	write(faiIndex, (fn_fa+".fai").c_str() ); 
    }

    idx = 0;
    assert (getIdByName(faiIndex, seq_name, idx));
    only_one_seq = true;
}

void FastaFileHandler::fetch_fasta_upper(int beginPos, int endPos, seqan::CharString &seq)
{
    assert (only_one_seq);
    assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
    seqan::toUpper(seq);
}

void FastaFileHandler::fetch_fasta_upper(int beginPos, int endPos, seqan::CharString &seq, string seq_name)
{
    assert (!only_one_seq);
    assert (getIdByName(faiIndex, seq_name, idx));
    assert (!readRegion(seq, faiIndex, idx, beginPos, endPos));
    seqan::toUpper(seq);
}

seqan::CharString FastaFileHandler::fasta_seq(string fa_input, string seq_name,int beginPos, int endPos)
{
    FastaFileHandler * fasta_fh = new FastaFileHandler(fa_input, seq_name);
    seqan::CharString fa_seq;
    fasta_fh -> fetch_fasta_upper(beginPos, endPos, fa_seq);
    delete fasta_fh;
    return fa_seq;
}


AluconsHandler::AluconsHandler(string fn_fa, string sn):
    FastaFileHandler(fn_fa)
{
    update_seq_name(sn);
    ifstream fin( (fn_fa+".fai").c_str() );
    string line, aname;
    stringstream ss;
    while (getline(fin, line))
    {
	ss.clear(); ss.str(line);
	ss >> aname;
	seq_names.push_back(aname);
    }
    fin.close();
}

void AluconsHandler::update_seq_name(string sn)
{
    seq_name = sn;
    seqs.clear();  
    unsigned idx = 0;
    seqan::CharString seq;
    assert (getIdByName(faiIndex, seq_name, idx));
    readSequence(seq, faiIndex, idx);
    seqan::toUpper(seq);    
    seq_len = length(seq);
    seqs[1] = seq;
    seqan::ModifiedString <seqan::CharString, seqan::ModReverse > myRev(seq);
    seqan::CharString rev_seq = myRev;
    seqs[4] = rev_seq;  // 1,4 are reversed
    reverseComplement(seq);
    seqs[2] = seq;
    rev_seq = myRev;
    seqs[3] = rev_seq; // 2,3 are reversed
}

seqan::CharString AluconsHandler::fetch_alucons(int key)
{
    assert ( key >=1 and key <= 4) ;
    return seqs[key];
}



void parse_reading_group(string file_rg, map<string, int> & rg_to_idx)
{
    ifstream fin (file_rg.c_str() );
    if (!fin)
    {
	cerr << "ERROR: " << file_rg << " missing\n";
	exit(0);
    }
    int idx = 0;
    string rg;
    while (fin >> rg) rg_to_idx[rg] = idx++;
    fin.close();
    //cout << rg_to_idx.size() << " reading groups exist\n";
}

int get_rgIdx(map<string, int> & rg_to_idx, seqan::BamAlignmentRecord & record)
{
    seqan::BamTagsDict tags(record.tags);
    unsigned idx_rg1;  // idx in bam file
    if (!findTagKey(idx_rg1, tags, "RG"))
	return 0;
    string rg_name = seqan::toCString(getTagValue(tags, idx_rg1));
    int idx_rg2 = 0; // idx in my file 
    map <string, int>::iterator rgi = rg_to_idx.find(rg_name);
    if (rgi != rg_to_idx.end()) idx_rg2 = rgi->second;
    return idx_rg2;
}

void parse_cigar(string cigar, list <char> & opts, list <int> & cnts)
{
    string cnt;
    int cnt_int;
    opts.clear();
    cnts.clear();
    for (size_t i = 0; i < cigar.size(); i++)
    {
	if ( !isdigit(cigar[i]) )
	{
	    opts.push_back(cigar[i]);
	    if (!cnt.empty())
	    {
		seqan::lexicalCast2(cnt_int, cnt);
		cnts.push_back(cnt_int);
	    }
	    cnt = "";
	}
	else
	{
	    cnt += cigar[i];
	}
    }
    seqan::lexicalCast2(cnt_int, cnt);
    cnts.push_back(cnt_int);  
}

string get_cigar(seqan::BamAlignmentRecord &record)
{
    stringstream ss;
    for (unsigned li = 0; li < length(record.cigar); li++) 
	ss << record.cigar[li].operation << record.cigar[li].count;
    return ss.str();
}

void debug_print_read(seqan::BamAlignmentRecord &record, ostream & os)
{
    os << record.qName << " " << record.beginPos << " " << record.beginPos + getAlignmentLengthInRef(record)  << " " ;
    os << get_cigar(record) << " " << record.pNext << " " << record.tLen << endl;
}

int numOfBestHits(seqan::BamAlignmentRecord &record)
{
    seqan::BamTagsDict tagsDict(record.tags);
    unsigned myIdx = 0, valInt = 1;  
    if ( findTagKey(myIdx, tagsDict, "X0")) extractTagValue(valInt, tagsDict, myIdx);
    return max((int)valInt, 1);
}

bool find_read(string &bam_input, string &bai_input, string &chrn, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region)
{ // flank_region = 0, print this read; otherwise, print its pair
    int that_begin = this_pos - max(flank_region, 10); 
    int that_end = this_pos + max(flank_region, 10);
    TNameStore      nameStore;
    TNameStoreCache nameStoreCache(nameStore);
    TBamIOContext   context(nameStore, nameStoreCache);
    seqan::BamHeader header;
    seqan::BamAlignmentRecord record;  
    seqan::Stream<seqan::Bgzf> inStream;  
    open(inStream, bam_input.c_str(), "r");
    seqan::BamIndex<seqan::Bai> baiIndex;
    assert( !read(baiIndex, bai_input.c_str())) ;
    int rID;
    assert ( !readRecord(header, context, inStream, seqan::Bam()) ); 
    assert ( getIdByName(nameStore, chrn, rID, nameStoreCache) ); // change context ??
    bool hasAlignments = false;
    jumpToRegion(inStream, hasAlignments, context, rID, that_begin, that_end, baiIndex);
    if (!hasAlignments) return false;
    while (!atEnd(inStream))
    {
	assert (!readRecord(record, context, inStream, seqan::Bam())); 
	if (record.rID != rID || record.beginPos >= that_end) break;
	if (record.beginPos < that_begin) continue;
	if (record.qName == this_qName)
	{
	    if ((flank_region > 0 and record.beginPos != this_pos) or (flank_region == 0) )
	    {
		that_record = record;
		return true;
	    }
	}
    }
    return false;
}

bool get_trim_length( seqan::BamAlignmentRecord & record, int & trimb, int & trime, int bpQclip)
{  
    int beg, end;
    for( beg = 0; (beg < (int)length( record.seq )) && ( qualToInt( record.qual[beg] ) < bpQclip ); beg++ )
    {}
    for( end = length( record.seq )-1; (end > 0) && ( qualToInt( record.qual[end] ) < bpQclip ) ; end-- )
    {}
    trimb = beg;
    trime = length( record.seq) - end - 1;  
    return end - beg > 20;  // at least 20 bp
}

bool trim_clip_soft_first( seqan::BamAlignmentRecord & record, int & clipLen, seqan::CharString & clipSeq, int bpQclip)
{  
    assert(record.cigar[0].operation == 'S');
    int beg, end;
    for( beg = 0; (beg < (int)length( record.seq )) && ( qualToInt( record.qual[beg] ) < bpQclip ); beg++ )
    {}
    for( end = length( record.seq )-1; (end > 0) && ( qualToInt( record.qual[end] ) < bpQclip ) ; end-- )
    {}
    if ( end - beg < CLIP_BP) return false;
    int cut_beg = beg;
    int cut_end = length( record.seq) - end - 1;
    clipLen = record.cigar[0].count - cut_beg; // positive, right clip 
    if ( end - beg >= CLIP_BP and clipLen >= CLIP_BP and (int)( length(record.seq) - record.cigar[0].count) >= cut_end + CLIP_BP )
    {
	clipSeq = infix( record.seq, beg, end+1 );
	return true;
    }
    return false;
}

bool trim_clip_soft_last( seqan::BamAlignmentRecord & record, int & clipLen, seqan::CharString & clipSeq, int bpQclip)
{  // *M*S
    int idx = length(record.cigar) - 1;
    assert(record.cigar[idx].operation == 'S');
    int beg, end;
    for( beg = 0; (beg < (int)length( record.seq )) && ( qualToInt( record.qual[beg] ) < bpQclip ); beg++ )
    {}
    for( end = length( record.seq )-1; (end > 0) && ( qualToInt( record.qual[end] ) < bpQclip ) ; end-- )
    {}
    if ( end - beg < CLIP_BP) return false;
    int cut_beg = beg;
    int cut_end = length( record.seq) - end - 1;
    clipLen = record.cigar[idx].count - cut_end; 
    if ( end - beg >= CLIP_BP and clipLen >= CLIP_BP and (int)(length(record.seq) - record.cigar[idx].count) >= cut_beg + CLIP_BP)
    {
	clipSeq = infix( record.seq, beg, end+1 );
	clipLen = - clipLen;  // negative, left clip
	//cout << "err " << cut_beg << " "  << cut_end << " " <<  length(record) << " " << record.cigar[idx].count << " " << cut_beg << endl;
	return true;
    }
    return false;
}

RepMaskPos::RepMaskPos(string file_rmsk, vector<string> &chrns, int join_len)
{
    for (vector<string>::iterator ci = chrns.begin(); ci != chrns.end(); ci ++ )
    {
	string file_alu = file_rmsk + "alu_" + *ci;
	ifstream fin( file_alu.c_str() );
	if (!fin)
	{
	    cerr << "file " << file_alu << " not exist!\n";
	    exit(1);
	}
	string line, chrn, repType;
	int beginPos, endPos, beginPos_pre, endPos_pre;
	stringstream ss;  
	map<string, std::set <RepDB1> > rep_db;
	while (getline(fin, line))
	{
	    ss.clear(); ss.str( line );
	    ss >> chrn >> beginPos >> endPos >> repType;
	    assert(repType.substr(0,3) == "Alu");
	    assert(chrn == *ci );
	    RepDB1 one_rep = RepDB1(beginPos, endPos, "");
	    rep_db[chrn].insert(one_rep);  // sorted
	}
	fin.close();    
	// combine two blocks if they are nearby
	std::set <RepDB1> ::iterator pi = rep_db[chrn].begin();
	beginPos_pre = (*pi).begin_pos;
	endPos_pre = (*pi).end_pos;
	pi ++;
	for ( ; pi != rep_db[chrn].end(); pi++)
	{
	    beginPos = (*pi).begin_pos;
	    endPos = (*pi).end_pos;
	    if (beginPos - endPos_pre > join_len)
	    {  // close this block, create new
		beginP[chrn].push_back(beginPos_pre);
		endP[chrn].push_back(endPos_pre);    
		beginPos_pre = beginPos;
	    }
	    endPos_pre = endPos;
	}
	beginP[chrn].push_back(beginPos_pre);
	endP[chrn].push_back(endPos_pre);    
	assert(beginP[chrn].size() == endP[chrn].size() );
    }

}

RepMaskPos::~RepMaskPos(void)
{
    for ( map<string, vector<int> >::iterator pi = beginP.begin(); pi != beginP.end(); pi++ ) 
	(pi->second).clear();
    beginP.clear();
    for ( map<string, vector<int> >::iterator pi = endP.begin(); pi != endP.end(); pi++ ) 
	(pi->second).clear();
    endP.clear();
}

AluRefPosRead::AluRefPosRead(string file_alupos, int minLen)
{
    ifstream fin( file_alupos.c_str());
    if (!fin) 
	try
	{
	    throw 1;    
	}
    catch(int e)
    {
	cerr << "#### ERROR #### file: "<< file_alupos << " not exists!!!" << endl;
    }

    string _tmp1, _tmp2, alu_type;
    char chain;
    int bp, ep;
    while (fin >> _tmp1 >> bp >> ep >> alu_type >> _tmp2 >> chain)
    {
	if (ep - bp < minLen) continue;
	if (beginP.empty()) minP = bp;
	beginP.push(bp);
	endP.push(ep);
	strandP.push(chain);
	aluType.push(alu_type);
    }
    maxP = ep;
    //cerr << "queue from " << file_alupos << " with " << beginP.size() << " loci, " << minP << " to " << maxP << endl;
    fin.close();
}

int AluRefPosRead::updatePos(int &beginPos, int &endPos, char &chain, string & alu_type)
{
    if (!beginP.size()) return 0;
    beginPos = beginP.front();
    endPos = endP.front();
    chain = strandP.front();
    alu_type = aluType.front();
    beginP.pop();
    endP.pop();
    strandP.pop();
    aluType.pop();
    return 1;
}


AluRefPos::AluRefPos(string fn, int minLen_alu)
{
    ifstream fin( fn.c_str());
    assert(fin);
    string chrn, at, line;
    int bp, ep;
    stringstream ss;
    while ( getline(fin, line))
    {
	ss.clear(); ss.str(line);
	ss >> chrn >> bp >> ep >> at;
	if ( minLen_alu > 0 and ep - bp < minLen_alu) continue;
	RepDB1 one_alu = RepDB1(bp, ep, at);
	alu_db.insert( one_alu );
    }
    fin.close();  
    adi = alu_db.begin(); // initialize pointer
    db_size = alu_db.size();
}

AluRefPos::AluRefPos(string fn, int minLen_alu, int minDist_neighbor)
{
    // assume input is sorted 
    ifstream fin( fn.c_str());
    assert(fin);
    string chrn, at, line;
    int bp, ep, bp_pre, ep_pre;
    stringstream ss;
    list <RepDB1> alu_db_list;
    getline(fin, line);
    ss.str(line);
    ss >> chrn >> bp_pre >> ep_pre;
    while ( getline(fin, line))
    {
	ss.clear(); ss.str(line);
	ss >> chrn >> bp >> ep >> at;
	bool flag = true;
	if ( !alu_db_list.empty() and abs( (alu_db_list.back()).end_pos - bp ) < minDist_neighbor )
	{
	    alu_db_list.pop_back();
	    flag = false;
	}
	else if ( minLen_alu > 0 and ep - bp < minLen_alu)
	{      
	    flag = false;
	}
	else if ( abs(ep_pre - bp) < minDist_neighbor )
	{
	    flag = false;
	}

	///cout << flag << " " << bp << " " << alu_db_list.size() << endl;
	if (flag )
	{
	    RepDB1 one_alu = RepDB1(bp, ep, at);
	    alu_db_list.push_back( one_alu );
	}
	bp_pre = bp;
	ep_pre = ep;
    }
    fin.close();  
    RepDB1 one_alu = RepDB1(bp_pre, ep_pre, at);
    alu_db_list.push_back( one_alu );

    for (list <RepDB1>::iterator al = alu_db_list.begin(); al != alu_db_list.end(); al++ ) 
	alu_db.insert(*al);
    alu_db_list.clear();

    adi = alu_db.begin(); // initialize pointer
    db_size = alu_db.size();
}

bool AluRefPos::within_alu( int pos)
{
    for (std::set <RepDB1>::iterator di = alu_db.begin(); di != alu_db.end(); di++ ) 
    {
	if ( pos >= (*di).begin_pos and pos <= (*di).end_pos )
	    return true;
	if ( pos > (*di).end_pos )
	    break;
    }
    return false;
}

void AluRefPos::debug_print()
{
    cout << "DB size " << db_size << endl;
    for (std::set <RepDB1>::iterator di = alu_db.begin(); di != alu_db.end(); di++ ) 
	cout << (*di).begin_pos << " " << (*di).end_pos << endl;
}

void get_RG_info(string fn, string pn, map<string, int > & isPairEnd)
{
    ifstream fin(fn.c_str());
    string line, _pn, _rg;
    int _ispairend;
    stringstream ss;
    while (getline(fin, line))
    {
	ss.clear(); ss.str(line);
	ss >> _pn;
	if (_pn != pn) continue;
	ss >> _rg >> _ispairend;
	isPairEnd[_rg] = _ispairend;
    }
    fin.close();
}

void read_pdf_pn( string prefix, string pn, string pdf_param,  map <int, EmpiricalPdf *> & pdf_rg, map<string, int> &isPairEnd)
{
    ifstream fin( get_name_rg(prefix, pn).c_str() );
    int idx = 0; 
    string rg;
    while (fin >> rg) 
	pdf_rg[idx++] = new EmpiricalPdf( get_name_rg_pdf(prefix, pn, rg, pdf_param), rg, isPairEnd);
    fin.close();
}


seqan::Pair<seqan::CigarElement<>::TCount> mappedInterval(seqan::String<seqan::CigarElement<> > & cigar)
{
    typedef seqan::CigarElement<>::TCount TSize;
    TSize len = 0;
    TSize beginPos = 0;
    TSize endPos = 0;
    seqan::Iterator<seqan::String<seqan::CigarElement<> > >::Type itEnd = end(cigar);
    for (seqan::Iterator<seqan::String<seqan::CigarElement<> > >::Type it = begin(cigar); it < itEnd; ++it)
    {
	len += (*it).count;
	switch ((*it).operation)
	{
	    case 'S': case 'H':
		if (it == begin(cigar))
		    beginPos += (*it).count;
		break;
	    case 'D':
		len -= (*it).count;
		break;
	    case 'M': case 'I':
		endPos = len;
		break;
	}
    }
    return seqan::Pair<TSize>(beginPos, endPos);
}

double avgQuality(seqan::CharString & qual, seqan::Pair<seqan::CigarElement<>::TCount> & interval)
{
    if (interval.i1 >= interval.i2) return 0;
    if (length(qual) == 0) return 50; // Accept undefined quality strings ('*' in sam format).    
    seqan::Iterator<seqan::CharString>::Type it = begin(qual);
    seqan::Iterator<seqan::CharString>::Type itEnd = begin(qual);
    it += interval.i1;
    itEnd += interval.i2;    
    unsigned totalQual = 0;
    while (it != itEnd)
    {
	totalQual += *it;
	++it;
    }
    return totalQual/(interval.i2-interval.i1);
}

bool QC_insert_read_qual( seqan::BamAlignmentRecord &record)
{  
    seqan::Pair<seqan::CigarElement<>::TCount> mi = mappedInterval(record.cigar);
    if (avgQuality(record.qual, mi) < 20)  // min average quality should be 20
	return false;
    return QC_insert_read(record);
};

string phred_scaled(float p0, float p1, float p2)
{
    float pmax = max(p0, max(p1, p2));
    //cout << p0 << " " << p1 << " " << p2 << " " << pmax << " " << p0/pmax << endl;
    string s0 =  ( p0 == pmax ) ?  "0" : phred_log(p0/pmax);
    string s1 =  ( p1 == pmax ) ?  "0" : phred_log(p1/pmax);
    string s2 =  ( p2 == pmax ) ?  "0" : phred_log(p2/pmax);
    return s0  + "," + s1 + "," + s2; 
}

string replace_str0_str(string fn, string chrn, string chr0)
{
    size_t found = fn.find(chr0);
    assert (found != string::npos);
    return fn.replace(found, chr0.size(), chrn);
}

void get_min_value(map <int, float> & m, float & min_val, int & min_key)
{
    map <int, float>::iterator mi = m.begin();
    min_val = mi->second;
    min_key = mi->first ;
    mi++;
    while ( mi != m.end() )
    {
	if ( min_val > mi->second)
	{
	    min_val = mi->second;
	    min_key = mi->first;
	}
	mi++;
    }
}

void log10P_to_P(float *log_gp, float *gp, int max_logp_dif)
{  
    assert( max_logp_dif > 0);
    map <int, float> idx_logp;
    int i, min_idx;
    for ( i = 0; i < 3; i++)
    {
	idx_logp[i] = log_gp[i];
	assert (log_gp[i] <= 0);
    }
    float min_logp;
    get_min_value(idx_logp, min_logp, min_idx);
    for ( i = 0; i < 3; i++)
    {
	if ( idx_logp[i] - min_logp > max_logp_dif )
	{
	    idx_logp[min_idx] = 1;   // set p[i] = 0 afterwards
	    min_logp = 1; 
	    break;
	}
    }

    if (min_logp > 0)
    { // check the two elements left
	get_min_value(idx_logp, min_logp, min_idx);
	for ( i = 0; i < 3; i++)
	{ 
	    if ( idx_logp[i] < 1 and idx_logp[i] - min_logp > max_logp_dif )
	    {
		idx_logp[min_idx] = 1;
		break;
	    }
	}
    }
    get_min_value(idx_logp, min_logp, min_idx);
    //cout << "log " << min_logp << " " << max_logp_dif << " " << endl;
    float ratio_sum = 0;
    for ( i = 0; i < 3; i++) 
	if (idx_logp[i] < 1) ratio_sum += pow(10, (idx_logp[i] - min_logp));
    for ( i = 0; i < 3; i++)
    {
	if (idx_logp[i] > 0) gp[i] = 0;
	else gp[i] = pow(10, (idx_logp[i] - min_logp)) / ratio_sum;
    }

}

void get_chrn(string fn, map<int, string> & rid_chrn)
{
    ifstream fin(fn.c_str());
    int rid;
    string chrn;
    while ( fin >> rid >> chrn ) 
	rid_chrn[rid] = chrn;
    fin.close();
}


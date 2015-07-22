// utils function for alu_insert 
#include "insert_utils.h"

GetIOPathInsert::GetIOPathInsert(string config_file):
    GetIOPath(config_file)
{
    path0_ = CheckFolder(output_path_ + "insert_alu0/");
    path1_ = CheckFolder(output_path_ + "insert_alu1/");
    pathClip_ = CheckFolder(path1_ + "clip/");
    pathCons_ = CheckFolder(path1_ + "cons/");
    pathDel0_ = CheckFolder(output_path_ + "insert_alu2/");
    check_folder_exists(pathDel0_ + "tmp1s/");
    check_folder_exists(pathDel0_ + "tmp2s/");
    bam_rid_chrn_ = output_path_ + "bam_rid_chrn";
}

bool AlumateINFO::sort_pos2(const AlumateINFO* a, const AlumateINFO* b)
{
    return a-> pos2 < b-> pos2;
}

void AlumateINFO::delete_list(list <AlumateINFO *> & alumate_list)
{
    for ( list <AlumateINFO *> ::iterator ai = alumate_list.begin(); ai != alumate_list.end(); ai++ ) 
	delete *ai;
    alumate_list.clear();
}

/** move read *S*M to left until can't move */
bool clipRight_move_left(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len)
{
    int i = 1;  // i = 0 is the last match by BWA
    while ( *cigar_cnts.begin() - i >= 0 and 
	    ref_fa[clipPos - refBegin - i] == read_seq[*cigar_cnts.begin() - i] )
	i++;
    if ( *cigar_cnts.begin() - (i-1) < CLIP_BP )
	return false;
    clipPos -= (i-1);  
    align_len = length(read_seq) - *cigar_cnts.begin() + (i-1);
    return true;  
}

/** move read *M*S to right until can't move */
bool clipLeft_move_right(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len)
{
    align_len = length(read_seq) - cigar_cnts.back();
    int i = 0;
    while (align_len + i < (int) length(read_seq) and ref_fa[clipPos - refBegin + i] == read_seq[align_len + i]) i++;
    if ( cigar_cnts.back() - i < CLIP_BP)
	return false;
    clipPos += i;
    align_len += i;
    return true;
}

// get inferred split positions
bool read_first2col(string fn, vector < pair<int, int> > & insert_pos, bool has_header)
{
    insert_pos.clear();
    ifstream fin(fn.c_str());
    if (!fin) return false;
    stringstream ss;
    string line, tmpv;
    if (has_header)  getline(fin, line); 
    int beginPos, endPos;
    int ni = 0;
    std::set < pair <int, int> > lefts_rights;
    while (getline(fin, line))
    {
	ni ++;
	ss.clear(); ss.str( line );  
	ss >> beginPos >> endPos;
	if ( abs (beginPos - endPos ) < MAX_POS_DIF )
	{
	    lefts_rights.insert( get_valid_pair(beginPos, endPos) ); 
	}
	else
	{
	    if (beginPos) lefts_rights.insert( make_pair(beginPos, beginPos) ); 
	    if (endPos) lefts_rights.insert( make_pair(endPos, endPos) ); 
	}
    }

    fin.close();
    if (lefts_rights.empty()) return false;
    std::set < pair <int, int> >::iterator si = lefts_rights.begin();
    int beginPre = (*si).first;
    int endPre =  (*si).second;
    si++;
    for ( ; si != lefts_rights.end(); si++)
    {
	beginPos = (*si).first;
	endPos = (*si).second;
	if ( abs(beginPre - beginPos) < MAX_POS_DIF )
	{
	    beginPre = (beginPre + beginPos) / 2 ;
	    endPre = (endPre + endPos) / 2 ;
	}
	else
	{
	    insert_pos.push_back( make_pair(beginPre, endPre) );
	    beginPre = beginPos;
	    endPre = endPos;
	}
    }
    insert_pos.push_back(make_pair(beginPre, endPre) );
    return !insert_pos.empty(); 
}

void normalize_llh(float *loggp, float th)
{
    th = - abs(th);
    float maxv = max(max(loggp[0], loggp[1]), loggp[2]);
    float minv = min(min(loggp[0], loggp[1]), loggp[2]);
    int id_max = -1, id_mid = -1, id_min = -1;
    for (int i = 0; i < 3; i ++ )
    {
	if ( abs(loggp[i] - maxv) < 0.00001 and id_max < 0) id_max = i;
	else if ( abs(loggp[i] - minv) < 0.00001 and id_min < 0) id_min = i;
	else id_mid = i;
    }
    //if ( ! (id_max > -1 and id_mid > -1 and id_min > -1) ) cout << loggp[0] << " " <<  loggp[1] << " " << loggp[2] << endl;
    assert(id_max > -1 and id_mid > -1 and id_min > -1);
    loggp[id_mid] += abs(loggp[id_max]);
    loggp[id_min] += abs(loggp[id_max]);
    loggp[id_max] = 0;
    if ( loggp[id_min] < th )
    {
	loggp[id_min] = th;
	float th_mid = loggp[id_max] - (loggp[id_max] - loggp[id_mid]) * 0.8;
	if ( loggp[id_mid] < th_mid ) loggp[id_mid] = th_mid;
    }
}

int parseline_ins(string line0, ostream & fout, map <int, EmpiricalPdf *> & pdf_rg, float logPE, int estimatedAluLen, int errCode, bool test_print, float *test_gp)
{
    //errCode1: 0/0 or .. ==> 0/1, reduce ph0 and useLen
    //errCode2: 1/1 ==> 0/1
    //errCode3: 0/1 or .. ==> 1/1, reduce ph0 ?? (0.2)
    const float ratioMax = 6.; 
    const float maxLh = 40.;

    float *log10_gp = new float[3];
    float *log10_gpu = new float[3];
    float *gp = new float[3];
    for (int i = 0; i < 3; i++)
    {
	log10_gp[i] = 0;
	log10_gpu[i] = 0;
	gp[i] = 0;
    }
    stringstream ss;
    string chrn, insertMid, debugInfo, token;
    int idx, insert_len, midCnt, midCnt_alu, midCnt_clip, clipCnt, unknowCnt;
    vector < pair <int, int> > unknowInfo;
    ss.str(line0); 
    ss >> chrn >> insertMid >> debugInfo >> midCnt_alu;
    string exact_left, exact_right;
    split_by_sep(debugInfo, exact_left, exact_right, ',');  
    //bool both_side = (exact_left != "0") and (exact_right != "0");
    if ( exact_left  == "0" ) exact_left = exact_right;
    if ( exact_right == "0" ) exact_right = exact_left;
    if ( midCnt_alu < 0 )
    {
	fout << chrn << " " << exact_left << " " << debugInfo << " " << midCnt_alu << endl;
	return 0;
    }

    ss >> midCnt_clip >> clipCnt >> unknowCnt;
    midCnt = midCnt_alu + midCnt_clip;   // sum of the two 

    int cntCnt = midCnt + clipCnt;
    int covCnt = cntCnt + unknowCnt;
    if ( covCnt < 3 or covCnt > 1000 ) return 0; // considered as missing  
    if ( cnts_confident(midCnt) and cnts_confident(clipCnt) )
    { //if ( midCnt * ratioMax < clipCnt) 
	fout << chrn << " " << exact_left << " " << debugInfo << " " << midCnt_alu << " " << midCnt_clip << " " << clipCnt << " " << unknowCnt 
	    << " 0 1 0 " << estimatedAluLen << endl;
	return 0;
    }
    float ph0 = 0.3;
    if ( midCnt <= 2 or clipCnt <= 2 ) 
	ph0 = (midCnt_clip == 0) ? 0.5 : 0.4;

    if ( errCode == 3 or errCode == 1 ) ph0 = 1. / (1. + ratioMax); 
    if (test_print) cout << "ph0 is " << ph0 << endl;

    logPE = - abs(logPE);
    float down_weight = abs(logPE) - 1.; // use less 
    float maxLhCnt = min( (float) ( abs(logPE) * 0.8 * cntCnt), maxLh);
    float maxLhLen = min( (float) ( down_weight * 0.8 * unknowCnt), maxLh);  // generally use less 
    if (cntCnt > 0 )
    {
	float logPM = log10 ( 1 - pow(10, logPE) );
	log10_gp[0] = clipCnt * logPE + midCnt * logPM;
	log10_gp[1] = midCnt * log10 (ph0) + clipCnt * log10 (1 - ph0) + (midCnt + clipCnt) * logPM  ;
	log10_gp[2] = midCnt * logPE + clipCnt * logPM;
	if (test_print)   cout << maxLhCnt << " ## " <<  log10_gp[2] << " " << log10_gp[1] << " " << log10_gp[0] << " , ";
	normalize_llh(log10_gp, maxLhCnt); 
	if (test_print)   cout << log10_gp[2] << " " << log10_gp[1] << " " << log10_gp[0] << endl;
    }
    int _unknowCnt = 0;
    if ( unknowCnt > 0)
    {
	if (errCode == 4 ) estimatedAluLen = 180; // reduce to min
	for (int i = 0; i < unknowCnt; i++)
	{
	    getline(ss, token, ':');
	    seqan::lexicalCast2(idx, token);
	    if ( pdf_rg.find(idx) == pdf_rg.end())   // fixme: check why it happens ?
		continue;
	    if (pdf_rg[idx]->RG_mate_pair)
		continue;
	    _unknowCnt++;
	    getline(ss, token, ' ');
	    seqan::lexicalCast2(insert_len, token);      
	    float p_y, p_z;
	    //if (test_print) cout << insert_len << " " << pdf_rg[idx]->pdf_obs(insert_len + estimatedAluLen) << " " << pdf_rg[idx]->pdf_obs(insert_len) << endl;   
	    pdf_rg[idx]->ratio_obs(insert_len + estimatedAluLen, insert_len, down_weight, p_y, p_z);
	    log10_gpu[0] += log10 (p_y);
	    log10_gpu[1] += log10 (ph0 * p_y + (1 - ph0) * p_z) ;
	    log10_gpu[2] += log10 (p_z);
	}
	if (test_print)   cout <<  maxLhLen << " ## " <<  log10_gpu[2] << " " << log10_gpu[1] << " " << log10_gpu[0] << " , ";
	normalize_llh(log10_gpu, maxLhLen);
	if (test_print)   cout << log10_gpu[2] << " " << log10_gpu[1] << " " << log10_gpu[0] << endl;
    }

    for (int i = 0; i < 3; i++) log10_gp[i] += log10_gpu[i];
    if (test_print)   cout << "logAll " <<  log10_gp[2] << " " << log10_gp[1] << " " << log10_gp[0] << endl;

    if (midCnt >= 3 and log10_gp[2] > max (log10_gp[1], log10_gp[0]) )
    { 
	if (errCode == -1 ) return 1;
	//else cerr << "##errCode1 " << line0 << endl;
    }
    if (clipCnt >= 3 and log10_gp[0] > max (log10_gp[1], log10_gp[2]) )
    {
	if (errCode == -1) return 2;
	//else cerr << "##errCode2 " << line0 << endl;
    }
    if ( clipCnt <= 2 and midCnt > 3 * max(clipCnt,2) and log10_gp[0] < max (log10_gp[1], log10_gp[2]) )
    {
	if (errCode == -1) return 3;
	//else cerr << "##errCode3 " << line0 << endl;
    }
    if ( log10_gp[0] > max(log10_gp[1], log10_gp[2]) and log10_gpu[2] > max(log10_gpu[0], log10_gpu[1]) )
    {
	if (errCode == -1) return 4;
	//else cerr  << "##errCode4 " << line0 << endl;
    }

    if ( !p11_is_dominant(log10_gp, - LOG10_GENO_PROB) )
    {
	log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1  
	fout << chrn << " " << exact_left << " " << debugInfo << " " << midCnt_alu << " " << midCnt_clip << " " << clipCnt << " " << _unknowCnt
	    << " " << setprecision(6) << gp[2] << " " << gp[1] << " " << gp[0] << " " << estimatedAluLen << endl;  // NB: switch 00 and 11
    }
    else 
	fout << chrn << " " << exact_left << " " << debugInfo << " " << -(midCnt + clipCnt + _unknowCnt) << endl;

    if ( test_print ) cout << gp[2] << " " << gp[1] << " " << gp[0] << endl;
    if ( test_gp != NULL )  for (int i = 0; i < 3; i++ ) test_gp[i] = gp[i];

    delete log10_gp;
    delete log10_gpu;
    delete gp;    
    return 0;
}

// used for debugging 
int parseline_cnt(string line0)
{
    stringstream ss;
    string chrn, insertMid, debugInfo, token;
    int midCnt, clipCnt, unknowCnt;
    vector < pair <int, int> > unknowInfo;
    ss.str(line0); 
    ss >> chrn >> insertMid >> debugInfo >> midCnt;
    string exact_left, exact_right;
    split_by_sep(debugInfo, exact_left, exact_right, ',');  
    if ( exact_left  == "0" ) exact_left = exact_right;
    if ( exact_right == "0" ) exact_right = exact_left;
    if ( midCnt < 0 ) return 0;
    ss >> clipCnt >> unknowCnt;
    int c1 = midCnt + clipCnt + unknowCnt;
    if ( c1 > 300 ) 
	cout << line0 << endl;
    if (midCnt == 0 or clipCnt == 0 )
	cout << c1 << " # " << midCnt << " " << clipCnt << endl;
    else 
	cout << c1 << " & " << (midCnt+0.1)/(clipCnt + 0.1) << endl;
    return 0;
}

void align_clip_to_consRef(string shortSeq, string longSeq, int & refBegin, int & refEnd,  int clipLen)
{
    const float dif_th = 0.2;
    int shortLen = shortSeq.size();
    TAlign align;
    seqan::Score<int> scoringScheme(1, -1, -2, -3); 
    resize(rows(align), 2);
    assignSource(row(align,0), shortSeq);  // 2,3
    assignSource(row(align,1), longSeq);   // 1,4, free gap at end
    globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());
    int l0 = toViewPosition(row(align, 1), 0);
    int l1 = toViewPosition(row(align, 1), longSeq.size());
    int s0 = toViewPosition(row(align, 0), 0);
    int s1 = toViewPosition(row(align, 0), shortLen);
    int align_len = min(s1, l1) - max(s0, l0);
    int align_len_shortSeq = shortSeq.size();

    if ( abs(align_len_shortSeq - align_len) <= dif_th * align_len_shortSeq 
	    and align_len >= CLIP_BP )
    {
	if ( clipLen > 0 ) // ??? 
	    refEnd = s0 + clipLen;  // might > longSeq.size() 
	else if ( clipLen < 0 ) 
	    refBegin = s1 + clipLen;  // might be negative
    }
}

bool align_alu_to_consRef(const string & shortSeq, const string & longSeq, float dif_th, string loginfo)
{
    TAlign align;
    seqan::Score<int> scoringScheme(1, -1, -2, -3); 
    resize(rows(align), 2);
    assignSource(row(align,0), shortSeq);  // 2,3
    assignSource(row(align,1), longSeq);   // 1,4, free gap at end
    globalAlignment(align, scoringScheme, seqan::AlignConfig<true, true, true, true>());
    int s0 = toViewPosition(row(align, 0), 0);
    int s1 = toViewPosition(row(align, 0), shortSeq.size());
    int l0 = toViewPosition(row(align, 1), 0); 
    int l1 = toViewPosition(row(align, 1), longSeq.size());  
    int align_len = min(s1, l1) - max(s0, l0);
    int align_len_shortSeq = shortSeq.size();
    if ( l0 > s0 )
	align_len_shortSeq -= (l0 - s0);
    else if ( s1 > l1 )
	align_len_shortSeq -= (s1 - l1);
    bool align_pass =  align_len_shortSeq >= CLIP_BP and abs(align_len_shortSeq - align_len) <= dif_th * align_len_shortSeq ;

    if ( !align_pass and !loginfo.empty()) 
    {
	cout << loginfo << ": align_len " << align_len << endl;
	cout << shortSeq << endl;
	cout << longSeq << endl;
	cout << align << endl;
    }

    return align_pass;
}


void filter_outlier_pn(string path_input, string fn_suffix, map<int, string> &ID_pn, string chrn, string file_pn_used_output, float percentage_pn_used)
{
    string line;
    int ni;
    map < string, int > pn_lineCnt;
    set <int> lineCnt;  
    for (map<int, string>::iterator pi = ID_pn.begin(); pi != ID_pn.end(); pi++)
    {
	string file_st = path_input + pi->second + "." + fn_suffix + "." + chrn;
	ifstream fin(file_st.c_str());
	ni = 0;
	while (fin >> line) ni++;
	fin.close();
	lineCnt.insert(ni);
	pn_lineCnt[pi->second] = ni;
    }

    set <int>::iterator li = lineCnt.begin();
    int cnt_th = 0;
    ni = 0;
    while ( ni++ < percentage_pn_used * lineCnt.size())
	if ( ++li != lineCnt.end() ) cnt_th = *li;
    ofstream fout(file_pn_used_output.c_str());
    for (map < string, int >::iterator pi = pn_lineCnt.begin(); pi != pn_lineCnt.end(); pi++)
	if ( pi->second <= cnt_th) fout << pi->first << endl;
    fout.close();
}

bool align_alu_cons(seqan::CharString &ref_fa, seqan::CharString alucons, float & sim_rate,float sim_th, bool read_is_clipped)
{
    TAlign align;
    int gapP = read_is_clipped ? 2 : 1;
    seqan::Score<int> scoringScheme(0, -1, - gapP, -2 * gapP);   
    resize(rows(align), 2);
    assignSource(row(align,0), ref_fa);  // 2,3
    assignSource(row(align,1), alucons);   // 1,4, free gap at end
    globalAlignment(align, scoringScheme, seqan::AlignConfig<true, false, false, true>());
    int l0 = toViewPosition(row(align, 1), 0);
    int l1 = toViewPosition(row(align, 1), length(alucons));
    int s0 = toViewPosition(row(align, 0), 0);
    int s1 = toViewPosition(row(align, 0), length(ref_fa));
    sim_rate = 0;
    int align_len = min(l1, s1) - max(l0, s0);
    //cout << align << " " << align_len << endl;
    if ( align_len <= CLIP_BP or align_len <= length(ref_fa) * sim_th or align_len >= length(ref_fa) / sim_th )
	return false;
    TRow &row0 = row(align,0);
    TRowIterator it0 = begin(row0);
    TRow &row1 = row(align,1);
    TRowIterator it1 = begin(row1);
    int i = 0, dif = 0;
    while ( i++ < max(l0, s0) )
    {
	it0++;
	it1++;
    }
    while ( i++ <= min(l1, s1) )
    {
	if ( (*it0) != (*it1) ) dif++;     ////if(isGap(it1))
	it0++;
	it1++;
    }
    sim_rate = 1 - dif / (float) align_len;

    //if (sim_rate) cout << dif << " " << align << endl;

    return  sim_rate >= sim_th;
}

string align_alu_cons_call(seqan::CharString & ref_fa, AluconsHandler *alucons_fh, float & sim_rate, float sim_th, bool read_is_clipped)
{
    for ( vector<string>::iterator si = (alucons_fh->seq_names).begin(); si != (alucons_fh->seq_names).end(); si++)
    {
	alucons_fh->update_seq_name(*si);
	for (int k = 1; k <= 4; k++) 
	    if (align_alu_cons(ref_fa, alucons_fh->fetch_alucons(k), sim_rate, sim_th, read_is_clipped))
		return *si;
    }

    return "";
}

bool covered_reads(BamFileHandler * bam_fh, string chrn, int p1, int p2, int minCnt)
{
    if (! bam_fh->jump_to_region(chrn, p1, p2) ) return false;
    seqan::BamAlignmentRecord record;
    int cnt = 0;
    while ( true )
    {
	string read_status = bam_fh->fetch_a_read(chrn, p1, p2, record);
	if (read_status == "stop" ) break;
	if (read_status == "skip" or !QC_insert_read(record)) continue; 
	if (record.rID != record.rNextId) continue;
	if (abs(record.tLen) > DISCORDANT_LEN )
	{
	    cnt++;
	}
	else
	{
	    int r1 = min ( record.beginPos, record.pNext);
	    int r2 = r1 + abs(record.tLen) ;
	    if ( min(p2, r2) - max (r1, p1) > 0 )
		cnt++;
	}
	if ( cnt >= minCnt) return true;
    }
    return cnt >= minCnt;
}


bool combine_pns_vcf(float read_dist_th, string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, map <string, std::set<int> > & chrn_aluBegin, float llh_th, string ref_name)
{
    vector <string>::iterator ci, pi;
    //seqan::VcfStream vcfout("-", seqan::VcfStream::WRITE);
    seqan::VcfStream vcfout(f_out.c_str(), seqan::VcfStream::WRITE);
    for ( ci = chrns.begin(); ci != chrns.end(); ci++)
	appendValue(vcfout.header.sequenceNames, *ci);
    for ( pi = pns.begin(); pi != pns.end(); pi++) 
	appendValue(vcfout.header.sampleNames, *pi);
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileformat", "VCFv4.1"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("fileDate", "201402"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("reference", ref_name));  // eg: hg19
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=NS,Number=1,Description=\"Number of Samples With Data\">"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=AF,Number=A, Description=\"Allele Frequency\">"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("INFO", "<ID=PA,Number=A, Description=\"Evidence of insertions based on alu mate (otherwise from clip reads)\">"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=highCov, Description=\"Region coverage too high\">"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FILTER", "<ID=chisq1, Description=\"likelihood low, failed chisq df = 1\">"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FORMAT", "<ID=GT,Number=1,Type=Integer, Description=\"Genotype\">"));
    appendValue(vcfout.header.headerRecords, seqan::VcfHeaderRecord("FORMAT", "<ID=PL,Number=3,Type=Integer, Description=\"Phred-scaled likelihoods for genotypes\">"));

    seqan::VcfRecord record;    
    record.ref = "0";
    record.alt = "1";
    record.format = "GT:PL";

    ifstream fin;
    stringstream ss;
    string line, chrn, aluEnd;
    int aluBegin, cii, flag;
    int _midCnt_alu, _midCnt_clip, clipCnt, unknowCnt;
    map <int, int> midCnts_alu, midCnts_clip; 
    float p0, p1, p2;      
    string phred_str_00 = "0,100,255";  // min allowed 
    string phred_str_missing = "0,0,0";
    map < pair<int, string>, map <string, GENO_PROB > > pos_pnInfo;
    for ( cii = 0, ci = chrns.begin(); ci != chrns.end(); ci++, cii++)
    {
	pos_pnInfo.clear();
	for ( pi = pns.begin(); pi != pns.end(); pi++)
	{
	    fin.open( (path0 + *pi + f_in_suffix).c_str() );
	    if ( !fin )
	    {
		cerr << "ERROR: " << path0 + *pi + f_in_suffix << "not exist!" << endl; 
		exit(0);
	    }
	    getline(fin, line); 
	    flag = 1;
	    while ( getline(fin, line) )
	    {
		ss.clear(); ss.str( line );
		ss >> chrn;
		if (chrn != *ci)
		{
		    if (flag) continue;
		    else break;
		}
		flag = 0; 
		ss >> aluBegin >> aluEnd >> _midCnt_alu;
		if ( _midCnt_alu < 0 )
		{  // evidence of G00
		    GENO_PROB one_gi = GENO_PROB(1, 0, 0, phred_str_00, 0); 
		    pos_pnInfo[ make_pair(aluBegin, aluEnd) ].insert ( std::pair<string, GENO_PROB>(*pi, one_gi) );
		}
		else
		{
		    ss >> _midCnt_clip >> clipCnt >> unknowCnt >> p0 >> p1 >> p2 ;
		    if ( !key_exists(midCnts_alu, aluBegin) )
		    {
			midCnts_alu[aluBegin] = 0;
			midCnts_clip[aluBegin] = 0;
		    }
		    midCnts_alu[aluBegin] += _midCnt_alu;
		    midCnts_clip[aluBegin] += _midCnt_clip;
		    string phred_scaled_str = phred_scaled(p0, p1, p2);
		    int _g = 0;
		    if (p1 > p0 or p2 > p0) _g =  ( p1 > p2) ? 1 : 2;
		    GENO_PROB one_gi = GENO_PROB(p0, p1, p2, phred_scaled_str, _g); 
		    pos_pnInfo[ make_pair(aluBegin, aluEnd) ].insert ( std::pair<string, GENO_PROB>(*pi, one_gi) );
		}
	    }
	    fin.close();
	}

	map < pair<int, string>, map<string, GENO_PROB> >::iterator pi3;
	map < string, GENO_PROB>::iterator pi2;
	for (pi3 = pos_pnInfo.begin(); pi3 != pos_pnInfo.end(); pi3++)
	{
	    int altCnt = 0;
	    for ( pi2 = (pi3->second).begin(); pi2 != (pi3->second).end(); pi2++ ) 
		altCnt += (pi2->second).geno;
	    if (!altCnt) continue;
	    int _pos = (pi3->first).first;
	    int _tcnt = midCnts_clip[_pos] + midCnts_alu[_pos];
	    //if ( _tcnt < 5 ) continue;  // if total evidence is too small, ignore this position. fixme: comment this line if sample size small 
	    float alu_evidence =  (float) midCnts_alu[_pos] / (float) _tcnt;
	    if ( alu_evidence < read_dist_th or 1 - alu_evidence < read_dist_th ) 
		continue;
	    ////cout << (pi3->first).first << " " << altCnt << " " << alu_evidence << " "  << _tcnt << endl;      
	    record.beginPos = (pi3->first).first;
	    record.rID = cii;   // change ???
	    string debugInfo = (pi3->first).second;
	    record.id = debugInfo;
	    int n_missing = 0;
	    for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++)
	    {
		pi2 = pi3->second.find(*pi);
		string ginfo;
		if ( pi2 == pi3->second.end() )
		{
		    ginfo = "0/0:" + phred_str_missing ;
		    n_missing += 1;
		}
		else
		{ 
		    ginfo = reformat_geno( (pi2->second).geno ) + ":" + (pi2->second).phredStr;
		}
		appendValue(record.genotypeInfos, ginfo);      
	    }

	    record.qual = 0;
	    record.filter = "PASS";
	    if ( (pi3->first).second.find(',') != string::npos )
	    {
		string exact_left, exact_right;
		split_by_sep(debugInfo, exact_left, exact_right, ',');  
		if (exact_left == "0" or exact_right == "0")
		    record.filter = "BreakpointOneside";
	    }

	    if ( chrn_aluBegin[*ci].find( (pi3->first).first ) != chrn_aluBegin[*ci].end() )
	    {
		record.filter = "highCov";
		record.info = "NS=.;AF=.;PA=.";
	    }
	    else
	    {
		stringstream ss_record_info;
		int n_pass = pns.size() - n_missing;
		float altFreq = altCnt / 2./ float ( n_pass );
		ss_record_info << "NS=" << n_pass << ";AF=" << setprecision(4) << altFreq << ";PA=" << setprecision(4) << alu_evidence;
		record.info = ss_record_info.str();
		float llh_alt = 0, llh_00 = 0;
		float freq0 = (1 - altFreq) * (1 - altFreq);
		float freq1 = 2 * altFreq * (1 - altFreq);
		float freq2 = altFreq * altFreq;
		for (vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++)
		{
		    pi2 = pi3->second.find(*pi);
		    if ( pi2 == pi3->second.end() ) // missing data 
			continue;
		    llh_00 +=  ( (pi2->second).g0 > 0 ? log( (pi2->second).g0 ) : ( - LOG10_GENO_PROB ) ) ;
		    llh_alt += log (freq0 * (pi2->second).g0 + freq1 * (pi2->second).g1 + freq2 * (pi2->second).g2);
		}
		record.qual = llh_alt - llh_00; 
		if (record.qual  <= llh_th)  
		    record.filter = "chisq1";
	    }

	    writeRecord(vcfout, record);
	    clear(record.genotypeInfos);
	}
	cout << *ci << " size " << pos_pnInfo.size() << endl;
    }
    // chrn finished
    clear(record);
    seqan::close(vcfout);
    return true;
}

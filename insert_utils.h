#define SEQAN_HAS_ZLIB 1
#include "utils.h"
#ifndef INSERT_UTILS_H
#define INSERT_UTILS_H

typedef map<int, ofstream* > MapFO;
typedef map<string, ofstream* > MapSFO;
typedef seqan::Dna5String TSeq;

struct ParametersAluInsert : public ParametersAlu
{
    int left_plus_right, coverage_cnt, scan_win_len, max_len_region,
	clip_bp_mid, alucons_len, clip_q;
    float consensus_freq, pa_th;
    ParametersAluInsert():
	left_plus_right(4), // if coverage is small, change it to 3
	coverage_cnt(5), // if coverage >= 5 AND no clip reads ==> force genotype 00
	scan_win_len(400),  // 0.99 quantile. 450 for 0.999 quantile. only approximate,different reading group differs
	max_len_region(100),  // max length of insert region
	clip_bp_mid(30), 
	alucons_len(300), //  default alu length
	clip_q(25),    //cut reads on both ends 
	consensus_freq(0.7),  
	pa_th(0.01)
    {}
    // explained in vcf info 
};

class GetIOPathInsert : public GetIOPath{
    public:
	string path0_, path1_, pathClip_, pathCons_, pathDel0_;
	string bam_rid_chrn_;
	GetIOPathInsert(string config_file);
};

inline void close_fhs(MapFO & fileMap, map<int, string> & rID_chrn)
{
    for (map<int, string>::iterator rc = rID_chrn.begin(); rc != rID_chrn.end(); rc++)   
	delete fileMap[rc->first];
}

inline int get_min_pair( int pa, int pb )
{
    if (pa <= 0) return pb;
    if (pb <= 0) return pa;
    return min(pa, pb);
}

inline pair<int, int> get_valid_pair(int pa, int pb)
{
    if ( pa<=0 ) return make_pair(pb, pb);
    if ( pb<=0 ) return make_pair(pa, pa);
    return  make_pair(pa, pb);
}

inline void add3key(map<int, int> & confidentPos, int k)
{
    confidentPos[k] = 1;
    confidentPos[k-1] = 1;
    confidentPos[k+1] = 1;
}

class AlumateINFO
{
    public:
	string qname;
	int len_read, rid1, rid2, pos1, pos2, rgIdx;
	bool RC1, RC2;
	AlumateINFO( string & qn, int lr, int id1, int id2, int p1, int p2, int rgIdx, bool r1, bool r2) 
	    : qname(qn), len_read(lr), rid1(id1), rid2(id2), pos1(p1), pos2(p2), rgIdx(rgIdx), RC1(r1), RC2(r2)
	{}
	static bool sort_pos2(const AlumateINFO* a, const AlumateINFO* b);
	static void delete_list(list <AlumateINFO *> & alumate_list);
};

class READ_INFO
{
    public:
	int beginPos, endPos;
	bool should_be_left; // inferred from RC flag
	string alu_type;
	READ_INFO(int p, int lr, bool sbl, string alu_type) : beginPos(p), endPos(p+lr-2), should_be_left(sbl), alu_type(alu_type)
    {}
};

class ALUREAD_INFO
{
    public:
	seqan::CharString qName;
	int this_pos, clipLeft, pos, readLen;
	bool sameRC;
	ALUREAD_INFO(seqan::CharString & qn, int tp, int cl, int p, int rlen, bool t) : qName(qn), this_pos(tp), clipLeft(cl), pos(p), readLen(rlen), sameRC(t)
    { }
};

bool clipRight_move_left(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len);
bool clipLeft_move_right(seqan::CharString & read_seq, seqan::CharString & ref_fa, list <int> & cigar_cnts, int refBegin, int & clipPos, int & align_len);

bool read_first2col(string fn, vector < pair<int, int> > & insert_pos, bool has_header);
void normalize_llh(float *loggp, float th); 
int parseline_ins(string line0, ostream & fout, map <int, EmpiricalPdf *> & pdf_rg, float logPE, int estimatedAluLen, int err_code, bool test_print, float *test_gp = NULL);
int parseline_cnt(string line0);
void align_clip_to_consRef(string shortSeq, string longSeq, int & refBegin, int & refEnd, int clipLen);
bool align_alu_to_consRef(const string & shortSeq, const string & longSeq, float dif_th, string loginfo) ;

bool align_alu_cons(seqan::CharString &ref_fa, seqan::CharString alucons, float & sim_rate,float sim_th, bool read_is_clipped);
string align_alu_cons_call(seqan::CharString & ref_fa, AluconsHandler *alucons_fh, float & sim_rate, float sim_th, bool read_is_clipped = false);
bool covered_reads(BamFileHandler * bam_fh, string chrn, int p1, int p2, int minCnt); 

// following, depreciated
void filter_outlier_pn(string path_input, string fn_suffix, map<int, string> &ID_pn, string chrn, string file_pn_used_output, float percentage_pn_used);
bool combine_pns_vcf(float read_dist_th, string path0, string f_in_suffix, string f_out, vector <string> &pns, vector <string> & chrns, map <string, std::set<int> > & chrn_aluBegin, float llh_th, string ref_name);
#endif /*INSERT_UTILS_H*/

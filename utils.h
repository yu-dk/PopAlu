// some data struct to read and write data 
#define SEQAN_HAS_ZLIB 1
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include <seqan/vcf_io.h>
#include "common.h"
#ifndef UTILS_H
#define UTILS_H

typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
typedef seqan::Align<seqan::CharString, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow; 
typedef seqan::Iterator<TRow>::Type TRowIterator;

struct MyUpper : public unary_function<char,char>
{
    inline char operator()(char x) const 
    {
	if (('a' <= x) && (x <= 'z')) return (x + ('A' - 'a'));
	return x; 
    }
};

class GENO_PROB
{
    public:
	float g0, g1, g2;
	string phredStr;
	int geno;
	GENO_PROB(float g0, float g1, float g2, string ps, int _g) : g0(g0), g1(g1), g2(g2), phredStr(ps), geno(_g)
    {}
};

class GetIOPath
{
    public:
	ConfigFileHandler * cf_fh;
	string bam_path_, output_path_, insertlen_, insertlen_dist_, insertlen_count_;
	GetIOPath(string config_file);
	~GetIOPath()
	{ delete cf_fh;}

	string GetBamFileName(string pn)
	{ return bam_path_ + pn + ".bam"; }
	string CheckFolder(string path);
};

struct ParametersAlu{
    string pdf_param;
    float log10_ratio_ub;
    int mid_cov_cnt, alu_flank;
    ParametersAlu(): pdf_param("150_1000_10"), 
    //ParametersAlu(): pdf_param("100_1000_5"), 
    log10_ratio_ub(3.5),
    mid_cov_cnt(10), // if clip reads >= this number, no need to align to consensus alu
    alu_flank(600)

    {}

}; 

inline string get_name_rg(string prefix, string pn)
{ return prefix + "RG." + pn;}
inline string get_name_rg_pdf(string prefix, string pn, string rg, string pdf_param)
{ return prefix + pn + ".count." + rg + "." + pdf_param; }
// !RC
inline bool left_read( seqan::BamAlignmentRecord &record)
{return (record.beginPos < record.pNext);}

inline bool QC_read( seqan::BamAlignmentRecord &record)
{  
    return (!hasFlagQCNoPass(record)) and (!hasFlagDuplicate(record)) and hasFlagMultiple(record) and (!hasFlagSecondary(record));
}

inline bool QC_delete_read( seqan::BamAlignmentRecord &record)
{  
    return QC_read(record) and (record.rID == record.rNextId) and hasFlagAllProper(record) and (abs(record.tLen) <= DISCORDANT_LEN) and (hasFlagRC(record) != hasFlagNextRC(record)) and ( left_read(record) != hasFlagRC(record)); 
};

inline bool QC_insert_read( seqan::BamAlignmentRecord &record)
{  
    return QC_read(record) and (!hasFlagUnmapped(record)) and (!hasFlagNextUnmapped(record));
};

inline bool qname_match( seqan::BamAlignmentRecord &record, string qn)
{  
    return qn == toCString(record.qName);
};

inline int qualToInt( char c )
{ return (int) c -33; };

inline bool has_soft_last(seqan::BamAlignmentRecord &record, unsigned min_bp)
{ 
    unsigned i = length(record.cigar) - 1;
    return (record.cigar[i].operation == 'S') and (record.cigar[i].count >= min_bp) ;
};

inline bool has_soft_first(seqan::BamAlignmentRecord &record, unsigned min_bp)
{ 
    return (record.cigar[0].operation == 'S') and (record.cigar[0].count >= min_bp); 
};

// check if the RC flag is valid, if a pair of reads span insertion region
inline bool aluclip_RC_match(seqan::BamAlignmentRecord &record, bool s1, bool s2)
{
    if ( (record.rID != record.rNextId ) or abs(record.tLen) > DISCORDANT_LEN )
    {
	if ( s1 and !hasFlagRC(record))   // left of breakpoint
	    return true;
	if ( s2 and hasFlagRC(record))
	    return true;
	return false;
    }
    return true;
}

inline bool is_aluflag_cnt(seqan::BamAlignmentRecord &record)
{ 
    int nc = length(record.cigar);
    if ( has_soft_last(record, CLIP_BP) and has_soft_first(record, CLIP_BP))
	return false;
    if ( record.cigar[0].operation == 'S' and record.cigar[0].count > length( record.seq ) - 50)
	return false;
    if ( record.cigar[nc-1].operation == 'S' and record.cigar[nc-1].count > length( record.seq ) - 50)
	return false;
    int non_match_len = 0;
    for (int i = 1; i < nc - 1 ; i++) 
	if (record.cigar[i].operation != 'M' ) 
	    non_match_len += record.cigar[i].count; 
    return non_match_len <= 5;
}

inline int count_non_match(seqan::BamAlignmentRecord &record)
{ 
    int non_match_len = 0;
    for (size_t i = 0; i < length(record.cigar); i++) 
	if (record.cigar[i].operation != 'M') non_match_len += record.cigar[i].count; 
    return non_match_len;
}

inline bool p00_is_dominant(float * log10_p, int min_log10p)
{ return  max( log10_p[2] - log10_p[0], log10_p[1] - log10_p[0]) <= min_log10p; }
inline bool p11_is_dominant(float * log10_p, int min_log10p)
{ return  max( log10_p[0] - log10_p[2], log10_p[1] - log10_p[2]) <= min_log10p; }
inline bool cnts_confident(int cnt)
{ return cnt >=5 ;}
inline string phred_log (float p)
{ 
    if (!p) return "255";
    //float tmpf = -(log10 (p) * 10); 
    float tmpf = max( (float) 1.0, - (log10 (p) * 10)) ; 
    return int_to_string( (int) tmpf );
}

class EmpiricalPdf
{
    map <int, float> prob_vec; 
    float min_prob;
    public:
    int min_len, max_len, bin_width;
    bool RG_mate_pair;
    EmpiricalPdf(string pdf_file, string rg_name, map <string, int> & isPairEnd); // if median insert length > len_ub, ignore this RG 
    float pdf_obs(int insertlen);
    void ratio_obs(int y, int z, float log10_ratio_ub, float & py, float & pz);
    static void delete_map(map <int, EmpiricalPdf *> & epdf_rg); 
};

class BamFileHandler{
    public:
	string fn_bai;
	string fn_bam_out;
	map<int, string> rID_chrn ;
	BamFileHandler(vector<string> &chrns, string bam_input, string bai_input, string bam_output="");
	~BamFileHandler(void);
	bool jump_to_region(int rid, int region_begin, int region_end);
	bool jump_to_region(string chrn, int region_begin, int region_end);
	bool fetch_a_read(seqan::BamAlignmentRecord & record);
	string fetch_a_read(string chrn, int region_begin, int region_end, seqan::BamAlignmentRecord & record);
	bool get_chrn(int query_rid, string & pairChrn);
	void print_mapping_rID2chrn();
	bool write_a_read(seqan::BamAlignmentRecord & record);  
	static BamFileHandler * openBam_24chr(string bam_input, string bai_input="", string bam_output=""); // chr1 - chrY

    private:
	TNameStore  nameStore;
	TNameStoreCache nameStoreCache;
	TBamIOContext context;
	seqan::Stream<seqan::Bgzf> inStream ;
	map<seqan::CharString, int> chrn_rID ;
	seqan::BamIndex<seqan::Bai> baiIndex ;    
	seqan::Stream<seqan::Bgzf> outStream ;
};

class FastaFileHandler
{
    public:
	seqan::FaiIndex faiIndex;
	unsigned idx;
	bool only_one_seq;
	FastaFileHandler(string fn_fa);
	FastaFileHandler(string fn_fa, string seq_name);
	void fetch_fasta_upper(int beginPos, int endPos, seqan::CharString & seq);
	void fetch_fasta_upper(int beginPos, int endPos, seqan::CharString & seq, string seq_name);
	static seqan::CharString fasta_seq(string fa_input, string seq_name,int beginPos, int endPos);
};


class AluconsHandler : public FastaFileHandler
{
    public:
	vector <string> seq_names;
	string seq_name;
	int seq_len;
	AluconsHandler(string fn_fa, string sn = "AluY");
	void update_seq_name(string sn);
	seqan::CharString fetch_alucons(int key);
    private:
	map <int, seqan::CharString> seqs;
};

class AluRefPosRead // used only by alu_hg18
{
    queue<int> beginP, endP;
    queue<char> strandP;
    queue<string> aluType;
    public:
    int minP, maxP;
    AluRefPosRead(string file_alupos, int minLen = 200);
    int updatePos(int &beginPos, int &endPos, char & chain, string & alu_type);
};

class RepDB1
{
    public:
	int begin_pos, end_pos;
	string alu_type;
	RepDB1(int bp, int ep, string at):begin_pos(bp), end_pos(ep), alu_type(at)
    {}
	bool operator<(const RepDB1 & other)const
	{ return begin_pos < other.begin_pos; }

};

class AluRefPos   // used by alu_delete, alu_insert
{
    public: 
	int db_size;
	AluRefPos(string fn, int minLen_alu);
	AluRefPos(string fn, int minLen_alu, int minDist_neighbor);   
	bool nextdb()
	{ adi++; return adi != alu_db.end();}
	int get_beginP() const
	{ assert(adi != alu_db.end()); return (*adi).begin_pos; }
	int get_endP() const
	{ assert(adi != alu_db.end()); return (*adi).end_pos; }
	string get_type() const
	{ assert(adi != alu_db.end()); return (*adi).alu_type; }
	bool within_alu(int pos);
	void debug_print();
	~AluRefPos(void)
	{ alu_db.clear();}
    private:
	std::set <RepDB1> alu_db;
	std::set <RepDB1>::iterator adi;  
};


class RepMaskPos
{
    public:
	map<string, vector<int> > beginP, endP;
	RepMaskPos(string file_rmsk, vector<string> &chrns, int join_len = 20);
	~RepMaskPos(void);
};

void get_RG_info(string fn, string pn, map<string, int > & isPairEnd);  /* wether this reading group is paired end or mate pair */
seqan::Pair<seqan::CigarElement<>::TCount> mappedInterval(seqan::String<seqan::CigarElement<> > & cigar);
double avgQuality(seqan::CharString & qual, seqan::Pair<seqan::CigarElement<>::TCount> & interval);
bool QC_insert_read_qual( seqan::BamAlignmentRecord &record);
string replace_str0_str(string fn, string chrn, string chr0="chr0");
void read_pdf_pn( string prefix, string pn, string pdf_param,  map <int, EmpiricalPdf *> & pdf_rg, map<string, int> &isPairEnd);
void get_min_value(map <int, float> & m, float & min_val, int & min_key);
void log10P_to_P(float *log_gp, float *gp, int max_logp_dif);
string phred_scaled(float p0, float p1, float p2);
void parse_reading_group(string file_rg, map<string, int> & rg_to_idx);

int get_rgIdx(map<string, int> & rg_to_idx, seqan::BamAlignmentRecord & record);
void parse_cigar(string cigar, list <char> & opts, list <int> & cnts);
string get_cigar(seqan::BamAlignmentRecord &record);
void debug_print_read(seqan::BamAlignmentRecord &record, ostream & os = cout);
bool find_read(string &bam_input, string &bai_input, string &chrn, string &this_qName, int this_pos, seqan::BamAlignmentRecord &that_record, int flank_region);

bool get_trim_length( seqan::BamAlignmentRecord & record, int & trimb, int & trime, int bpQclip);
bool trim_clip_soft_first( seqan::BamAlignmentRecord & record, int & clipLen, seqan::CharString & clipSeq, int bpQclip);  
bool trim_clip_soft_last( seqan::BamAlignmentRecord & record, int & clipLen, seqan::CharString & clipSeq, int bpQclip);

int numOfBestHits(seqan::BamAlignmentRecord &record);
void get_chrn(string fn, map<int, string> & rid_chrn);


#endif /*UTILS_H*/

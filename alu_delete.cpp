#define SEQAN_HAS_ZLIB 1
#include "delete_utils.h"

ParametersAluDelete p_alu;

void count_reads(map <seqan::CharString, T_READ> &qName_info, map < T_READ, int > &readCnt)
{
	readCnt.clear();
	for (map <seqan::CharString, T_READ>::iterator rt = qName_info.begin(); rt != qName_info.end(); rt++) 
		addKey(readCnt, rt->second); 
}

int check_one_pos(BamFileHandler* bam_fh, FastaFileHandler *fasta_fh, map <string, int> &rg_to_idx, string chrn, int aluBegin, int aluEnd, int &totalCnt, map <seqan::CharString, T_READ> &qName_info,  map<seqan::CharString, string> & rg_str)
{  
	qName_info.clear();
	int reads_cnt = 0;
	seqan::BamAlignmentRecord record;  
	while ( true )
	{
		string read_status = bam_fh->fetch_a_read(chrn, aluBegin - ::p_alu.alu_flank, aluEnd + ::p_alu.alu_flank, record);
		if (read_status == "stop" ) break;
		if (read_status == "skip" or !QC_delete_read(record)) continue;
		reads_cnt ++;  
		map <seqan::CharString, T_READ >::iterator qItr = qName_info.find(record.qName);
		if ( qItr != qName_info.end() and qItr->second != unknow_read) 
			continue;
		T_READ iread = classify_read( record, aluBegin, aluEnd, fasta_fh);   
		if ( iread == useless_read) continue;    
		//    if ( iread == clip_read) 
		//      cout << "debug " << aluBegin << " " << get_cigar(record) << " " << record.qName << endl;

		qName_info[record.qName] = iread;
		if ( qItr == qName_info.end() and iread == unknow_read)
		{
			int rgIdx = get_rgIdx(rg_to_idx, record);
			if (rgIdx >= 0)
			{
				stringstream rg_ss;
				rg_ss << rgIdx << ":" << abs(record.tLen);
				rg_str[record.qName] = rg_ss.str();
			}
		}

	}
	totalCnt = reads_cnt;
	if ( length(record.seq) * (float) reads_cnt / (aluEnd - aluBegin + 2 * ::p_alu.alu_flank) > ::p_alu.coverage_high) return ::p_alu.coverage_high;
	return 1;
}

int delete_search(BamFileHandler *bam_fh, string file_fa_prefix, vector<string> &chrns, string &f_out, string &f_log, string &file_alupos_chr0, map<string, int> &rg_to_idx)
{    
	map < seqan::CharString, T_READ> qName_info;  
	ofstream f_tmp1( f_out.c_str()); 
	f_tmp1 << "chr aluBegin aluEnd totalCnt midCnt clipCnt unknowCnt unknowStr\n";
	ofstream f_log1( f_log.c_str());  // print out info for clip reads 

	for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++)
	{
		string chrn = *ci;
		string file_alupos = replace_str0_str(file_alupos_chr0, chrn, "chr0");
		//cout << "debug " << chrn << " " << file_alupos << endl;
		AluRefPos *alurefpos = new AluRefPos(file_alupos, ::p_alu.minLen_alu_del, 300);  // min distance to neighbor is 300 bp 
		FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + chrn + ".fa", chrn);    
		int aluBegin, aluEnd;
		for (int count_loci = 0; ; count_loci++)
		{
			int totalCnt = 0;
			if ( !alurefpos->nextdb() ) break;      
			aluBegin = alurefpos->get_beginP();
			aluEnd = alurefpos->get_endP();
			//cout << "debug " << aluBegin << " " << aluEnd << endl;

			if (aluBegin <= ::p_alu.alu_flank or !bam_fh->jump_to_region(chrn, aluBegin-::p_alu.alu_flank, aluEnd + ::p_alu.alu_flank))
				continue;
			map<seqan::CharString, string> rg_str;
			int check_alu = check_one_pos(bam_fh, fasta_fh, rg_to_idx, chrn, aluBegin, aluEnd, totalCnt, qName_info, rg_str);
			if ( check_alu == 0 ) continue;
			if ( check_alu == ::p_alu.coverage_high)
			{
				f_log1 << "COVERAGE_HIGH " << chrn << " " << aluBegin << " " << aluEnd << " " << totalCnt << endl;
				continue;
			}
			map < T_READ, int > readCnt;
			count_reads(qName_info, readCnt); 
			if ( readCnt[clip_read] or readCnt[unknow_read] )
			{// not interesting if all are mid_reads
				f_tmp1 << chrn << " " << aluBegin << " " << aluEnd << " " << totalCnt << " " <<  readCnt[mid_read]
					<< " " << readCnt[clip_read] << " " << readCnt[unknow_read];
				for (map<seqan::CharString, string>::iterator ri = rg_str.begin(); ri != rg_str.end(); ri++) 
					if (qName_info[ri->first] == unknow_read)
						f_tmp1 << " " << ri->second ;
				f_tmp1 << endl;
			}
			else if (totalCnt > 0)
				f_tmp1 << chrn << " " << aluBegin << " " << aluEnd << " " <<  -totalCnt << endl;  // diff from missing data 
		}
		delete alurefpos;
		delete fasta_fh;
		cout << "file_alupos:done  " << file_alupos << endl;  
	}
	f_tmp1.close();
	f_log1.close();
	return 0;
}

int parseline_del(ostream & fout, string &line, map <int, EmpiricalPdf *> & pdf_rg, float logPE, int err_code, bool test_print = false)
{
	const int read_len = 100;
	float *log10_gp = new float[3];
	stringstream ss;
	string chrn;
	int aluBegin, aluEnd, totalCnt, midCnt, clipCnt, unknowCnt;
	ss.clear(); ss.str(line); 
	ss >> chrn >> aluBegin >> aluEnd >> totalCnt;
	if (totalCnt < 0)
	{    // G00 > 1e-5, not missing data
		fout << chrn << " " << aluBegin << " " << aluEnd << " " << totalCnt << endl;  
		return 0;
	}

	// clipCnt is evidence if H1
	ss >> midCnt >> clipCnt >> unknowCnt ;
	if ( clipCnt >= ::p_alu.mid_cov_cnt and midCnt >= ::p_alu.mid_cov_cnt)
	{
		fout << chrn << " " << aluBegin << " " << aluEnd << " " << midCnt << " " << clipCnt << " " << unknowCnt << " 0 1 0\n";
		return 0;
	}

	logPE = - abs(logPE);
	float logPM = log10 ( 1 - pow(10, logPE) );
	float bias1 = (aluEnd - aluBegin + 3. * read_len) / (3. * read_len); 
	float ph0 = bias1 / (bias1 + 1);
	if (err_code == 1) ph0 = 0.85;   // 0/0 should be 0/1
	else if (err_code == 2) ph0 = 0.5;   // 1/1 should be 0/1
	else if (err_code == 3) logPE = -2;  // 0/1 should be 1/1 
	if (test_print) cout << "ph0 "<< ph0 << endl;

	float *gp = new float[3];
	for (int i = 0; i < 3; i++) log10_gp[i] = 0;
	if (midCnt + clipCnt > 0)
	{
		log10_gp[0] = clipCnt * logPE + midCnt * logPM;
		log10_gp[1] = midCnt * log10 (ph0) + clipCnt * log10 (1 - ph0) + (midCnt + clipCnt) * logPM  ; 
		log10_gp[2] = midCnt * logPE + clipCnt * logPM;
	}
	if (test_print ) cout << "debug1# " << log10_gp[0] << " " << log10_gp[1] << " " << log10_gp[2] << endl;   

	int _unknowCnt = 0;
	if (unknowCnt)
	{ 
		int insert_len, idx;
		string token;
		for (int i = 0; i < unknowCnt; i++)
		{
			getline(ss, token, ':');
			seqan::lexicalCast2(idx, token);
			getline(ss, token, ' ');
			seqan::lexicalCast2(insert_len, token);      
			if (pdf_rg[idx]->RG_mate_pair)
				continue;

			_unknowCnt++;
			float p_y, p_z;
			float down_weight = abs(logPE) - 0.5;
			if ( err_code == 3 )  down_weight = abs(logPE) - 1.5;
			pdf_rg[idx]->ratio_obs(insert_len, insert_len - aluEnd + aluBegin, down_weight, p_y, p_z);      
			log10_gp[0] += log10 (p_y);
			log10_gp[1] += log10 (ph0 * p_y + (1 - ph0) * p_z) ;
			log10_gp[2] += log10 (p_z);
		}
	}

	if (test_print ) cout << "debug2# " << log10_gp[0] << " " << log10_gp[1] << " " << log10_gp[2] << endl;   
	if ( clipCnt >= 3 and log10_gp[0] > max(log10_gp[1], log10_gp[2]) )
	{ // call 0/0 when clipCnt exists
		if (err_code == 0 ) return 1;
		cerr << "##error1 " << line << endl;
	}

	if ( midCnt >= 3 and log10_gp[2] > max(log10_gp[1], log10_gp[0]) )
	{ // call 1/1 when midCnt exists
		if (err_code == 0 ) return 2;
		cerr << "##error2 " << line << endl;
	}

	if ( clipCnt >= 3 and (clipCnt + 0.1) / (midCnt + 0.1) > 3 and log10_gp[2] < max ( log10_gp[1], log10_gp[0] ) )
	{ // should be 1/1
		if (err_code == 0 ) return 3;
		cerr << "##error3 " << line << endl;    
	}

	if ( !p00_is_dominant(log10_gp, - LOG10_GENO_PROB) )
	{
		fout << chrn << " " << aluBegin << " " << aluEnd << " " << midCnt << " " << clipCnt << " " << _unknowCnt ;
		log10P_to_P(log10_gp, gp, LOG10_GENO_PROB);  // normalize such that sum is 1
		fout << " " << setprecision(6) << gp[0] << " " << gp[1] << " " << gp[2] << endl;   
	}
	else
	{
		fout << chrn << " " << aluBegin << " " << aluEnd << " " << - totalCnt + unknowCnt - _unknowCnt << endl;  
	}
	delete gp;    
	delete log10_gp;
	return 0;
}

void read_highCov_region(string f_input, map < string, set<int> > & chrn_aluBegin)
{
	ifstream fin(f_input.c_str());
	string tmp1, tmp2, chrn;
	int aluBegin, aluEnd;
	while ( fin >> tmp1>> chrn >> aluBegin >> aluEnd >> tmp2 ) chrn_aluBegin[chrn].insert(aluBegin);
	fin.close();
}

void write_rm1(string f_output, map < string, std::set<int> > & chrn_aluBegin)
{
	ofstream fout(f_output.c_str());
	for (map < string, std::set<int> >::iterator it = chrn_aluBegin.begin(); it != chrn_aluBegin.end(); it++) 
		for ( std::set<int>::iterator i2 = (it->second).begin(); i2 != (it->second).end(); i2++ )
			fout << it->first << " " << *i2 << " highCov\n";
	fout.close();
}

int main( int argc, char* argv[] )
{
	::p_alu = ParametersAluDelete();
	string config_file = argv[1];
	string opt = argv[2];
	if (argc < 3) exit(1);

	boost::timer clocki;    
	clocki.restart();  
	ConfigFileHandler cf_fh = ConfigFileHandler(config_file);
	GetIOPathDelete *io_path = new GetIOPathDelete(config_file);

	vector<string> chrns;
	string s_chrns = cf_fh.get_conf("chrns");
	parse_chrns(s_chrns, chrns);
	string file_fa_prefix = cf_fh.get_conf("file_fa_prefix");
	string file_alupos_prefix = cf_fh.get_conf("file_alupos_prefix"); 
	check_folder_exists(file_alupos_prefix);

	string path0 = io_path->path0_;
	string file_dist_prefix = io_path->insertlen_dist_;

	if ( opt == "0" )
	{  
		// sort alu files 
		for (vector<string>::iterator ci = chrns.begin(); ci!= chrns.end(); ci++)   // combine Alu if less than 10 
			sort_file_by_col<int> (file_alupos_prefix + "alu_" + *ci, 2, false);
		check_folder_exists(path0 + "log1s/");
		check_folder_exists(path0 + "tmp1s/");
		check_folder_exists(path0 + "tmp2s/");
	}
	else if ( opt == "1" )
	{  
		// step 1: detect each alu deletions for each individual 
		// USAGE: opt3/alu_delete [configuration file] 1 [sample name]

		if (argc != 4) 
			ErrorMsg1("[sample name] is missing", "alu_delete");

		string pn = argv[3];  //// string pn = ID_pn[ seqan::lexicalCast<int> (argv[3])];
		map<string, int> rg_to_idx;
		parse_reading_group( get_name_rg(file_dist_prefix, pn), rg_to_idx );

		string fn_tmp1 = path0 + pn + ".tmp1";
		string fn_log1 = path0 + pn + ".log1";
		string bam_input = io_path->GetBamFileName(pn);
		string bai_input = bam_input + ".bai";      
		BamFileHandler *bam_fh = BamFileHandler::openBam_24chr(bam_input, bai_input);
		string file_alupos = file_alupos_prefix + "alu_chr0";

		delete_search(bam_fh, file_fa_prefix, chrns, fn_tmp1, fn_log1, file_alupos, rg_to_idx);
		delete bam_fh;


		move_files(path0+"log1s/", path0 + pn + ".log1") ;
		move_files(path0+"tmp1s/", path0 + pn + ".tmp1") ;

		// now try to calculate prob of genotypes 
		string fn_tmp2 = path0 + "tmp2s/" + pn + ".tmp2";
		map <int, EmpiricalPdf *> pdf_rg;    
		map<string, int> isPairEnd;
		string file_seq_info = cf_fh.get_conf("file_seq_info", false);
		if (!file_seq_info.empty())
			get_RG_info(file_seq_info, pn, isPairEnd);
		read_pdf_pn(file_dist_prefix, pn, ::p_alu.pdf_param, pdf_rg, isPairEnd);

		ofstream fout(fn_tmp2.c_str());
		fout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";  
		string line;
		fn_tmp1 = path0 + "tmp1s/" + pn + ".tmp1";
		ifstream fin(fn_tmp1.c_str());
		assert(fin);
		getline(fin, line); // read header
		while (getline(fin, line))
		{
			int flagInt = parseline_del(fout, line, pdf_rg, ::p_alu.log10_ratio_ub, 0);
			if ( flagInt > 0) parseline_del(fout, line, pdf_rg, ::p_alu.log10_ratio_ub, flagInt);
		}
		fin.close();
		fout.close();
		EmpiricalPdf::delete_map(pdf_rg);
	}
	else if (opt == "2")
	{   
		// write vcf for all pn
		// step 2: create a big vcf file for genotype calls of all samples
		// USAGE: opt3/alu_delete [configuration file] 2
		vector <string> pns;
		string pns_to_use = cf_fh.get_conf("pn_del_vcf", false);
		if ( pns_to_use.empty() ) 
			pns_to_use = cf_fh.get_conf("file_pn");

		GetUsedSampleName(pns_to_use, pns); 

		map <string, std::set<int> > chrn_aluBegin;
		for ( vector <string>::iterator pi = pns.begin(); pi != pns.end(); pi++) 
			read_highCov_region(path0 + "log1s/" + *pi + ".log1",chrn_aluBegin);

		string path_input = path0 + "tmp2s/";
		string fn_prefix = path0 + int_to_string( pns.size()) ;
		string ref_name = cf_fh.get_conf("ref_name");
		float llh_th = 1.92; //qchisq(0.95, df=1)  [1] 3.841459 
		combine_pns_vcf(path_input, ".tmp2", fn_prefix + ".vcf", pns, chrns, chrn_aluBegin, llh_th, ref_name);  

	}
	else if (opt == "debug1")
	{ 
		FastaFileHandler *fasta_fh = new FastaFileHandler(file_fa_prefix + "chr1.fa", "chr1");
		seqan::CharString fa_seq;
		fasta_fh->fetch_fasta_upper(33465, 33495, fa_seq);
		cout << fa_seq << endl;
		delete fasta_fh;

	}
	else if (opt == "debug2")
	{ // debugging and manually check some regions 

		string pn = "1006-01";
		string line, output_line;
		map <int, EmpiricalPdf *> pdf_rg;    
		map<string, int> isPairEnd;
		read_pdf_pn(file_dist_prefix, pn, ::p_alu.pdf_param, pdf_rg, isPairEnd);

		cout << "chr aluBegin aluEnd midCnt clipCnt unknowCnt 00 01 11\n";    
		//line = "chr21 21903875 21904163 26 48 5 2 0:525 0:534";
		//line = "chr19 22768154 22768470 487 35 4 4 0:513 0:540 2:713 2:743";
		//line = "chr1 119558781 119559091 674 82 3 1 3:740";
		line = "chr2 157077122 157077335 493 1 5 0";
		int flagInt = parseline_del(cout, line, pdf_rg, ::p_alu.log10_ratio_ub, 0, true);
		if ( flagInt > 0) parseline_del(cout, line, pdf_rg, ::p_alu.log10_ratio_ub, flagInt, true);

	}
	else
	{
		cout << "unknown options !\n";
	}

	delete io_path;
	cout << "time used " << clocki.elapsed() << endl;
	return 0;  
}

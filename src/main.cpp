#include <stdio.h>
#include<algorithm>
#include<fstream>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <omp.h>

#include "Tachyon/Tachyon.h"
#include "Utils/util.h"
#include "Utils/Bioinformatics/Writer.hpp"


using namespace std;


void printHelp(boost::program_options::options_description general,
               boost::program_options::options_description makedb,
               boost::program_options::options_description aligner);

OutputType strToOutputType(const std::string &str);

ScoreMatrixType strToScorerType(const std::string &str);

AlignmentType strToAlignmentType(const std::string &str);


int main(int argc, const char *argv[]) {

	std::string database_path;

	std::string reduced_database;

	string command_;
	Command command;
	namespace po = boost::program_options;
	po::options_description general("General options");

	long threads;
	long kmer_len;
	string kmol_high, kmol_low;

	long seg_window;
	string seg_low_cut_;
	string seg_high_cut_;

	general.add_options()
			("help,h", "produce help message")
			("threads,p", po::value<long>(&threads)->default_value(8), "number of CPU threads")
			("db,d", po::value<string>(&database_path), "path to original nr file")
			("in,i", po::value<string>(&reduced_database), "path to reduced database")
			("kmer-len,l", po::value<long>(&kmer_len)->default_value(5), "k-mer length\n");

	po::options_description makedb("Makedb options");

	makedb.add_options()
			("kmol-high", po::value<string>(&kmol_high)->default_value("5,0"),
			 "most common k-mer over which length, e.g \n"
					 "5,0 means 5 k-mers over full length\n"
					 "3,60 means 3 kmer for every 60 AA \n")
			("kmol-low", po::value<string>(&kmol_low)->default_value("7,0"),
			 "least common k-mer over which length, e.g \n"
					 "5,0 means 5 k-mers over full length\n"
					 "3,60 means 3 kmer for every 60 AA \n")
			("seg-window", po::value<long>(&seg_window)->default_value(12), "Seg window")
			("seg-low-cut", po::value<string>(&seg_low_cut_)->default_value("2.2"), "seg lowCut")
			("seg-high-cut", po::value<string>(&seg_high_cut_)->default_value("2.5"), "Seg highCut");

	po::options_description aligner("Aligner options");

	std::string queries_path;

	int gap_open;
	int gap_extend;

	string matrix;
	std::string out_path;
	std::string out_format_string;

	double max_evalue;
	int max_alignments;
	string algorithm;
	int high_match;
	int low_match;

	aligner.add_options()
			("query,q", po::value<string>(&queries_path), "input query file")
			("high-match", po::value<int>(&high_match)->default_value(3), "minimum number of common kmer match")
			("low-match", po::value<int>(&low_match)->default_value(2), "minimum number of least common kmer match")
			("gapopen,g", po::value<int>(&gap_open)->default_value(10), "gap open penalty, default=10")
			("gapext,e", po::value<int>(&gap_extend)->default_value(1), "gap extend penalty, default=1")
			("matrix,m", po::value<string>(&matrix)->default_value("BLOSUM_62"), "score matrix")
			("out,o", po::value<string>(&out_path), "output file")
			("out-format", po::value<string>(&out_format_string)->default_value("bm9"), "outfmt <string>\n"
					"out format must be one of the following:\n"
					"bm0 \t- blast m0 output format\n"
					"bm8 \t- blast m8 tabular output format\n"
					"bm9 \t- blast m9 commented tabular output format")
			("evalue,v", po::value<double>(&max_evalue)->default_value(0.001), "evalue, default=0.001")
			("max-aligns,a", po::value<int>(&max_alignments)->default_value(10), "max aligns, default=10")
			("algorithm,A", po::value<string>(&algorithm)->default_value("SW"), "algorithm used for alignment.\n"
					"Must be on of the following:\n"
					"SW - Smith-Waterman local alignment\n"
					"NW - Needleman-Wunsch global alignment\n"
					"HW - semiglobal alignment\n"
					"OV - overlap alignment\n");


	po::options_description hidden("Hidden options");
	hidden.add_options()("command", po::value<string>(&command_));

	po::options_description all("Command line options");

	all.add(general).add(hidden).add(makedb).add(aligner);

	po::positional_options_description positional;
	positional.add("command", -1);

	po::variables_map vm;

	try {
		po::store(
				po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
		po::notify(vm);

		if (vm.count("help")) {
			printHelp(general, makedb, aligner);
			exit(-1);
		}

		command = strToCommand(command_);

		ScoreMatrixType scorer_type = strToScorerType(matrix);
		ScoreMatrix scorer(scorer_type, gap_open, gap_extend);

		AlignmentType align_type = strToAlignmentType(algorithm);
		Type input_type;

		if (command == Command::makedb) {
			int high_numb, low_numb, high_olen, low_olen;

			size_t p = kmol_high.find(',');
			high_numb = stoi(kmol_high.substr(0, p));
			high_olen = stoi(kmol_high.substr(p + 1));

			p = kmol_low.find(',');
			low_numb = stoi(kmol_low.substr(0, p));
			low_olen = stoi(kmol_low.substr(p + 1));

			Base base(database_path.c_str(), reduced_database.c_str(), kmer_len, high_numb, high_olen, low_numb,
			          low_olen, seg_window, stod(seg_low_cut_), stod(seg_high_cut_));
			base.read();
			base.makeIndexes();
		} else {

			ChainSet queries;
			readFastaFile(queries_path.c_str(), queries, 0);
			int query_size = queries.size();

			OutSet results;
			results.clear();
			results.resize(query_size);

			if (command == Command::blastp) input_type = Type::PROTEINS;
			else if (command == Command::blastx) input_type = Type::NUCLEOTIDES;

			Tachyon tachyon(database_path.c_str(), reduced_database.c_str());
			RunParams runParams(Type::PROTEINS, high_match, low_match, gap_open, gap_extend, max_evalue, max_alignments, scorer_type, align_type);
			tachyon.search(queries, runParams, results);

			OutputType out_format = strToOutputType(out_format_string);
			Writer writer(out_path, out_format, scorer);


			for (const auto &it: results) {
				writer.write_alignments(it, queries, tachyon.getBase());
			}
		}
	} catch (invalid_argument &e) {
		printf("%s\n", e.what());
		printHelp(general, makedb, aligner);
		exit(-1);
	}

}

OutputType strToOutputType(const std::string &str) {

	if (str.compare("bm0") == 0) {
		return OutputType::kBm0;
	} else if (str.compare("bm8") == 0) {
		return OutputType::kBm8;
	} else if (str.compare("bm9") == 0) {
		return OutputType::kBm9;
	} else {
		throw invalid_argument("Unrecognised output format!");
	}

}

ScoreMatrixType strToScorerType(const std::string &str) {

	if (str.compare("BLOSUM_45") == 0) {
		return ScoreMatrixType::kBlosum45;
	} else if (str.compare("BLOSUM_50") == 0) {
		return ScoreMatrixType::kBlosum50;
	} else if (str.compare("BLOSUM_62") == 0) {
		return ScoreMatrixType::kBlosum62;
	} else if (str.compare("BLOSUM_80") == 0) {
		return ScoreMatrixType::kBlosum80;
	} else if (str.compare("BLOSUM_90") == 0) {
		return ScoreMatrixType::kBlosum90;
	} else if (str.compare("PAM_30") == 0) {
		return ScoreMatrixType::kPam30;
	} else if (str.compare("PAM_70") == 0) {
		return ScoreMatrixType::kPam70;
	} else if (str.compare("PAM_250") == 0) {
		return ScoreMatrixType::kPam250;
	} else {
		throw invalid_argument("Unrecognised matrix!");
	}
}

AlignmentType strToAlignmentType(const std::string &str) {

	if (str.compare("NW") == 0) {
		return AlignmentType::kNW;
	} else if (str.compare("HW") == 0) {
		return AlignmentType::kHW;
	} else if (str.compare("OV") == 0) {
		return AlignmentType::kOV;
	} else if (str.compare("SW") == 0) {
		return AlignmentType::kSW;
	} else {
		throw invalid_argument("Unrecognised alignment algorithm!");
	}

}

void printHelp(boost::program_options::options_description general,
               boost::program_options::options_description makedb,
               boost::program_options::options_description aligner) {
	cout << endl << "Syntax:" << endl;
	cout << "  tachyon COMMAND [OPTIONS]" << endl << endl;
	cout << "Commands:" << endl;
	cout << "  makedb\tCreate indexes from FASTA file" << endl;
	cout << "  blastp\tAlign amino acid query sequences against a protein reference database" << endl;
	cout << "  blastx\tAlign DNA query sequences against a protein reference database" << endl;
	cout << endl;
	cout << general << endl << makedb << endl << aligner << endl;
}

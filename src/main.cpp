#include <stdio.h>
#include<algorithm>
#include<fstream>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <omp.h>

#include "tachyon.h"
#include "util.h"
#include "writer.hpp"


using namespace std;


void printHelp(boost::program_options::options_description general,
               boost::program_options::options_description makedb,
               boost::program_options::options_description aligner);

OutputType strToOutputType(const std::string& str);

ScoreMatrixType strToScorerType(const std::string& str);

AlignmentType strToAlignmentType(const std::string& str);
Command strToCommand(const std::string& str);

int main(int argc, const char *argv[]) {

	std::string database_path;

	std::string reduced_database;

	string command_;
	Command command;
	namespace po = boost::program_options;
	po::options_description general("General options");

	long threads;
	general.add_options()
			("help,h", "produce help message")
			("threads,p", po::value<long>(&threads)->default_value(8), "number of CPU threads")
			("db,d", po::value<string>(&database_path), "path to original nr file");

	po::options_description makedb("Makedb options");
	makedb.add_options()
			("red,r", po::value<string>(&reduced_database), "path to reduced database\n");

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
	aligner.add_options()
			("in,i", po::value<string>(&reduced_database), "path to reduced database")
			("query,q", po::value<string>(&queries_path), "input query file")
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

	try{
		po::store(
				po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
		po::notify(vm);


		if (vm.count("help")) {
			printHelp(general,makedb,aligner);
			exit(-1);
		}

		command = strToCommand(command_);

		ScoreMatrixType scorer_type = strToScorerType(matrix);
		ScoreMatrix scorer(scorer_type,gap_open, gap_extend);

		AlignmentType  align_type = strToAlignmentType(algorithm);
		Type input_type;

		omp_set_dynamic(0);
		omp_set_num_threads(threads);

		Base base(database_path.c_str(), reduced_database.c_str());
		if(command == Command::makedb) {
			base.read();
			base.make_indexes();
		}else {
			base.dump_in_memory();
			EValue evalue_params(base.database_size(), scorer);
			ChainSet queries;
			readFastaFile(queries_path.c_str(), queries, 0);
			int query_size = queries.size();

			OutSet results;
			results.clear();
			results.resize(query_size);

			if(command == Command::blastp) input_type = Type::PROTEINS;
			else if(command == Command::blastx) input_type = Type::NUCLEOTIDES;

			Tachyon tachyon(base);
			tachyon.search(queries, input_type, results, align_type, max_evalue, max_alignments, evalue_params, scorer);

			OutputType  out_format = strToOutputType(out_format_string);
			Writer writer (out_path, out_format, scorer);


			for (const auto &it: results) {
				writer.write_alignments(it, queries, base);
			}
		}
	}catch(invalid_argument& e) {
		printf("%s\n",e.what());
		printHelp(general,makedb,aligner);
		exit(-1);
	}

}

OutputType strToOutputType(const std::string& str) {

	if (str.compare("bm0") == 0) {
		return OutputType::kBm0;
	} else if (str.compare("bm8") == 0) {
		return OutputType::kBm8;
	} else if (str.compare("bm9") == 0) {
		return OutputType::kBm9;
	}else{
		throw invalid_argument("Unrecognised output format!");
	}

}

ScoreMatrixType strToScorerType(const std::string& str) {

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
	}else {
		throw invalid_argument("Unrecognised matrix!");
	}
}

AlignmentType strToAlignmentType(const std::string& str) {

	if (str.compare("NW") == 0) {
		return AlignmentType::kNW;
	} else if (str.compare("HW") == 0) {
		return AlignmentType::kHW;
	} else if (str.compare("OV") == 0) {
		return AlignmentType::kOV;
	} else if (str.compare("SW") == 0) {
		return AlignmentType::kSW;
	}else {
		throw invalid_argument("Unrecognised alignment algorithm!");
	}

}

Command strToCommand(const std::string& str) {
	if(str == "makedb") return Command::makedb;
	else if(str == "blastp") return Command::blastp;
	else if(str == "blastx") return Command::blastx;
	else{
		throw invalid_argument("Unrecognised command!");
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
	cout<<general<<endl<<makedb<<aligner<<endl;
}
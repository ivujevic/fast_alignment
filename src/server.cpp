#include <stdio.h>
#include<algorithm>
#include<fstream>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <omp.h>
#include <boost/asio/io_service.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "tachyon.h"
#include "util.h"
#include "writer.hpp"
#include "util/json.hpp"

using namespace std;
using boost::asio::ip::tcp;
using json = nlohmann::json;

void printHelp(boost::program_options::options_description general,
               boost::program_options::options_description makedb);

OutputType strToOutputType(const std::string &str);

ScoreMatrixType strToScorerType(const std::string &str);

AlignmentType strToAlignmentType(const std::string &str);


typedef boost::shared_ptr<tcp::socket> socket_ptr;

void alg(string message, Base& base_) {
	json params = json::parse(message);

	string queries_path = params["q"];
	string type = params["type"];
	int high_match = params["hm"];
	int low_match = params["lm"];
	int gap_open = params["g"];
	int gap_extend = params["e"];
	string matrix = params["m"];
	string out_path = params["o"];
	string out_format_string = params["of"];
	double max_evalue = params["v"];
	int max_alignments = params["a"];
	string algorithm = params["A"];

	ScoreMatrixType scorer_type = strToScorerType(matrix);
	ScoreMatrix scorer(scorer_type, gap_open, gap_extend);
	AlignmentType align_type = strToAlignmentType(algorithm);

	Type input_type;

	if (type == "blastp") input_type = Type::PROTEINS;
	else if (type == "blastx") input_type = Type::NUCLEOTIDES;

	ChainSet queries;
	readFastaFile(queries_path.c_str(), queries, 0);

	int query_size = queries.size();

	OutSet results;
	results.clear();
	results.resize(query_size);

	EValue evalue_params(base_.database_size(), scorer);

	Tachyon tachyon(base_, high_match, low_match, base_.kmer_len());
	tachyon.search(queries, input_type, results, align_type, max_evalue, max_alignments, evalue_params, scorer);

	OutputType out_format = strToOutputType(out_format_string);
	Writer writer(out_path, out_format, scorer);

	for (const auto &it: results) {
		writer.write_alignments(it, queries, base_);
	}
}

void syn_session(socket_ptr socket, Base& base)
{
	boost::array<char, 1024> buffer;

	boost::system::error_code error;

	size_t len = socket->read_some(boost::asio::buffer(buffer));

	string message;
	copy(buffer.begin(), buffer.begin() + len, back_inserter(message));

	alg(message, base);
	boost::system::error_code ignored_error;
	boost::asio::write(*socket.get(), boost::asio::buffer("Finished"));
	cout<<"Finish"<<endl;
}

void syn_server(boost::asio::io_service& io_service, short port, Base& base) {
	tcp::acceptor a(io_service, tcp::endpoint(tcp::v4(), port));
	for (;;)
	{
		socket_ptr sock(new tcp::socket(io_service));
		a.accept(*sock);
		boost::thread t(boost::bind(syn_session, sock, base));
	}
}

int main(int argc, const char *argv[]) {

	std::string database_path;

	std::string reduced_database;

	string command_;
	Command command;
	namespace po = boost::program_options;
	po::options_description general("Server options");

	long threads;
	unsigned short port = 9341;
	long kmer_len;

	string kmol_high, kmol_low;

	long seg_window;
	string seg_low_cut_;
	string seg_high_cut_;

	general.add_options()
			("help,h", "produce help message")
			("threads,n", po::value<long>(&threads)->default_value(8), "max number of CPU threads per query")
			("in,i", po::value<string>(&reduced_database), "path to reduced database")
			("db,d", po::value<string>(&database_path), "path to original nr file")
			("port,p", po::value<unsigned short>(&port)->default_value(9341), "port on which server listening")
			("kmer-len,l", po::value<long>(&kmer_len)->default_value(5), "k-mer length\n");;

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


	po::variables_map vm;

	po::options_description hidden("Hidden options");
	hidden.add_options()("command", po::value<string>(&command_));

	po::positional_options_description positional;
	positional.add("command", -1);

	po::options_description all("Command line options");
	all.add(general).add(hidden).add(makedb);


	try {
		po::store(po::command_line_parser(argc, argv).options(all).positional(positional).run(), vm);
		po::notify(vm);

		if (vm.count("help")) {
			printHelp(general, makedb);
			exit(-1);
		}
		command = strToCommand(command_);

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
			base.make_indexes();
		} else {
			omp_set_dynamic(0);
			omp_set_num_threads(threads);

			Base base(database_path.c_str(), reduced_database.c_str());
			base.dump_in_memory();

			if (kmer_len != base.kmer_len()) {
				printf("Error: This database was reduced with different kmer length!\n");
				exit(-1);
			}

			boost::asio::io_service io_service;

			tcp::endpoint endpoint(tcp::v4(), port);

			syn_server(io_service,port, base);
		}

	} catch (invalid_argument &e) {
		printf("%s\n", e.what());
		printHelp(general, makedb);
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


void printHelp(boost::program_options::options_description general, boost::program_options::options_description makedb) {
	cout << endl << "Syntax:" << endl;
	cout << "  tachyon COMMAND [OPTIONS]" << endl << endl;
	cout << "Commands:" << endl;
	cout << "  makedb\tCreate indexes from FASTA file" << endl;
	cout << "  server\tStart server" << endl;
	cout << endl;
	cout << general << endl << makedb <<endl;
}
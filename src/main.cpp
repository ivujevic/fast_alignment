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
#include <boost/enable_shared_from_this.hpp>
#include "tachyon.h"
#include "util.h"
#include "writer.hpp"
#include "util/json.hpp"

using namespace std;
using boost::asio::ip::tcp;
using json = nlohmann::json;

void printHelp(boost::program_options::options_description general);

OutputType strToOutputType(const std::string &str);

ScoreMatrixType strToScorerType(const std::string &str);

AlignmentType strToAlignmentType(const std::string &str);

Command strToCommand(const std::string &str);

void serverThread(string &settings, tcp::socket* so,int i);


class session
		:public boost::enable_shared_from_this<session>
{
public:
	session(boost::asio::io_service& io_service, Base& base)
			:socket_(io_service),
			 base_(base){}

	tcp::socket& socket(){
		return socket_;
	}

	void start() {

		boost::array<char, 1024> buffer;

		boost::system::error_code error;
		size_t len = socket_.read_some(boost::asio::buffer(buffer));

		string message;
		copy(buffer.begin(), buffer.begin() + len, back_inserter(message));
		cout<<"Primio " <<message<<endl;
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

		Tachyon tachyon(this->base_, high_match, low_match, this->base_.kmer_len());
		tachyon.search(queries, input_type, results, align_type, max_evalue, max_alignments, evalue_params, scorer);

		OutputType out_format = strToOutputType(out_format_string);
		Writer writer(out_path, out_format, scorer);

		for (const auto &it: results) {
			writer.write_alignments(it, queries, base_);
		}
		boost::asio::write(socket_, boost::asio::buffer("Finished"));
	}

private:
	tcp::socket socket_;
	Base& base_;
};

typedef boost::shared_ptr<session> session_ptr;


class server
{
public:
	server(boost::asio::io_service& io_service, const tcp::endpoint& endpoint, Base& base)
			:io_service_(io_service),
			 acceptor_(io_service, endpoint),
	         base_(base)
	{
		cout<<"Server is started and ready for use\n";

		start_accpet();
	}

	void start_accpet() {
		session_ptr new_session(new session(io_service_, base_));
		Tachyon tachyon(base_, base_.kmer_len(), 5, 6);
		acceptor_.async_accept(new_session->socket(),
		                       boost::bind(&server::handle_accept, this, new_session, boost::asio::placeholders::error()));
	}

	void handle_accept(session_ptr session,
	                   const boost::system::error_code& error)
	{
		if (!error) {
			session->start();
		}

		start_accpet();
	}

private:
	boost::asio::io_service& io_service_;
	tcp::acceptor acceptor_;
	Base& base_;
};

typedef boost::shared_ptr<server> server_ptr;

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

	general.add_options()
			("help,h", "produce help message")
			("threads,n", po::value<long>(&threads)->default_value(8), "max number of CPU threads per query")
			("in,i", po::value<string>(&reduced_database), "path to reduced database")
			("db,d", po::value<string>(&database_path), "path to original nr file")
			("port,p", po::value<unsigned short>(&port)->default_value(9341), "port on which server listening")
			("kmer-len,l", po::value<long>(&kmer_len)->default_value(5), "k-mer length\n");;

	po::variables_map vm;

	try {
		po::store(po::command_line_parser(argc, argv).options(general).run(), vm);
		po::notify(vm);
		if (vm.count("help")) {
			printHelp(general);
			exit(-1);
		}
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

		server_ptr server_(new server(io_service, endpoint, base));

		io_service.run();
	} catch (invalid_argument &e) {
		printf("%s\n", e.what());
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

Command strToCommand(const std::string &str) {
	if (str == "makedb") return Command::makedb;
	else if (str == "blastp") return Command::blastp;
	else if (str == "blastx") return Command::blastx;
	else {
		throw invalid_argument("Unrecognised command!");
	}
}

void printHelp(boost::program_options::options_description general) {
	cout << endl << "Syntax:" << endl;
	cout << "  tachyon [OPTIONS]" << endl << endl;
	cout << general << endl;
}
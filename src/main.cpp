#include <stdio.h>
#include<algorithm>
#include<fstream>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <omp.h>
#include <boost/asio/io_service.hpp>
#include <boost/asio.hpp>

#include "tachyon.h"
#include "util.h"
#include "writer.hpp"


using namespace std;
using boost::asio::ip::tcp;

void printHelp(boost::program_options::options_description general,
               boost::program_options::options_description makedb,
               boost::program_options::options_description aligner);

OutputType strToOutputType(const std::string &str);

ScoreMatrixType strToScorerType(const std::string &str);

AlignmentType strToAlignmentType(const std::string &str);

Command strToCommand(const std::string &str);

void serverThread(string &settings, tcp::socket* so,int i);

int main(int argc, const char *argv[]) {

	std::string database_path;

	std::string reduced_database;

	string command_;
	Command command;
	namespace po = boost::program_options;
	po::options_description general("General options");

	long threads;
	long max_queries;
	unsigned short port = 9341;
	general.add_options()
			("help,h", "produce help message")
			("threads,n", po::value<long>(&threads)->default_value(8), "max number of CPU threads per query")
			("max-query", po::value<long>(&max_queries)->default_value(100), "max number of accepted queries")
			("in,i", po::value<string>(&reduced_database), "path to reduced database")
			("port,p", po::value<unsigned short>(&port)->default_value(9341), "port on which server listening")
			("db,d", po::value<string>(&database_path), "path to original nr file");

	cout << "Listening on " << port << endl;
	boost::asio::io_service io_service;

	tcp::acceptor acceptor(io_service, tcp::endpoint(tcp::v4(), port));


	#pragma omp parallel
	{
		#pragma omp single
		{
			for (int i =0; ;i++) {
				tcp::socket socket(io_service);

				acceptor.accept(socket);
				boost::array<char, 512> buffer;

				boost::system::error_code error;
				size_t len = socket.read_some(boost::asio::buffer(buffer));

				string message;
				copy(buffer.begin(), buffer.begin() + len, back_inserter(message));
				#pragma omp task
				{
					serverThread(message, &socket,i);
				}
				break;
			}
		}
	}


}

void serverThread(string &settings, tcp::socket* socket,int i) {
	//boost::asio::write(socket, boost::asio::buffer("Finished " + to_string(i)));

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
	cout << general << endl << makedb << aligner << endl;
}
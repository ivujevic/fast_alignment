#include <stdio.h>

#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/asio/io_service.hpp>
#include <boost/asio.hpp>
#include "util.h"

using namespace std;


using boost::asio::ip::tcp;

int main() {

	string port = "9341";

	std::string database_path;

	std::string reduced_database;

	string command_;
	Command command;
	namespace po = boost::program_options;
	po::options_description aligner("General options");

	string queries_path,matrix,out_path,out_format_string,algorithm;
	int gap_open,gap_extend,max_alignments;
	double max_evalue;

	aligner.add_options()
			("port,p", po::value<string>(&port), "port on which server listening")
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


	boost::asio::io_service io_service;

	tcp::resolver resolver(io_service);

	tcp::resolver::query query(ip_adress, port.c_str());
	tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);

	tcp::socket socket(io_service);
	boost::asio::connect(socket, endpoint_iterator);

	boost::array<char, 128> buffer;
	boost::system::error_code error;
	boost::system::error_code ignored_error;

	//send options to server

	string message = "Poslao sam poruku ! ";
	/*string message = queries_path + " " + to_string(gap_open) + to_string(gap_extend)
	                 + " " + matrix + " "  + out_path + " " + out_format_string + " "
	                 + to_string(max_evalue) + " " + to_string(max_alignments) + " " + algorithm;
*/
	cout<<"Saljem poruku"<<endl;
	boost::asio::write(socket, boost::asio::buffer(message), ignored_error);

	cout << "Query processing!" << endl;
	string received_message = "";
	while(message != "Finished") {
		try{
		size_t len = socket.read_some(boost::asio::buffer(buffer));
			cout.write(buffer.data(),len);
			break;
		}catch(...) {

		}
	}
	cout<<"Finised!"<<endl;
}
#include <cstdlib>
#include <deque>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/asio.hpp>
using boost::asio::ip::tcp;

#include "util/json.hpp"

using namespace std;
using json = nlohmann::json;

void printHelp(boost::program_options::options_description general);


int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

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
    string server_port;

    string type;

    aligner.add_options()
            ("help,h", "produce help message")
            ("query,q", po::value<string>(&queries_path), "input query file")
            ("server-port,p", po::value<string>(&server_port)->default_value("9341"), "port of server")
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

    po::variables_map vm;

    try {
        po::store(
                po::command_line_parser(argc, argv).options(aligner).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            printHelp(aligner);
            exit(-1);
        }

        json params;
        params["q"] = queries_path;
        params["type"] = type ;
        params["hm"] = high_match;
        params["lm"] = low_match;
        params["g"] = gap_open;
        params["e"] = gap_extend;
        params["m"] = matrix;
        params["o"] = out_path;
        params["of"] = out_format_string;
        params["v"] = max_evalue;
        params["a"] = max_alignments;
        params["A"] = algorithm;

        boost::asio::io_service io_service;

        tcp::resolver resolver(io_service);

        tcp::resolver::query query("localhost", server_port.c_str());
        tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);

        tcp::socket socket(io_service);
        boost::asio::connect(socket, endpoint_iterator);

        boost::array<char, 128> buffer;
        boost::system::error_code error;
        boost::system::error_code ignored_error;



        boost::asio::write(socket, boost::asio::buffer(params.dump()), ignored_error);

        cout << "Query processing!" << endl;
        string message;
        while(message !="Finished") {
            size_t len = socket.read_some(boost::asio::buffer(buffer));
            string message;
            copy(buffer.begin(),buffer.begin() + len,back_inserter(message));
            cout<<message<<endl;
        }
        cout<<"Finished"<<endl;
    } catch (invalid_argument &e) {
        printf("%s\n", e.what());
        exit(-1);
    }
    return 0;
}

void printHelp(boost::program_options::options_description general) {
    cout << endl << "Syntax:" << endl;
    cout << " runquery [OPTIONS]" << endl << endl;
    cout << general << endl;
}
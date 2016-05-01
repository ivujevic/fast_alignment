#include <iostream>
#include <fstream>
#include <string>
#include <boost/asio.hpp>
#include "Constants.h"
#include<boost/array.hpp>

using boost::asio::ip::tcp;
using namespace std;

void printHelp() {
    cout << "\nUsage: runquery [-q query] [-c Score_Cutoff] [-oDir output_directory] [-maxDepth max_depth]\n\
                 -q   \n\
                       filename containing sequence(s) in fasta format.\n\
                 -c \n\
                       OPTIONAL  Only hits with scores below the chosen value will be displayed.\n\
                                          (default 0.015)\n\
                 -oDir\n \
                       path to directory for output\n\
                 -maxDepth \n\
                       Max depth of tachyon iterative. First step of iteration has depth 1 \n\
                 \n\n";
}

int main(int argc, char * argv[]) {
  string queryFile = "";
  string cut_off = "0.015";
  string out_dir = "";
  string out_format = "-d";
  string maxDepth ="0";

  int nec = 0;

  if (argc < 7) {
    printHelp();
    return 0;
  } else {
    for (int i = 1; i < argc; i++) {
      string argument = argv[i];

      if (argument == "-q") {
        if (i + 1 < argc) {
          queryFile = argv[i + 1];
          nec++;
        } else {
          printHelp();
          return 0;
        }
      } else if (argument == "-c") {
        if (i + 1 < argc) {
          cut_off = argv[i + 1];
        } else {
          printHelp();
          return 0;
        }
      } else if (argument == "-oDir") {
        if (i + 1 < argc) {
          out_dir = argv[i + 1];
          nec++;
        } else {
          printHelp();
          return 0;
        }
      }else if( argument == "-maxDepth"){
        if( i+ 1 < argc){
          maxDepth = argv[i+1];
          nec++;
        }else{
          printHelp();
          return 0;	
        }
      }
    }
  }
  if( nec != 3) {
    printHelp();
    return 0;
  }
  try {
    boost::asio::io_service io_service;

    tcp::resolver resolver(io_service);

    tcp::resolver::query query(ip_adress, to_string(server_port).c_str());
    tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);

    tcp::socket socket(io_service);
    boost::asio::connect(socket, endpoint_iterator);


    boost::array<char, 128> buffer;
    boost::system::error_code error;
    boost::system::error_code ignored_error;

    //send options to server


    string message = queryFile + ";;" + cut_off + ";;" + out_dir + ";;" + maxDepth; //use ;; for separate
    boost::asio::write(socket, boost::asio::buffer(message), ignored_error);

    cout << "Query processing!" << endl;

    while(message !="Finished") {
      size_t len = socket.read_some(boost::asio::buffer(buffer));
      string message;
      copy(buffer.begin(),buffer.begin() + len,back_inserter(message));
      cout<<message<<endl;
    }
    cout<<"Finished"<<endl;
  }catch(std::exception& e) {
  }

  return 0;
}

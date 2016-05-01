/* 
 * File:   DatabaseElement.h
 * Author: vujevic
 *
 * Created on April 25, 2014, 11:43 PM
 */
#pragma once
#include<iostream>
#include<string>
#include<boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include<boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace std;


/**
 * Class which represents one element from database in fasta format.
 */
class DatabaseElement {
public:


	DatabaseElement(uint32_t id, char* name, uint32_t name_length,
	                char* data, uint32_t data_length);

	DatabaseElement(string name, int name_len, string sequence, int seq_len,int id)
			:name_(name),name_len_(name_len),sequence_(sequence),sequence_len_(seq_len),id_(id){}
    /**
     * Getters and setters for private members.
     */
    
    long id() const{return id_;}
    const string& name() const{return name_;}
	int name_len() const{return name_len_;}

	const string& sequence() const{return sequence_;}
	int sequence_len() const{return sequence_len_;};
		
private:

    string name_;
	int name_len_;
	string sequence_;
	int sequence_len_;

    long id_;
    
    
};

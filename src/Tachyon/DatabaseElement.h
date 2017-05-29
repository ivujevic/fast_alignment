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


/**
 * Class which represents one element from a database in fasta format.
 */
class DatabaseElement {
public:


	DatabaseElement();
	DatabaseElement(int id, char* name, int nameLen,
	                char* sequence, int sequenceLen);


    /**
     * Getters and setters for private members.
     */

    long id() const;
    const std::string& getName() const;
	int getNameLen() const;

	const std::string& getSequence() const;
	int getSequenceLen() const;

private:
	long id_;
	std::string sequence_;
	int sequenceLen_;
    std::string name_;
	int nameLen_;
};

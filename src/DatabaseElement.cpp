/* 
 * File:   DatabaseElement.cpp
 * Author: vujevic
 * 
 */

#include "DatabaseElement.h"

std::vector<char> kCoder = {
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,   0,   1,   2,   3,   4,
        5,   6,   7,   8,   9,  10,  11,  12,  13,  14,
        15,  16,  17,  18,  19,  20,  21,  22,  23,  24,
        25,  -1,  -1,  -1,  -1,  -1,  -1,   0,   1,   2,
        3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
        13,  14,  15,  16,  17,  18,  19,  20,  21,  22,
        23,  24,  25,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
        -1,  -1,  -1,  -1,  -1
};

DatabaseElement::DatabaseElement(uint32_t id, char *name, uint32_t name_length, char *data, uint32_t data_length) {
    // remove trailing white spaces
    while (isspace(name[name_length - 1])) {
        --name_length;
    }


    std::string data_;
    data_.reserve(data_length);
    uint32_t data_ptr = 0, valid_data_length = 0;

    for (uint32_t i = 0; i < data_length; ++i) {
        auto c = kCoder[data[i]];
        if (c != -1) {
            data[data_ptr++] = c +'A';
            ++valid_data_length;
        }
    }

    id_ = id;
    name_ = std::string(name,name_length);
    name_len_ = name_length;
    sequence_ = std::string(data,valid_data_length);
    sequence_len_ = valid_data_length;
}
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

DatabaseElement::DatabaseElement() {

}

DatabaseElement::DatabaseElement(int id, char* name, int nameLen, char* sequence,
                                 int sequenceLen) {

    // remove trailing white spaces
    while (isspace(name[nameLen - 1])) {
        nameLen--;
    }


    std::string data_;
    data_.reserve(sequenceLen);
    uint32_t data_ptr = 0, valid_data_length = 0;

    for (uint32_t i = 0; i < sequenceLen; ++i) {
        auto c = kCoder[sequence[i]];
        if (c != -1) {
            sequence[data_ptr++] = c +'A';
            ++valid_data_length;
        }
    }

    id_ = id;
    name_ = std::string(name, 0, nameLen);
    nameLen_ = nameLen;
    sequence_ = std::string(sequence, 0, valid_data_length);
    sequenceLen_ = valid_data_length;

}

const std::string& DatabaseElement::getName() const { return name_; }
int DatabaseElement::getNameLen() const { return nameLen_; }
const std::string& DatabaseElement::getSequence() const { return sequence_; }
int DatabaseElement::getSequenceLen() const { return sequenceLen_; }
long DatabaseElement::id() const { return id_; }

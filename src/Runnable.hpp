//
// Created by vujevic on 4/24/17.
//

#pragma once

class Runnable {
public:
    Runnable() {}

    virtual ~Runnable() {}
    virtual  void run() = 0;
};

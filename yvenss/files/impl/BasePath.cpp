/*
 * BasePath.cpp
 *
 *  Created on: Feb 25, 2017
 *      Author: ary
 */

#include "../headers/BasePath.h"

#include <stdio.h>
#include <errno.h>

BasePath::BasePath() {
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))) {
        throw errno;
    }

    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0';
    string str(cCurrentPath);
    path = str;
}

string BasePath::getBasePath() {
    return path;
}

BasePath::~BasePath() {
}


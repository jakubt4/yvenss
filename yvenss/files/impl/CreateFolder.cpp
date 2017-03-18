/*
 * CreateFolder.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#include "../headers/CreateFolder.h"

CreateFolder::CreateFolder(string path) {
    struct stat sb;
    if (stat(path.c_str(), &sb) != 0 && !S_ISDIR(sb.st_mode)) {
        mkdir(path.c_str(), 0777);
        cout << "Creation of " << path << endl;
    }
}

CreateFolder::~CreateFolder() {
}

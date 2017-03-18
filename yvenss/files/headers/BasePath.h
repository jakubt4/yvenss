/*
 * BasePath.h
 *
 *  Created on: Feb 25, 2017
 *      Author: ary
 */

#ifndef EVENT_SHAPE_SORT_UTILS_FILES_BASEPATH_H_
#define EVENT_SHAPE_SORT_UTILS_FILES_BASEPATH_H_

#include <string>
#include <stdio.h>
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

class BasePath {
    private:
        string path;
    public:
        BasePath();
        string getBasePath();
        virtual ~BasePath();
};

#endif /* EVENT_SHAPE_SORT_UTILS_FILES_BASEPATH_H_ */

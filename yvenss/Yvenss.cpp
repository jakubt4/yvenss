//============================================================================
// Name        : Yvenss.cpp
// Author      : Jakub Toth
// Version     :
// Copyright   : open source
// Description : Sorting of events
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "generator/Generator.h"
#include "sorter/Sorter.h"

using namespace std;

/*
 * Main method with parameter for generating specific number of events
 */
int main(int argc, char* argv[]) {
    bool rotation = true;
    string path = "";
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-gen") == 0 || strcmp(argv[i], "--generate") == 0) {
            if ((i + 1) < argc) {
                Generator eveGen(atoi(argv[i + 1]));
                i++;
                eveGen.generate();
            } else {
                Generator eveGen(1000);
                eveGen.generate();
            }
        } else if (strcmp(argv[i], "--without-rotations") == 0 || strcmp(argv[i], "-wr") == 0) {
            rotation = false;
        } else if (strcmp(argv[i], "--path") == 0 || strcmp(argv[i], "-p") == 0) {
            if ((i + 1) < argc) {
                path = string(argv[i + 1]);
                i++;
            } else {
                cout << "WARNING >> Path wasn't entered. Loading from /generated-files/." << endl;
            }
        }
    }
    Sorter s = Sorter(rotation, path);
    s.sort();
    cout << "SORTING DONE" << endl;

	return EXIT_SUCCESS;
}

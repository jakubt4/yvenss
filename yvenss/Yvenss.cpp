//============================================================================
// Name        : Yvenss.cpp
// Author      : Jakub Toth
// Version     :
// Copyright   : open source
// Description : main class of Yvenss
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "generator/Generator.h"
#include "sorter/Sorter.h"

using namespace std;

int main(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-gen") == 0) {
            if ((i + 1) < argc) {
                Generator eveGen(atoi(argv[i + 1]));
                i++;
                eveGen.generate();
            } else {
                Generator eveGen(1000);
                eveGen.generate();
            }
        } else {
            cout << "Bad parameter !!" << endl;
            cout << "Try: -gen int" << endl;
            return 0;
        }
    }

    Sorter s;
    s.sort();
    cout << "SORTING DONE" << endl;

	return EXIT_SUCCESS;
}

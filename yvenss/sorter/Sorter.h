/*
 * Sorter.h
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#ifndef SORTER_H_
#define SORTER_H_

#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <iomanip>

#include "../files/headers/PathConstants.h"
#include "../files/headers/BasePath.h"
#include "../files/headers/CreateFolder.h"

    class Sorter {
    private:
        int multiplicitySorterFields = 3;
        int angleBins = 20;
        int eventFields = angleBins + multiplicitySorterFields;
        int multiplicityField = 20;
        int sorterField = 21;
        int sumAngleBinsByIFields = 21;
        int myBinId = 22;
        int bins = 10;
        long double PI = 3.14159265358979323846;
        long double baseAngle = PI / 10.0;
        int eventsPerBin;
        int events;
        long double ONE_DEGREE = PI / 180;

        typedef vector<long double> Row; // One row of the matrix
        typedef vector<Row> Events; //  a vector of rows
        typedef vector<Events> Bins;

        Bins binsMatrixOrig;
    public:
        Sorter();
        void sort();
        void printAngleBins();
        void printEvents();
        virtual ~Sorter();
    };

#endif /* SORTER_H_ */

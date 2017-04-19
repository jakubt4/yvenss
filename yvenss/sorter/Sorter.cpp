/*
 * Sorter.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#include "Sorter.h"

/*
 * Comparator method for comparison of two objects (used like parameter for sort() method) under an array (vector).
 */
bool compare(vector<long double> row_a, vector<long double> row_b) {
    return row_a[21] > row_b[21];
}

/*
 * Constructor of sorter prepare data (read them a make initial sorting) for sorting.
 */
Sorter::Sorter() {

    //remove old results
    BasePath bp;
    string removeOld = "rm -r " + bp.getBasePath() + RESULT;
    system(removeOld.c_str());

    cout << "--------------------------------------------------" << endl;
    cout << "Loading data.." << endl;

    //path to data for loading
    string path = bp.getBasePath() + GEN_EVENTS_FILE_PATH;
    std::ifstream ifile_ev(path.c_str(), std::ios::in);

    //read number of events
    events = 0;
    int actualM = 0;
    std::vector<int> m;
    while (ifile_ev >> events) {
        ifile_ev >> actualM;
        m.push_back(actualM);
    }
    ifile_ev.close();

    //events per bin
    eventsPerBin = events / bins;
    cout << "EVENTS = " << events << "  ::  EVENTS PER BIN = " << eventsPerBin << endl;

    //2D vector for empirical distributions
    Events eventsMatrix(events, Row(eventFields));
    Events eventsMatrixTest(events, Row());
    //setup initial values for 2D vector
    for (int i = 0; i < events; i++) {
        for (int j = 0; j < eventFields; j++) {
            eventsMatrix[i][j] = 0.0;
            eventsMatrixTest[i].push_back(0.0);
        }
        eventsMatrix[i][multiplicityField] = m[i];
        eventsMatrixTest[i][multiplicityField] = m[i];
    }
    m.clear();

    int event = 1;
    int min = 1;
    int max = 1000;

    long double angle = 0.0;

    //Values for sorting by q2
    //long double q_x = 0.0;
    //long double q_y = 0.0;
    //long double q_2 = 0.0;

    //Sorting by q2
//    for (int i = 0; i < events; i++) {
//        if (events >= 1000) {
//            if (i % 1000 == 0) {
//                min = i + 1;
//                max = i + 1000;
//            }
//            path = bp.getBasePath() + GEN_PATH + "generated_particles_" + to_string(min) + "_" + to_string(max) + "/"
//                    + GEN_PARTICLES + to_string(i + 1);
//        } else {
//            path = bp.getBasePath() + GEN_PARTICLES_FILE_PATH + to_string(i + 1);
//        }
//        std::ifstream ifile(path.c_str(), std::ios::in);
//
//        while (ifile >> event) {
//            ifile >> angle;
//            q_x += cos(2.0 * angle);
//            q_y += sin(2.0 * angle);
//            for (int j = 0; j < 20; j++) {
//                if (angle >= (baseAngle * j) && angle < (baseAngle * (j + 1.0))) {
//                    eventsMatrix[i][j] += 1;
//                    break;
//                }
//            }
//        }
//
//        ifile.close();
//
//        q_2 = sqrt(pow(q_x, 2) + pow(q_y, 2)) / sqrt(eventsMatrix[i][multiplicityField]);
//        eventsMatrix[i][sorterField] = q_2;
//        q_x = 0;
//        q_y = 0;
//        if ((i + 1) % 1000 == 0) {
//            cout << event << '/' << events << endl;
//        }
//    }

    //read data from files and prepare empirical distributions
    for (int i = 0; i < events; i++) {
        //prepare path for file to read data from
        if (events >= 1000) {
            if (i % 1000 == 0) {
                min = i + 1;
                max = i + 1000;
            }
            path = bp.getBasePath() + GEN_PATH + "generated_particles_" + to_string(min) + "_" + to_string(max) + "/"
                    + GEN_PARTICLES + to_string(i + 1);
        } else {
            path = bp.getBasePath() + GEN_PARTICLES_FILE_PATH + to_string(i + 1);
        }

        //read data from file
        std::ifstream ifile(path.c_str(), std::ios::in);

        Row actualEvent;
        while (ifile >> event) {
            ifile >> angle;
            actualEvent.push_back(angle);
        }

        ifile.close();

        std::sort(actualEvent.begin(), actualEvent.end());
        //empirical distribution
        for (int k = 0; k < actualEvent.size(); k++) {
            for (int j = 0; j < 20; j++) {
                if (actualEvent[k] >= (baseAngle * j) && actualEvent[k] < (baseAngle * (j + 1.0))) {
                    eventsMatrixTest[i][j] += 1;
                    eventsMatrixTest[i].push_back(actualEvent[k]);
                    break;
                }
            }
        }

        if ((i + 1) % 1000 == 0) {
            cout << event << '/' << events << endl;
        }
    }

    cout << "Data loaded." << endl;
    cout << "--------------------------------------------------" << endl;

    //variance of events by empirical distributions
    cout << "Initial sorting." << endl;

    long double x[20];
    long double av_x[events];

    //prepare average from intervals for every interval
    for (int j = 0; j < 20; j++) {
        x[j] = (baseAngle * j) + (((baseAngle * (j + 1.0)) - (baseAngle * j)) / 2);
    }

    //compute normalized average of the interval classification
    for (int i = 0; i < events; i++) {
        long double sum_av_x = 0.0;
        for (int j = 0; j < 20; j++) {
//            sum_av_x += (eventsMatrix[i][j] * x[j]);
            sum_av_x += (eventsMatrixTest[i][j] * x[j]);
        }
//        av_x[i] = sum_av_x / eventsMatrix[i][multiplicityField];
        av_x[i] = sum_av_x / eventsMatrixTest[i][multiplicityField];
    }

    //compute variance
    Row vec_res(events);

    for (int i = 0; i < events; i++) {
        long double result_x = 0.0;
        for (int j = 0; j < 20; j++) {
//            result_x += (pow((x[j] - av_x[i]), 2) * eventsMatrix[i][j]);
            result_x += (pow((x[j] - av_x[i]), 2) * eventsMatrixTest[i][j]);
        }
//        result_x = result_x / eventsMatrix[i][multiplicityField];
        result_x = result_x / eventsMatrixTest[i][multiplicityField];
        vec_res[i] = result_x;
    }

    //set results of variance like sorting value
    std::sort(vec_res.begin(), vec_res.end());

    for (int i = 0; i < events; i++) {
//        eventsMatrix[i][sorterField] = vec_res[i];
        eventsMatrixTest[i][sorterField] = vec_res[i];
    }

//    std::sort(eventsMatrix.begin(), eventsMatrix.end(), compare);
    std::sort(eventsMatrixTest.begin(), eventsMatrixTest.end(), compare);

    Events newAngles(events, Row());

    for (int i = 0; i < events; i++) {
        for (int j = eventFields; j < eventsMatrixTest[i].size(); j++) {
            newAngles[i].push_back(eventsMatrixTest[i][j]);
        }
    }

    Row previousEvent;
    previousEvent.assign(newAngles[0].begin(), newAngles[0].end());
    long double kolmogorov_smirnovov_previousEv = (2 * PI) - previousEvent[previousEvent.size() - 1];
    for (int i = 0; i < (previousEvent.size() - 1); i++) {
        kolmogorov_smirnovov_previousEv += ((i + 1) * (previousEvent[i + 1] - previousEvent[i])) / previousEvent.size();
    }

    long double kolmogorov_smirnovov_actualEv = 0.0;

    for (int i = 1; i < events; i++) {
        kolmogorov_smirnovov_actualEv = (2 * PI) - newAngles[i][newAngles[i].size() - 1];
        for (int j = 0; j < (newAngles[i].size() - 1); j++) {
            kolmogorov_smirnovov_actualEv += ((j + 1) * (newAngles[i][j + 1] - newAngles[i][j])) / newAngles[i].size();
        }

        if (kolmogorov_smirnovov_actualEv == kolmogorov_smirnovov_previousEv) {
            break;
        }

        long double baseDiff = 0.0;
        if (kolmogorov_smirnovov_actualEv > kolmogorov_smirnovov_previousEv) {
            baseDiff = kolmogorov_smirnovov_actualEv - kolmogorov_smirnovov_previousEv;
            while (true) {
                for (int j = 0; j < (newAngles[i].size()); j++) {
                    newAngles[i][j] -= ONE_DEGREE;
                    if (newAngles[i][j] < 0) {
                        newAngles[i][j] += (2 * PI);
                    }
                }
                kolmogorov_smirnovov_actualEv = (2 * PI) - newAngles[i][newAngles[i].size() - 1];
                for (int j = 0; j < (newAngles[i].size() - 1); j++) {
                    kolmogorov_smirnovov_actualEv += ((j + 1) * (newAngles[i][j + 1] - newAngles[i][j]))
                            / newAngles[i].size();
                }
                long double actualDiff = kolmogorov_smirnovov_actualEv - kolmogorov_smirnovov_previousEv;
                if (baseDiff > actualDiff) {
                    baseDiff = actualDiff;
                } else {
                    break;
                }
            }
        } else {
            while (true) {
                for (int j = 0; j < (newAngles[i].size()); j++) {
                    newAngles[i][j] += ONE_DEGREE;
                    if (newAngles[i][j] > (2 * PI)) {
                        newAngles[i][j] -= (2 * PI);
                    }
                }
                kolmogorov_smirnovov_actualEv = (2 * PI) - newAngles[i][newAngles[i].size() - 1];
                for (int j = 0; j < (newAngles[i].size() - 1); j++) {
                    kolmogorov_smirnovov_actualEv += ((j + 1) * (newAngles[i][j + 1] - newAngles[i][j]))
                            / newAngles[i].size();
                }
                long double actualDiff = kolmogorov_smirnovov_actualEv - kolmogorov_smirnovov_previousEv;
                if (baseDiff > actualDiff) {
                    baseDiff = actualDiff;
                } else {
                    break;
                }
            }
        }
        kolmogorov_smirnovov_previousEv = kolmogorov_smirnovov_actualEv;
        previousEvent.clear();
        previousEvent.assign(newAngles[i].begin(), newAngles[i].end());
    }

    for (int i = 0; i < events; i++) {
        for (int k = 0; k < newAngles[i].size(); k++) {
            for (int j = 0; j < 20; j++) {
                if (newAngles[i][k] >= (baseAngle * j) && newAngles[i][k] < (baseAngle * (j + 1.0))) {
                    eventsMatrix[i][j] += 1;
                    break;
                }
            }
        }
    }

    //3D vector for empirical distributions of events in event bins - distribution of events from 2D vector to event
    //bins in 3D vector
    Bins binsMatrix(bins, Events(eventsPerBin, Row(eventFields)));

    int actualEvent = 0;
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < eventsPerBin; j++) {
            for (int k = 0; k < eventFields; k++) {
                binsMatrix[i][j][k] = eventsMatrix[actualEvent][k];
//                binsMatrix[i][j][k] = eventsMatrixTest[actualEvent][k];
                if (k == myBinId) {
                    binsMatrix[i][j][k] = i;
                }
            }
            actualEvent += 1;
        }
    }

    //assign 3D vector to global value of 3D vector
    binsMatrixOrig.assign(binsMatrix.begin(), binsMatrix.end());

    cout << "Initial sorting DONE." << endl;
    cout << "--------------------------------------------------" << endl;

    //clean up of local objects to release memory
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < eventsPerBin; j++) {
            binsMatrix[i][j].clear();
        }
        binsMatrix[i].clear();
    }
    binsMatrix.clear();
    for (int i = 0; i < events; i++) {
        eventsMatrix[i].clear();
    }
    eventsMatrix.clear();
}

/*
 * Sorting method for sorting empirical distributions.
 */
void Sorter::sort() {
    //sorting timer - START
    clock_t begin = clock();

    bool sorter = true;

    int cyc = 0;
    int diff = 0;
    int previousDiff = -1;
    int probablySorted = 0;
    int probablySortedMax = 10;

    //cycle for sorting
    while (sorter) {

        sorter = false;

        //number of cycles
        cyc++;

        //cycle timer - START
        clock_t begin_cyc = clock();

        /*STARTING OF ALGHORITM FOR SORTING*/
        /*********************************************************************************************************************/
        //local 2D vector for computing of probability for each angle bin i and event bin u that particle is in the ith bin
        //given the event is in the event bin u(probability_i_in_u)
        Events probability_i_in_u(bins, Row(sumAngleBinsByIFields));

        //setup all values of probability_i_in_u to 0
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            std::fill(probability_i_in_u[i].begin(), probability_i_in_u[i].end(), 0);
        }

        //computing of probability_i_in_u
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int k = 0; k < multiplicityField; k++) {
                    probability_i_in_u[i][k] += binsMatrixOrig[i][j][k];
                }
                probability_i_in_u[i][sumAngleBinsByIFields - 1] += binsMatrixOrig[i][j][multiplicityField];
            }
        }

        /*********************************************************************************************************************/

        //local 3D vector for Calculate the probability that an event in bin u is described by a set of numbers {n} denoting
        //ni particles belong to ith bin (probability_ni_u)
        Bins probability_ni_u(bins, Events(eventsPerBin, Row(bins)));

        //setup all values of probability_ni_u to 1
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                std::fill(probability_ni_u[i][j].begin(), probability_ni_u[i][j].end(), 1);
            }
        }

        //computing of probability_ni_u by use of values from probability_i_in_u
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int l = 0; l < bins; l++) {
                    for (int k = 0; k < multiplicityField; k++) {
                        long double p_i_u = 0;
                        if (probability_i_in_u[l][sumAngleBinsByIFields - 1] != 0) {
                            p_i_u = probability_i_in_u[l][k] / probability_i_in_u[l][sumAngleBinsByIFields - 1];
                        }
                        long double pow_p_i_u = pow(p_i_u, binsMatrixOrig[i][j][k]);
                        probability_ni_u[i][j][l] *= pow_p_i_u;
                    }
                }
            }
        }

        /*********************************************************************************************************************/

        //prepare local 3D vector for calculating of probability for each event with record {n} that it belongs to the
        //bin u (probability_u_ni)
        Bins probability_u_ni(bins, Events(eventsPerBin, Row(bins)));

        //computing of probability_u_ni by use of values from probability_ni_u
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int l = 0; l < bins; l++) {
                    long double denominator = 0.0;
                    for (int k = 0; k < bins; k++) {
                        denominator += (probability_ni_u[i][j][k] * (1.0 / bins));
                    }
                    probability_u_ni[i][j][l] = (probability_ni_u[i][j][l] * (1.0 / bins)) / denominator;
                }
            }
        }

        /*********************************************************************************************************************/

        //local 2D vector for calculating of average bin number for every event (a_u)
        Events a_u(bins, Row(eventsPerBin));

        //computing of a_u be use of values from probability_u_ni
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                a_u[i][j] = 0;
                for (int k = 0; k < bins; k++) {
                    a_u[i][j] += ((k + 1) * probability_u_ni[i][j][k]);
                }
            }
        }

        /*********************************************************************************************************************/

        //assign average bin numbers to events in event bins in global 3D vector - new sorting value
#pragma omp parallel for num_threads(10)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                binsMatrixOrig[i][j][sorterField] = a_u[i][j];
            }
        }

        //end of algorithm for computing average bin numbers
        /*********************************************************************************************************************/

        //prepare 2D vector from global 3D vector for sorting by sorter value
        Events eventsMatrix(events, Row(eventFields));
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int k = 0; k < eventFields; k++) {
                    eventsMatrix[(i * eventsPerBin) + j][k] = binsMatrixOrig[i][j][k];
                }
            }
        }

        //sort by comparator method
        std::sort(eventsMatrix.begin(), eventsMatrix.end(), compare);

        //local 3D vector for distribution of sorted events to event bins
        Bins binsMatrix(bins, Events(eventsPerBin, Row(eventFields)));

        int actualEvent = 0;
        //distributing of events to event bins & checking of change of bin id for event from previous sorted events
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int k = 0; k < eventFields; k++) {
                    binsMatrix[i][j][k] = eventsMatrix[actualEvent][k];
                }
                if (binsMatrix[i][j][myBinId] != i) {
                    binsMatrix[i][j][myBinId] = i;
                    sorter = true;
                    diff += 1;
                }
                actualEvent += 1;
            }
        }

        //assigning of distributed events in event bins to global 3D vector
        binsMatrixOrig.assign(binsMatrix.begin(), binsMatrix.end());

        //clean up of local objects to release memory
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                binsMatrix[i][j].clear();
            }
            binsMatrix[i].clear();
        }
        binsMatrix.clear();
        for (int i = 0; i < events; i++) {
            eventsMatrix[i].clear();
        }
        eventsMatrix.clear();
        //end of clean up

        if (events < 10000) {
            if (cyc % 1000 == 0) {
                cout << cyc << ". cycle done. " << endl;
            }
        } else {
            if (events >= 50000) {
                cout << cyc << ". cycle done. " << endl;
            } else {
                if (cyc % 100 == 0) {
                    cout << cyc << ". cycle done. " << endl;
                }
            }
        }
        if (cyc == 1) {
            previousDiff = diff;
        } else {
            if (probablySorted == probablySortedMax) {
                cyc -= 10;
                break;
            }
            if (previousDiff == diff) {
                probablySorted += 1;
            } else {
                probablySorted = 0;
            }
            previousDiff = diff;
        }

        diff = 0;

        if (cyc == 5000) {
            break;
        }

        //cycle timer - END
        clock_t endcyc = clock();
        //computing of elapsed time of one cycle (in seconds)
        double elapsed_cyc = double(endcyc - begin_cyc) / CLOCKS_PER_SEC;
        //cout << "ELAPSED TIME PER CYC: " << elapsed_cyc << " :: Cycle= " << cyc << endl;
    }

    //sorting timer - END
    clock_t end = clock();
    //computing of elapsed time of sorting (in seconds)
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "ELAPSED TIME : " << elapsed_secs << " :: NUMBRE OF CYCLES = " << cyc << endl;

    //writing results to file like average of empirical distributions of events of every event bin and like all empirical
    //distributions of events of every event bin
    cout << "Writing results to file.." << endl;

    //prepare folder for results
    BasePath bp;
    string basePath = bp.getBasePath();
    string resultFolder = basePath + RESULT;
    CreateFolder cf(resultFolder);

    ofstream result_f;
    ofstream result_a_f;
    ofstream result_oa_f;
    string binPathOvAv = basePath + RESULT + "av_bin";
    result_oa_f.open(binPathOvAv.c_str());

    for (int i = 0; i < bins; i++) {

        //prepare file for average of empirical distributions of events of every event bin
        string binPath = basePath + RESULT + "bin_" + to_string(i + 1);
        //prepare file for all empirical distributions of events of every event bin
        string binPathAv = basePath + RESULT + "av_bin_" + to_string(i + 1);

        int av_u[multiplicityField];
        for (int m = 0; m < multiplicityField; m++) {
            av_u[m] = 0;
        }

        result_f.open(binPath.c_str());

        for (int j = 0; j < eventsPerBin; j++) {
            result_f << j + 1 << "  ";
            for (int k = 0; k < multiplicityField; k++) {
                av_u[k] += binsMatrixOrig[i][j][k];
                result_f << binsMatrixOrig[i][j][k] << "  ";
            }
            result_f << binsMatrixOrig[i][j][sorterField] << endl;
        }

        result_a_f.open(binPathAv.c_str());

        for (int l = 0; l < multiplicityField; l++) {
            result_a_f << av_u[l] / eventsPerBin << setw(4);
            result_oa_f << av_u[l] / eventsPerBin << setw(4);
        }

        result_a_f << endl;
        result_oa_f << endl;
        result_a_f.close();
        result_f.close();
    }
    result_oa_f.close();
}

/*
 * Print out angle bins of every event
 */
void Sorter::printAngleBins() {
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < eventsPerBin; j++) {
            for (int k = 0; k < eventFields; k++) {
                cout << binsMatrixOrig[i][j][k] << " | ";
            }
            cout << endl;
        }
    }
}

/*
 * Print out sorter of every event of event bins
 */
void Sorter::printEvents() {
    int actualEvent = 0;
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < eventsPerBin; j++) {
            if ((j + 1) % 10 == 0) {
                cout << "EVENT " << actualEvent + 1 << " :: BIN : " << binsMatrixOrig[i][j][myBinId]
                        << " :: SORTED BY : " << binsMatrixOrig[i][j][sorterField] << endl;
            }
            actualEvent += 1;
        }
    }
}

Sorter::~Sorter() {
}

/*
 * Sorter.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#include "Sorter.h"

bool compare(vector<long double> row_a, vector<long double> row_b) {
    return row_a[21] > row_b[21];
}

Sorter::Sorter() {
    BasePath bp;
    string removeOld = "rm -r " + bp.getBasePath() + RESULT;
    system(removeOld.c_str());

    cout << "--------------------------------------------------" << endl;
    cout << "Loading data.." << endl;

    string path = bp.getBasePath() + GEN_EVENTS_FILE_PATH;
    std::ifstream ifile_ev(path.c_str(), std::ios::in);

    events = 0;
    int actualM = 0;
    std::vector<int> m;
    while (ifile_ev >> events) {
        ifile_ev >> actualM;
        m.push_back(actualM);
    }
    ifile_ev.close();

    eventsPerBin = events / bins;

    cout << "EVENTS = " << events << "  ::  EVENTS PER BIN = " << eventsPerBin << endl;
    Events eventsMatrix(events, Row(eventFields));

    for (int i = 0; i < events; i++) {
        for (int j = 0; j < eventFields; j++) {
            eventsMatrix[i][j] = 0.0;
        }
        eventsMatrix[i][multiplicityField] = m[i];
    }
    m.clear();

    int event = 1;
    long double angle = 0.0;
    long double q_x = 0.0;
    long double q_y = 0.0;
    long double q_2 = 0.0;

    int min = 1;
    int max = 1000;

    for (int i = 0; i < events; i++) {
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
        std::ifstream ifile(path.c_str(), std::ios::in);

        while (ifile >> event) {
            ifile >> angle;
            q_x += cos(2.0 * angle);
            q_y += sin(2.0 * angle);
            for (int j = 0; j < 20; j++) {
                if (angle >= (baseAngle * j) && angle < (baseAngle * (j + 1.0))) {
                    eventsMatrix[i][j] += 1;
                    break;
                }
            }
        }

        ifile.close();

        q_2 = sqrt(pow(q_x, 2) + pow(q_y, 2)) / sqrt(eventsMatrix[i][multiplicityField]);
        eventsMatrix[i][sorterField] = q_2;
        q_x = 0;
        q_y = 0;
        if ((i + 1) % 1000 == 0) {
            cout << event << '/' << events << endl;
        }
    }

    cout << "Data loaded." << endl;
    cout << "--------------------------------------------------" << endl;

    std::sort(eventsMatrix.begin(), eventsMatrix.end(), compare);

    Bins binsMatrix(bins, Events(eventsPerBin, Row(eventFields)));

    int actualEvent = 0;
    for (int i = 0; i < bins; i++) {
        for (int j = 0; j < eventsPerBin; j++) {
            for (int k = 0; k < eventFields; k++) {
                binsMatrix[i][j][k] = eventsMatrix[actualEvent][k];
                if (k == myBinId) {
                    binsMatrix[i][j][k] = i;
                }
            }
            actualEvent += 1;
        }
    }

    binsMatrixOrig.assign(binsMatrix.begin(), binsMatrix.end());

    //clean up
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
    //clean end
}

void Sorter::sort() {
    //int iam, np;
    //unsigned long long sum = 0;
    //printAngleBins();

    clock_t begin = clock();
    int cyc = 0;
    int diff = 0;
    int previousDiff = -1;
    int probablySorted = 0;
    int probablySortedMax = 10;
    bool sorter = true;
    while (sorter) {
        sorter = false;
        cyc++;

        Events probability_i_in_u(bins, Row(sumAngleBinsByIFields));
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            std::fill(probability_i_in_u[i].begin(), probability_i_in_u[i].end(), 0);
        }
//        long double probability_i_in_u[bins][sumAngleBinsByIFields];
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int k = 0; k < multiplicityField; k++) {
//                    if (j == 0) {
//                        probability_i_in_u[i][k] = 0;
//                        probability_i_in_u[i][sumAngleBinsByIFields - 1] = 0;
//                    }
                    probability_i_in_u[i][k] += binsMatrixOrig[i][j][k];
                }
                probability_i_in_u[i][sumAngleBinsByIFields - 1] += binsMatrixOrig[i][j][multiplicityField];
            }
        }
        /*********************************************************************************************************************/
//        long double probability_ni_u[bins][eventsPerBin][bins];
        Bins probability_ni_u(bins, Events(eventsPerBin, Row(bins)));
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                std::fill(probability_ni_u[i][j].begin(), probability_ni_u[i][j].end(), 1);
            }
        }
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int l = 0; l < bins; l++) {
//                    probability_ni_u[i][j][l] = 1;
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
        Bins probability_u_ni(bins, Events(eventsPerBin, Row(bins)));
//        long double probability_u_ni[bins][eventsPerBin][bins];
#pragma omp parallel for schedule(dynamic, 2)
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
//        long double a_u[bins][eventsPerBin];
        Events a_u(bins, Row(eventsPerBin));
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                a_u[i][j] = 0;
                for (int k = 0; k < bins; k++) {
                    a_u[i][j] += ((k + 1) * probability_u_ni[i][j][k]);
                }
            }
        }
        /*********************************************************************************************************************/
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                binsMatrixOrig[i][j][sorterField] = a_u[i][j];
            }
        }

        Events eventsMatrix(events, Row(eventFields));
#pragma omp parallel for schedule(dynamic, 2)
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
                for (int k = 0; k < eventFields; k++) {
                    eventsMatrix[(i * eventsPerBin) + j][k] = binsMatrixOrig[i][j][k];
                }
            }
        }

//        cout << "*********************************************" << endl;
//        for (int i = 0; i < events; i++) {
//            cout << i + 1 << ". :: BIN " << eventsMatrix[i][myBinId] << " :: SORTER " << eventsMatrix[i][sorterField]
//                    << endl;
//        }
        std::sort(eventsMatrix.begin(), eventsMatrix.end(), compare);
//        cout << "*********************************************" << endl;
//        for (int i = 0; i < events; i++) {
//            cout << i + 1 << ". :: BIN " << eventsMatrix[i][myBinId] << " :: SORTER " << eventsMatrix[i][sorterField]
//                    << endl;
//        }
//        cout << "*********************************************" << endl;

        Bins binsMatrix(bins, Events(eventsPerBin, Row(eventFields)));

        int actualEvent = 0;
        for (int i = 0; i < bins; i++) {
            for (int j = 0; j < eventsPerBin; j++) {
//                cout << eventsMatrix[actualEvent][myBinId] << " : " << binsMatrixOrig[i][j][myBinId] << endl;
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

        binsMatrixOrig.assign(binsMatrix.begin(), binsMatrix.end());

        //clean up
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
        //clean end

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
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout << "ELAPSED TIME : " << elapsed_secs << " :: NUMBRE OF CYCLES = " << cyc << endl;
    cout << "Writing results to file.." << endl;
    BasePath bp;
    string basePath = bp.getBasePath();
    string resultFolder = basePath + RESULT;
    CreateFolder cf(resultFolder);

    ofstream result_f;
    ofstream result_a_f;
    for (int i = 0; i < bins; i++) {
        string binPath = basePath + RESULT + "bin_" + to_string(i + 1);
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
        }
        result_a_f << endl;
        result_a_f.close();
        result_f.close();
    }
}

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

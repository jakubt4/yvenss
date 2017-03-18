/*
 * Generator.cpp
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#include "Generator.h"

/**
 * Prepare object for generating events
 * @param _events - number of events to generate
 */
Generator::Generator(int _events) {
    events = _events;
    BasePath bp;
    basePath = bp.getBasePath() + GEN_PATH;
    string removeOld = "rm -r " + basePath;
    system(removeOld.c_str());
    CreateFolder cr(basePath);
    srand(time(NULL));
}

int seed() {
    srand(time(0));
    int s = rand();
    return s;
}

struct random gen(seed());

double generateMultiplicity(int events) {
    int max = sqrt(events);
    int min = events / 10;
    if (min == 0) {
        min = 1;
    }
    return 300 + gen.doub() * 2700;
}

double gauss(double mu) {
    double u, v, x, y, q;
    do {
        u = gen.doub();
        v = 1.7156 * (gen.doub() - 0.5);
        x = u - 0.449871;
        y = abs(v) + 0.386595;
        q = pow(x, 2) + y * (0.19600 * y - 0.25472 * x);
    } while (q > 0.27597 && (q > 0.27846 || pow(v, 2) > -4 * log(u) * pow(u, 2)));
    return mu + 0.025 * v / u;
}

double phi() {
    return 2 * 3.14159265358979323846 * gen.doub();
}

double run() {
    return gen.doub();
}

double gen_fi(double fn[2], double vn[2]) {

    int i = 0;
    double p, f;
    double x = 0;
    double fi = 0;

    f = 1;
    for (i = 0; i < 2; i++) {
        f = f + 2 * vn[i];
    }

    do {
        x = run();
        p = 1;
        fi = phi();
        for (i = 0; i < 2; i++) {
            p = p + 2 * vn[i] * cos((i + 1) * (fi - fn[i]));
        }
    } while (x > p / f);

    return fi;
}

/**
 * Generate events:
 * events -> /generated_files/generated_events
 * particles -> /generated_files/generated_particles
 */
void Generator::generate() {
    cout << "--------------------------------------------------" << endl;
    cout << "Generating events..." << endl;

    if (events >= 1000) {
        string newPath = basePath + "generated_particles_1_1000/";
        CreateFolder cr(newPath);
    }

    ofstream events_f;
    string eventsPath = basePath + GEN_EVENTS;
    events_f.open(eventsPath.c_str());
    int min = 1;
    int max = 1000;

    for (int i = 0; i < events; i++) {
        int multiplicity = round(generateMultiplicity(events));
        events_f << i + 1 << setw(15) << multiplicity << endl;
        for (int j = 0; j < 2; j++) {
            vn0[j] = params[j][0] * (pow((75.0 / 2700.0), 2)) * multiplicity * multiplicity
                    + (-params[j][1] * 75.0 / 2700.0 - params[j][0] * 125.0 / 27.0) * multiplicity + params[j][2]
                    + params[j][0] * (250.0 * 250.0 / 9.0) + params[j][1] * (250.0 / 3.0);
            vn[j] = gauss(vn0[j]);
            fn[j] = phi();
        }

        ofstream particles_f;
        string particlesPath;
        if (events >= 1000) {
            if (i < 1000) {
                particlesPath = basePath + "generated_particles_1_1000/" + GEN_PARTICLES + to_string(i + 1);
                particles_f.open(particlesPath.c_str());
            } else {
                if ((i % 1000) == 0) {
                    min = i + 1;
                    max = i + 1000;
                    string newPath = basePath + "generated_particles_" + to_string(min) + "_" + to_string(max) + "/";
                    CreateFolder cr(newPath);
                }
                particlesPath = basePath + "generated_particles_" + to_string(min) + "_" + to_string(max) + "/"
                        + GEN_PARTICLES + to_string(i + 1);
                particles_f.open(particlesPath.c_str());
            }
        } else {
            particlesPath = basePath + GEN_PARTICLES + to_string(i + 1);
            particles_f.open(particlesPath.c_str());
        }

        for (int k = 0; k < multiplicity; k++) {
            fi = gen_fi(fn, vn);
            particles_f << i + 1 << setw(15) << fi << endl;
        }
        particles_f.close();

        if ((i + 1) % 1000 == 0) {
            cout << i + 1 << '/' << events << endl;
        }

    }
    events_f.close();

    cout << "All events generated." << endl;
    cout << "--------------------------------------------------" << endl;
}

Generator::~Generator() {
}

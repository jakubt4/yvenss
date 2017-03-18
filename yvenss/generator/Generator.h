/*
 * Generator.h
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#ifndef GENERATOR_H_
#define GENERATOR_H_

#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <fstream>
#include <iomanip>

#include "../files/headers/BasePath.h"
#include "../files/headers/CreateFolder.h"


    class Generator {
    private:
        int events;
        string basePath;
        double fn[2];
        double vn[2];
        double vn0[2];
        double fi = 0;
        double params[2][3] = { { 0, -0.000006, -0.00018 }, { -0.000092, 0.00811, 0.0418 } };
    public:
        Generator(int _events);
        void generate();
        virtual ~Generator();
    };

typedef unsigned long long int Ullong;

struct random {
        Ullong u, v, w;
        random(Ullong j) :
                v(4101842887655102017LL), w(1)
        //constructor
        {
            u = j ^ v;
            int64();
            v = u;
            int64();
            w = v;
            int64();
        }
        inline Ullong int64() //inline; compiler inserts complete function body when the function is called
        {
            u = u * 2862933555777941757LL + 7046029254386353087LL;
            v ^= v >> 17;
            v ^= v << 31;
            v ^= v >> 8;
            w = 4294957665U * (w & 0xffffffff) + (w >> 32);
            Ullong x = u ^ (u << 21);
            x ^= x >> 35;
            x ^= x << 4;
            return (x + v) ^ w;
        }
        inline double doub() {
            return 5.42101086242752217E-20 * int64();
        }
};

#endif /* GENERATOR_H_ */

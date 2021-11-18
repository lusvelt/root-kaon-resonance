#include <string>

using namespace std;

#ifndef PARAMETERS_H
#define PARAMETERS_H

const int N_PARTICLE_TYPES = 7;
const int N_ITERATIONS = 1E5;
const int N_PARTICLES_PER_ITERATION = 100;
const int MAX_PRODUCTS = 200;
const double AVG_P = 1.0;
const int N_BINS = 50;
const int N_BINS_INV_MASS = 100;
const double MAX_MOMENTUM = 5.0;
const double MAX_ENERGY = 8.0;
const double MIN_INVARIANT_MASS = 0.5;
const double MAX_INVARIANT_MASS = 1.5;

const int PION_PLUS_BIN = 1;
const int PION_MINUS_BIN = 2;
const int KAON_PLUS_BIN = 3;
const int KAON_MINUS_BIN = 4;
const int PROTON_PLUS_BIN = 5;
const int PROTON_MINUS_BIN = 6;
const int KAON_STAR_BIN = 7;

const double PION_PLUS_PROB = 0.4;
const double PION_MINUS_PROB = 0.4;
const double KAON_PLUS_PROB = 0.05;
const double KAON_MINUS_PROB = 0.05;
const double PROTON_PLUS_PROB = 0.045;
const double PROTON_MINUS_PROB = 0.045;
const double KAON_STAR_PROB = 0.01;

const double PROBABILITIES[] = {
    PION_PLUS_PROB,
    PION_MINUS_PROB,
    KAON_PLUS_PROB,
    KAON_MINUS_PROB,
    PROTON_PLUS_PROB,
    PROTON_MINUS_PROB,
    KAON_STAR_PROB
};

const string PION_PLUS_LABEL = "π+";
const string PION_MINUS_LABEL = "π-";
const string KAON_PLUS_LABEL = "K+";
const string KAON_MINUS_LABEL = "K-";
const string PROTON_PLUS_LABEL = "p+";
const string PROTON_MINUS_LABEL = "p-";
const string KAON_STAR_LABEL = "K*";

const string LABELS[] = {
    PION_PLUS_LABEL,
    PION_MINUS_LABEL,
    KAON_PLUS_LABEL,
    KAON_MINUS_LABEL,
    PROTON_PLUS_LABEL,
    PROTON_MINUS_LABEL,
    KAON_STAR_LABEL
};

const double PION_PLUS_CUMULATIVE = PION_PLUS_PROB;
const double PION_MINUS_CUMULATIVE = PION_PLUS_CUMULATIVE + PION_MINUS_PROB;
const double KAON_PLUS_CUMULATIVE = PION_MINUS_CUMULATIVE + KAON_PLUS_PROB;
const double KAON_MINUS_CUMULATIVE = KAON_PLUS_CUMULATIVE + KAON_MINUS_PROB;
const double PROTON_PLUS_CUMULATIVE = KAON_MINUS_CUMULATIVE + PROTON_PLUS_PROB;
const double PROTON_MINUS_CUMULATIVE = PROTON_PLUS_CUMULATIVE + PROTON_MINUS_PROB;
const double KAON_STAR_CUMULATIVE = PROTON_MINUS_CUMULATIVE + KAON_STAR_PROB;

#endif
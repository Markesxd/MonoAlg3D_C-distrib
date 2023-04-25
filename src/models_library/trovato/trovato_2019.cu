#include "trovato_2019.h"
#include <stddef.h>
#include <stdint.h>

__global__ void kernel_set_model_initial_conditions(real *sv, int num_volumes, size_t pitch, bool use_adpt_dt, real min_dt) {
    int threadID = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadID < num_volumes) {

        real STATES[NEQ];
        
        // Steady-state 40 pulses (BCL=1000ms)
        STATES[0] = -86.7099;
        STATES[1] = 0.005431;
        STATES[2] = 0.000104;
        STATES[3] = 8.25533;
        STATES[4] = 8.25502;
        STATES[5] = 8.25503;
        STATES[6] = 143.743;
        STATES[7] = 143.744;
        STATES[8] = 143.744;
        STATES[9] = 4.4e-05;
        STATES[10] = 0.000103;
        STATES[11] = 1.26947;
        STATES[12] = 1.25254;
        STATES[13] = 1.27103;
        STATES[14] = 1.1e-05;
        STATES[15] = 0;
        STATES[16] = 0.006303;
        STATES[17] = 0.789469;
        STATES[18] = 0.789392;
        STATES[19] = 0.791301;
        STATES[20] = 0.580955;
        STATES[21] = 0.791719;
        STATES[22] = 0.000241;
        STATES[23] = 0.463851;
        STATES[24] = 0.239936;
        STATES[25] = 0.000272;
        STATES[26] = 0.646362;
        STATES[27] = 0.98999;
        STATES[28] = 0;
        STATES[29] = 1;
        STATES[30] = 0.926919;
        STATES[31] = 1;
        STATES[32] = 1;
        STATES[33] = 0.999976;
        STATES[34] = 1;
        STATES[35] = 1;
        STATES[36] = 0.005885;
        STATES[37] = 0.000303;
        STATES[38] = 0.994251;
        STATES[39] = 0.000367;
        STATES[40] = 0.566131;
        STATES[41] = 0.189842;
        STATES[42] = 0.000222;
        STATES[43] = 0.233515;
        STATES[44] = 0.997077;
        STATES[45] = 0.471259;

        for (int i = 0; i < NEQ; i++)
            *((real * )((char *) sv + pitch * i) + threadID) = STATES[i];


        if(use_adpt_dt) {
            *((real *)((char *)sv + pitch * NEQ) + threadID) = min_dt; // dt
            *((real *)((char *)sv + pitch * (NEQ + 1)) + threadID) = 0.0;    // time_new
            *((real *)((char *)sv + pitch * (NEQ + 2)) + threadID) = 0.0;    // previous dt
        }
    }
}

inline __device__ void RHS_gpu(real *sv, real *rDY_, real stim_current, int threadID_, real dt, size_t pitch, bool use_adpt_dt) {

    // Get the stimulus current from the current cell
    real calc_I_stim = stim_current;

    // State variables
    real V;
    real CaMKt;
    real cass;
    real nai;
    real nasl;
    real nass;
    real ki;
    real kss;
    real ksl;
    real cai;
    real casl;
    real cansr;
    real cajsr;
    real cacsr;
    real Jrel1;
    real Jrel2;
    real m;
    real hf;
    real hs;
    real j;
    real hsp;
    real jp;
    real mL;
    real hL;
    real hLp;
    real a;
    real i1;
    real i2;
    real d;
    real ff;
    real fs;
    real fcaf;
    real fcas;
    real jca;
    real ffp;
    real fcafp;
    real nca;
    real b;
    real g;
    real xrf;
    real xrs;
    real xs1;
    real xs2;
    real y;
    real xk1;
    real u;   

    if (use_adpt_dt)
    {
        V = sv[0];
        CaMKt = sv[1];
        cass = sv[2];
        nai = sv[3];
        nasl = sv[4];
        nass = sv[5];
        ki = sv[6];
        kss = sv[7];
        ksl = sv[8];
        cai = sv[9];
        casl = sv[10];
        cansr = sv[11];
        cajsr = sv[12];
        cacsr = sv[13];
        Jrel1 = sv[14];
        Jrel2 = sv[15];
        m = sv[16];
        hf = sv[17];
        hs = sv[18];
        j = sv[19];
        hsp = sv[20];
        jp = sv[21];
        mL = sv[22];
        hL = sv[23];
        hLp = sv[24];
        a = sv[25];
        i1 = sv[26];
        i2 = sv[27];
        d = sv[28];
        ff = sv[29];
        fs = sv[30];
        fcaf = sv[31];
        fcas = sv[32];
        jca = sv[33];
        ffp = sv[34];
        fcafp = sv[35];
        nca = sv[36];
        b = sv[37];
        g = sv[38];
        xrf = sv[39];
        xrs = sv[40];
        xs1 = sv[41];
        xs2 = sv[42];
        y = sv[43];
        xk1 = sv[44];
        u = sv[45];
    }
    else
    {
        V = *((real *)((char *)sv + pitch * 0) + threadID_);
        CaMKt = *((real *)((char *)sv + pitch * 1) + threadID_);
        cass = *((real *)((char *)sv + pitch * 2) + threadID_);
        nai = *((real *)((char *)sv + pitch * 3) + threadID_);
        nasl = *((real *)((char *)sv + pitch * 4) + threadID_);
        nass = *((real *)((char *)sv + pitch * 5) + threadID_);
        ki = *((real *)((char *)sv + pitch * 6) + threadID_);
        kss = *((real *)((char *)sv + pitch * 7) + threadID_);
        ksl = *((real *)((char *)sv + pitch * 8) + threadID_);
        cai = *((real *)((char *)sv + pitch * 9) + threadID_);
        casl = *((real *)((char *)sv + pitch * 10) + threadID_);
        cansr = *((real *)((char *)sv + pitch * 11) + threadID_);
        cajsr = *((real *)((char *)sv + pitch * 12) + threadID_);
        cacsr = *((real *)((char *)sv + pitch * 13) + threadID_);
        Jrel1 = *((real *)((char *)sv + pitch * 14) + threadID_);
        Jrel2 = *((real *)((char *)sv + pitch * 15) + threadID_);
        m = *((real *)((char *)sv + pitch * 16) + threadID_);
        hf = *((real *)((char *)sv + pitch * 17) + threadID_);
        hs = *((real *)((char *)sv + pitch * 18) + threadID_);
        j = *((real *)((char *)sv + pitch * 19) + threadID_);
        hsp = *((real *)((char *)sv + pitch * 20) + threadID_);
        jp = *((real *)((char *)sv + pitch * 21) + threadID_);
        mL = *((real *)((char *)sv + pitch * 22) + threadID_);
        hL = *((real *)((char *)sv + pitch * 23) + threadID_);
        hLp = *((real *)((char *)sv + pitch * 24) + threadID_);
        a = *((real *)((char *)sv + pitch * 25) + threadID_);
        i1 = *((real *)((char *)sv + pitch * 26) + threadID_);
        i2 = *((real *)((char *)sv + pitch * 27) + threadID_);
        d = *((real *)((char *)sv + pitch * 28) + threadID_);
        ff = *((real *)((char *)sv + pitch * 29) + threadID_);
        fs = *((real *)((char *)sv + pitch * 30) + threadID_);
        fcaf = *((real *)((char *)sv + pitch * 31) + threadID_);
        fcas = *((real *)((char *)sv + pitch * 32) + threadID_);
        jca = *((real *)((char *)sv + pitch * 33) + threadID_);
        ffp = *((real *)((char *)sv + pitch * 34) + threadID_);
        fcafp = *((real *)((char *)sv + pitch * 35) + threadID_);
        nca = *((real *)((char *)sv + pitch * 36) + threadID_);
        b = *((real *)((char *)sv + pitch * 37) + threadID_);
        g = *((real *)((char *)sv + pitch * 38) + threadID_);
        xrf = *((real *)((char *)sv + pitch * 39) + threadID_);
        xrs = *((real *)((char *)sv + pitch * 40) + threadID_);
        xs1 = *((real *)((char *)sv + pitch * 41) + threadID_);
        xs2 = *((real *)((char *)sv + pitch * 42) + threadID_);
        y = *((real *)((char *)sv + pitch * 43) + threadID_);
        xk1 = *((real *)((char *)sv + pitch * 44) + threadID_);
        u = *((real *)((char *)sv + pitch * 45) + threadID_);
    }

    #include "trovato_2019_common.inc"
}

//Include the default solver used by all models.
#include "../default_solvers.cu"

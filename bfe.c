#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bfe.h"

const int T[] = {5, 10, 20, 40, 80, 160};

const double e_alpha = -0.0114344987526;
const double e_beta= 1.3953731963283147;
const double e_gamma = 1.98E-6;

// weigth of the exponential moving average
// rule is 2/(n+1) so to keep account last ~10 values, p_avg = 0.18
static const double p_avg = 0.18; 

// before MIN_INTERVALS_FOR_ESTIMATION we keep a fast duty cycle to speed-up learning
#define MIN_INTERVALS_FOR_ESTIMATION 4

static double n[N][N]; // BFE matrix

static double bw[N][N]; // nodes count data to compute the backward estimation
static double fw[N][N]; // nodes count data to compute the forward estimation

static double m_hat[N]; // estimation of m(T)
static double p_hat[N]; // estimation of p_j

static double p_count[N]; // number of nodes with a given period

static double m_count[N]; // nodes discovered by me o na given period

static int ints[N]; // how many intervals passed available for average

static double energy_in_period(duty_cycle_idx t);
static void estimation_nt(double n_hat[N], duty_cycle_idx a);
static duty_cycle_idx find_maximum_T(double n_hat[N]);
static duty_cycle_idx adjust_T(duty_cycle_idx T_current, duty_cycle_idx T_max);

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// initialize data structures
void init_BFE() {
    memset(n, 0, N*N*sizeof(double));
    memset(bw, 0, N*N*sizeof(double));
    memset(fw, 0, N*N*sizeof(double));
    memset(m_hat, 0, N*sizeof(double));
    memset(p_hat, 0, N*sizeof(double));
    memset(m_count, 0, N*sizeof(double));
    memset(p_count, 0, N*sizeof(double));
    memset(ints, 0, N*sizeof(int));
}

/* Parameter: my (a) and other (b) duty cycle indexes
 * Number of packet sent by me (m_a) and by the other node during the encounter
 */
void onContactFinished(duty_cycle_idx a, duty_cycle_idx b, int m_a, int m_b) {
    int i;
    for (i=a; i<N; i++) {
        if (i >= b) {
            // Forward
            fw[i][b] += 1 - MAX(0, (1 - m_a/pow(2, i-a))) * MAX(0, (1- m_b/pow(2, i-b)));
        } else {
            // Backward
            if (m_b > 0)
                // B would discover me in any case
                bw[i][b] += 1;
            else
                // If I change my duty cycle from "a" to "i", 
                // would I have discover B the same?
                bw[i][b] += 1 - MAX(0, 1 - m_a/pow(2, i-a));
        }

        // number of nodes that I would have discover just because of my beaconing, 
        // disregarding the beaconing coming from other nodes. Used for BW estimation.
        if (m_a > 0) { 
            m_count[i] += 1 - MAX(0, 1 - m_a/pow(2, i-a));
        }

       
    }

    // calculate the distribution of the nodes in the environment according to their different duty cycles
    if (m_a > 0) {
        p_count[b] += 1;
    }
}

/*
 * Returns the energy consumed in period T according to eq. 6
 */
double energy_in_period(duty_cycle_idx t) {
    return e_alpha + e_beta * sqrt(T[t]) + e_gamma * T[t];
}

/* Estimate n(T)
 *  
 *
 */
void estimation_nt(double n_hat[N], duty_cycle_idx a) {
    int i, j, cointained_intervals;
    double c, p_hat_tot = 0.0;

    // update the number of intervals elapsed since the beginning
    for (i=a; i<N; i++) {
        ints[i] += (int)pow(2, N_MAX-i);
    }
   
    
    // calculate m_hat
    for (i=a; i<N; i++) {
        cointained_intervals = (int)pow(2, N_MAX - i); // how many intervals of type T_i are in T_{N_MAX} 

        //m_hat[i] =  (m_hat[i] * (ints[i] - cointained_intervals) + m_count[i]) / (double)ints[i];
        //
        // exponential moving average
        if (ints[N_MAX] == 1) {
            // first time!
            m_hat[i] = m_count[i]/cointained_intervals;
        } else {
            m_hat[i] = (1.0-p_avg) * m_hat[i]  + p_avg * m_count[i]/cointained_intervals;
        }

        m_count[i] = 0;
    }

    // calculate p_hat
    for (i=0; i<N; i++) {
        //p_hat[i] = (p_hat[i] * (ints[N_MAX] - 1) + p_count[i]) / (double)ints[N_MAX];
        // exponential moving average
        if (ints[N_MAX] == 1) {
            p_hat[i] = p_count[i];
        } else {
            p_hat[i] = (1.0-p_avg) * p_hat[i] + p_avg * p_count[i];
        }

        p_hat_tot += p_hat[i];
        p_count[i] = 0;
    }

    /*
    for (i=0; i<N; i++) 
        printf("p_hat[%d]=%lf ", i, p_hat[i]);
    printf("\n");

    for (i=0; i<N; i++) 
        printf("m_hat[%d]=%lf ", i, m_hat[i]);
    printf("\n");
    */

    // populate n using backward and forward data
    for (i=a; i<N; i++) {
        cointained_intervals = (int)pow(2, N_MAX - i); // how many intervals of type T_i are in T_{N_MAX} 
        for (j=0; j<N; j++) {
            if (i>=j) {
                // FW
                //n[i][j] = (n[i][j] * (ints[i] - cointained_intervals) + fw[i][j]) / (double)ints[i];
                // exponential moving average
                if (ints[N_MAX] == 1)
                    n[i][j] = fw[i][j]/cointained_intervals;
                else
                    n[i][j] = (1.0-p_avg) * n[i][j] + p_avg * fw[i][j]/cointained_intervals;

            } else {
                // BW   
                c = bw[i][j]/(double)cointained_intervals * pow(2, j-i) - (m_hat[i] * (pow(2, j-i) - 1)) * p_hat[j] / p_hat_tot;

                //n[i][j] = (n[i][j] * (ints[i] - cointained_intervals) + c*(double)cointained_intervals) / (double)ints[i];
                // exponential moving average
                if (ints[N_MAX] == 1)
                    n[i][j] = c;
                else
                    n[i][j] = (1.0-p_avg) * n[i][j] + p_avg * c;

            }
        }
    }

    // reset bw and fw
    memset(bw, 0, N*N*sizeof(double));
    memset(fw, 0, N*N*sizeof(double));

    // populate n_hat
    for (i=0; i<N; i++) {
        n_hat[i] = 0;
        for (j=0; j<N; j++) 
            n_hat[i] += n[i][j];
    }

}


/*
 * Returns the duty cycle index corresponding to the maximum value
 *
 *
 * */
duty_cycle_idx find_maximum_T(double n_hat[N]) {
    duty_cycle_idx i = 0;
    double n_max = -1, value;
    duty_cycle_idx n_max_idx = -1;
    for (i=0; i<N; i++)  {
        value = n_hat[i] / energy_in_period(i);
        if (n_max < value) {
            n_max = value;
            n_max_idx = i;
        }
    }
    return n_max_idx;
}

/*
 * Returns the next duty cycle index
 *
 * T_current: the current duty cycle
 * T_max: the maximum (best) duty cycle
 *
 * returns: the next duty cycle index
 */
duty_cycle_idx adjust_T(duty_cycle_idx T_current, duty_cycle_idx T_max) {
    if (T_max == T_current + 1 )
        return T_current;
    else if (T_max > T_current + 1)
        return MIN(N_MAX, T_current + 1);
    else
        return MAX(0, T_current - 1);
}

/*
 * Run BFE and calculate the next period T
 * Called every T_N
 *
 * */
duty_cycle_idx BFE_get_T(duty_cycle_idx current_T) {
    double n_hat[N] = {0.0}; // estimation of n(t)
    duty_cycle_idx max_t = -1;
    duty_cycle_idx next_T = -1, i, j;

    //if (ints[N_MAX]< MIN_INTERVALS_FOR_ESTIMATION) {
    //    return 0; // return the shortest interval
    //}

    estimation_nt(n_hat, current_T);
    //for (i=0; i<N; i++) 
    //    printf("n_hat[%d]=%lf ", i, n_hat[i]);
    //printf("\n");
    max_t = find_maximum_T(n_hat);
    //printf("find max_t = %d\n", max_t);
    next_T = adjust_T(current_T, max_t);
    //printf("find next_t = %d\n", next_T);

    return next_T;
}

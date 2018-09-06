#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bfe.h"

/* CONTANTS */
const int T[] = {5, 10, 20, 40, 80, 160};

const double e_alpha = -0.0114344987526;
const double e_beta= 1.3953731963283147;
const double e_gamma = 1.98E-6;

// weigth of the exponential moving average
// rule is 2/(n+1) so to keep account last ~10 values, p_avg = 0.18
static const double p_avg = 0.18; 


/* Internal data structures */

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

/* initialize data structures */
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

/* 
 * Function to be called every time an encounter finished
 *
 * a: my (node A) duty cycle index
 * b: other node (node B) duty cycle index
 * m_a: number of discovery packet emitted by node A during the encounter
 * m_b: number of discovery packet emitted by node B during the encounter
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

        // number of nodes that I would have discovered just because of my beaconing, 
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

/* Poulate n(T)
 *  
 * n_hat: empty vector to populate with the estimation of the homogeneous case
 * a: Our current duty cycle index
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

        // standard average
        //m_hat[i] =  (m_hat[i] * (ints[i] - cointained_intervals) + m_count[i]) / (double)ints[i];
        
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
        // exponential moving average
        if (ints[N_MAX] == 1)
            p_hat[i] = p_count[i];
        else 
            p_hat[i] = (1.0-p_avg) * p_hat[i] + p_avg * p_count[i];

        p_hat_tot += p_hat[i];
        p_count[i] = 0;
    }


    // populate n using backward and forward data
    for (i=a; i<N; i++) {
        cointained_intervals = (int)pow(2, N_MAX - i); // how many intervals of type T_i are in T_{N_MAX} 
        for (j=0; j<N; j++) {
            if (i>=j) {
                // Forward estimation - FW

                // standard average
                //n[i][j] = (n[i][j] * (ints[i] - cointained_intervals) + fw[i][j]) / (double)ints[i];

                // exponential moving average
                if (ints[N_MAX] == 1)
                    n[i][j] = fw[i][j]/cointained_intervals;
                else
                    n[i][j] = (1.0-p_avg) * n[i][j] + p_avg * fw[i][j]/cointained_intervals;

            } else {
                // Backward estimation - BW   
                c = bw[i][j]/(double)cointained_intervals * pow(2, j-i) - (m_hat[i] * (pow(2, j-i) - 1)) * p_hat[j] / p_hat_tot;

                // standard average
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
 * n_hat: vector of the estimated average number of nodes discovered 
 *        for all the duty cycles and in an homogeneous case 
 * returns the index corresponding to the maximum value of the ratio r
 */
duty_cycle_idx find_maximum_T(double n_hat[N]) {
    duty_cycle_idx i = 0;
    double n_max = -1, r;
    duty_cycle_idx n_max_idx = -1;
    for (i=0; i<N; i++)  {
        r = n_hat[i] / energy_in_period(i);
        if (n_max < r) {
            n_max = r;
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
 * Must be called every longest duty cycle 
 *
 * T_current: the current duty cycle
 * returns: the next duty cycle index
 */
duty_cycle_idx BFE_get_T(duty_cycle_idx current_T) {
    double n_hat[N] = {0.0}; // estimation of n(t)
    duty_cycle_idx max_t = -1;
    duty_cycle_idx next_T = -1, i, j;

    // estimate n_hat(T_i) for all the T_i
    estimation_nt(n_hat, current_T);
    // find the maximum duty cycle index 
    max_t = find_maximum_T(n_hat);
    // adjust the current duty cycle
    next_T = adjust_T(current_T, max_t);

    // return the next duty cycle that must be used
    // to converge to the optimum
    return next_T;
}

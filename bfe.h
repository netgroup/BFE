#ifndef BFE_H
#define BFE_H

#define N 6 // duty cycle indexes, goes from 0 to N-1 
#define N_MAX N - 1 // maximum duty cycle index

typedef int duty_cycle_idx;

/* 
 * Callback.
 * Must be colled when every contact finish
 *
 * a: our duty cycle index (node A)
 * b: other node duty cycle index (node B)
 * m_a : number of discovery packets sent during this encounter by node A 
 * m_b : number of discovery packets sent during this encounter by node B
 * */
void onContactFinished(duty_cycle_idx a, duty_cycle_idx b, int m_a, int m_b);

/* 
 * Main estimation routine.
 * Estimate the next duty cycle index we should use to maximize the contact discovered 
 * 
 * current_T: the currently used duty cycle period
 * returns: the next duty cycle index according to BFE
 */
duty_cycle_idx BFE_get_T(duty_cycle_idx current_T);

/* initialize internal data structure */
void init_BFE();

#endif

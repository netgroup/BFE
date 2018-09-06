#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bfe.h"


struct event {
   int type;
   double start_time;
   double end_time;
   int m_a;
   int m_b;
   int idx_b;
};

typedef struct event Event;

int main() {

    // initialize data structures
    init_BFE();

    duty_cycle_idx my_period = 0, next_period;

    // read the file with the events
    FILE* stream = fopen("scenarios/events_dynamic.csv", "r");
    Event e;

    while (!feof(stream)) {
        fscanf(stream, "%d,%lf,%lf,%d,%d, %d", &e.type, &e.start_time, &e.end_time, &e.m_a, &e.m_b, &e.idx_b);

        if (e.type == 1) {
            onContactFinished(my_period, e.idx_b-1, e.m_a, e.m_b);
        } else if (e.type == 2) {
            // this is for debugging pourposes
            // we call the evaluation function "BFE_get_T" and compare the value with some reference values
            /*
            e.m_a // new duty cycle index reference value
            e.m_b // old duty cycle index reference value
            e.end_time // current time
            */
            next_period = BFE_get_T(my_period); // call the evaluation function
            //printf("Current time: %lf\t% Reference period: d->%d\t BFE: %d->%d\n", e.end_time, e.m_b-1, e.m_a-1, my_period, next_period);
            printf("Current time: %lf\tCurrent period: %d\tNext period: %d\n", e.end_time, my_period, next_period);
            my_period = next_period;
        }
    }

    fclose(stream);

    return 0;

}

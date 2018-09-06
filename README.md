# BFE
Backward forward estimator 

This is the code of the BFE algorithm.

In bfe.h/.c there is the reference implementation with two main function: a callback that must be invoked every time a contact ends (onContactFinished) and the main estimation function (BFE_get_T) that returns the next duty cycle to use.

## Simulating a network scenario

In main.c there is the main loop. The algorithm can be tested against a set of different network scenarios whose events are reported in three files included in this repository:

* events_static.csv: 6 different duty cycles, nodes are distributed evenly across all the duty cycles
* events_dynamic.csv: the environment change the statistics
* events_homogeneous.csv: all the nodes have the maximum duty cycle (5) 

Events are reported as a tuples on a csv file format, where each line is a different event.
The first value of each event is a number that represents the event type.
Events of type 1 describes terminated contacts.
Events of type 2 are for debugging pourposes and reports some reference value of BFE together with the time in which we must call the evaluation function;


For the event type "1", the values are:
* event type (equals to 1)
* start time of the encounter
* end time of the encounter
* number of discovery packets that we sent to the other node
* number of discovery packets that the other node sent to us
* duty cycle index used by the other node


For the event type "2", the values are:
* event type (equals to 2)
* counter
* current time
* new duty cycle value id (computed with a different BFE simulator)
* old duty cycle value id (computed with a different BFE simulator)
* not used

Running the program, the algorithm print the duty cycle id in the range [0 - N) that optimize the maximum amount of discovered nodes with a fixed energy budget.

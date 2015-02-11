// Modified by Victor Fragoso <vfragoso@cs.ucsb.edu>
// 03/19/14  Adding Include header guards and namespace
#ifndef STATS_LIBRARIES_BURKARDT_SPEC_FUNCS_ASA121_H_
#define STATS_LIBRARIES_BURKARDT_SPEC_FUNCS_ASA121_H_
namespace asa121 {
void timestamp (void);
double trigamma (double x, int* ifault);
void trigamma_values (int* n_data, double* x, double* fx);
}  // asa121
#endif  // STATS_LIBRARIES_BURKARDT_SPEC_FUNCS_ASA121_H_

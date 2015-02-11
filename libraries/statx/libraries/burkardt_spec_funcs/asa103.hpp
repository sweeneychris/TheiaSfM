// Modified by Victor Fragoso <vfragoso@cs.ucsb.edu>
// 03/19/14  Adding Include header guards and namespace
#ifndef STATS_LIBRARIES_BURKARDT_SPEC_FUNCS_ASA103_H_
#define STATS_LIBRARIES_BURKARDT_SPEC_FUNCS_ASA103_H_
namespace asa103 {
double digamma (double x, int* ifault );
void psi_values (int* n_data, double* x, double* fx );
void timestamp (void);
}  // asa103
#endif  // STATS_LIBRARIES_BURKARDT_SPEC_FUNCS_ASA103_H_

//==========================================================================================
// Author Jens Chluba Aug/Sept 2010
// purpose: compute transitions rates for np, ns & nd-states of hydrogen
// last modification: Aug 2012
//==========================================================================================

#ifndef HI_TRANSITION_DATA_H
#define HI_TRANSITION_DATA_H

namespace HI_Transition_Data 
{
    //======================================================================================
    // Gamma_n=A_tot/(4 pi) for the first ten p-states (vacuum)
    //======================================================================================
    double Get_Gamma_np(int n);

    //======================================================================================
    // access to Ly-n rates (vacuum)
    //======================================================================================
    double Get_A_np1s(int n);

    //======================================================================================
    // A_npks transition rates for the first nmax p-states (vacuum)
    // k==1: just Ly-n Anp->1s rates
    // k>1 : rates into (!) the np-state
    //       n<k --> Aks-->np
    //       n>k --> 3*Anp->ks
    //======================================================================================
    double Get_A_npks(int k, int n);

    //======================================================================================
    // A_npkd transition rates for the first nmax p-states (vacuum)
    // k>1 : rates into (!) the np-state
    //       n<k --> Akd-->np
    //       n>k --> 3/5*Anp->kd
    //======================================================================================
    double Get_A_npkd(int k, int n);
}   

#endif

//==========================================================================================
//==========================================================================================
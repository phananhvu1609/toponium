/********************************************************************************
 *									  *										    *
 * 	Contains physical formulas:                                                 *
 * 	    - convert between kinematical variables									*
 * 	    - Differential cross sections for tree-level and higher order			*
 * 																				*
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/

#ifndef guard_Formulas_h
#define guard_Formulas_h

#include <iostream>
#include <math.h>
#include <complex>
#include "LHAPDF/LHAPDF.h"
#include "Utilities.h"
using namespace std;

// Compute Jacobian and costhetaCM from pT for gg > ttbar
// Warning: assume diffxsection is symmetric t <-> u
void Jacobian_and_Conversion_ggttbar_pT(double pT, double s_parton, double mt, double &costhetaCM, double &Jacobian);
// Compute Jacobian and costhetaCM from y for gg > ttbar
void Jacobian_and_Conversion_ggttbar_rapidity(double rapidity, double s_parton, double mt, double &costhetaCM, double &Jacobian);

// Cuts on xsection. True means pass cuts, false means exclude or fail cuts.
bool kinematical_cut(double costhetaCM, double s_parton, double mt);

// Compute dynamic scale
double DynamicScale(double s, double mt, double costhetaCM, int DynamicScaleChoice);

double ALP_DecayWidth(double Ctt, double fa, double mt, double Ma);

// Compute t,u from s, masses, costheta (gg > ttbar)
void costhetaCM_to_t_u_ggtt(double costhetaCM, double &t, double &u, double s, double mt);

// Tree level dxsect/dy for g + g -> t + tbar (in GeV^-2)
double dxsection_dy_ggttTree(double alphaS, double mt, double s, double y);

// Tree level dxsection/dcosthetaCM for g + g -> t + tbar. Unit GeV^-2. pCM = sqrt(s-4*mt2)
double dxsection_dcosthetaCM_ggttTree(double alphaS, double mt, double s, double costhetaCM);
double dxsection_dcosthetaCM_ggttTree_components(double alphaS, double mt, double s, double costhetaCM, double ts = 1, double tt = 1, double tu = 1);

// Toponium
double dxsect_dmToponium_partonic(int pid1, int pid2, double s_parton, double mtt, double ScaleFactorization, double ScaleRenormalization, double mt, double TopDecayWidth, double boundstateSpin, double boundstateJ, int boundstateColorConfig, LHAPDF::PDF* pdf);
double ImGreenFunction(double mtt, double mt, double TopDecayWidth, double alphaS, double CF, double CA, int boundstateColorConfig);

#endif

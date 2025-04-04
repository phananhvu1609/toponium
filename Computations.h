/********************************************************************************
 *									  *										    *
 * 	Contains computational procedures:                                          *
 * 	    - Cross section (integrate over differential ones)						*
 * 																				*
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/

#ifndef guard_Computations_h
#define guard_Computations_h

#include <iostream>
#include <cmath>
#include <fstream>
#include <gsl/gsl_errno.h>
#include "LHAPDF/LHAPDF.h"
#include "Utilities.h"
#include "Formulas.h"
using namespace std;

// Compute xsection of g + g -> t + tbar at tree-level (in GeV^-2)
double xsection_ggttTree(gsl_integration_workspace* workspace, double alphaS, double mt, double s, double &error, double costhetamin = -1, double costhetamax = 1, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);
// Compute xsection of g + g -> t + tbar at tree-level (in GeV^-2). No error returned
double xsection_ggttTree_noerr(gsl_integration_workspace* workspace, double alphaS, double mt, double s, double costhetamin = -1, double costhetamax = 1, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);
// Compute xsection of g + g -> t + tbar at tree-level (in GeV^-2) using rapidity as integration
double xsection_ggttTree_y(gsl_integration_workspace* workspace, double alphaS, double mt, double s, double &error, double ymin = HUGE_VALF, double ymax = HUGE_VALF, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);


// GSL conventions for 1/x1*PDF(x1)*PDF(s_parton/(s_hadron*x1))
double f_PDF_product(double x1, void * params);

// GSL conventions for 1/x1*PDF(x1)*PDF(s_parton/(s_hadron*x1)), u = ln(x)
double f_PDF_product_ln(double u, void *args);

// GSL convention for dxsection_dcosthetaCM_ggttTree. params[] = {alphaS, mt, s}
double f_dxsection_dcosthetaCM_ggttTree(double costhetaCM, void * params);
// GSL convention for dxsection_dy_ggttTree. params[] = {alphaS, mt, s}
double f_dxsection_dy_ggttTree(double y, void * params);

// GSL convention for Integrate over s_parton of integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection
// Warning: second integration variable must be independent of s_parton
double f_PDF_xsection(double s_parton, void *args);
// GSL convention for Integrate over costhetaCM of integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection
// Warning: second integration variable must be independent of s_parton
double f_PDF_xsection_costhetaCM(double costhetaCM, void *args);
// GSL convention for Integrate over s_parton of integrated( integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection, costhetaCM)
// Warning: second integration variable must be independent of s_parton
double f_dhadronicxsection_dsparton(double s_parton, void *args);

// ============ Data structures (of parameters passed to gsl integrand functions) ============
struct ParamPDF {
    int pid1;
    int pid2;
    double ScaleFactorization; // GeV
    double s_parton; // GeV
    double s_hadron; // GeV
};
struct IntegratePDF {
    LHAPDF::PDF* pdf;
    struct ParamPDF param;
    gsl_integration_workspace* workspace;
    double IntegrationRelativeError = NAN;    // if NAN, the default IntegrationRelativeError is used
    size_t NumberOfSubintervals = 0;    // if 0, the default NumberOfSubintervals is used
};
struct ParamXsection {
    int pid1;
    int pid2;
    double s_parton;    // GeV^2
    double alphaS;
    double mt;  // GeV
    double ScaleRenormalization; // GeV
    bool is_diffxsection = false;
    double SecondKinematicVariable = NAN;   // if NAN, SecondKinematicVarialbe is integrated over. If not, compute at that value of SecondKinematicVariable
    double MZ = 91.2; // GeV
    double vev = 246.22; // GeV
    double mtt = sqrt(s_parton);  // GeV
    double TopDecayWidth = 1.326; // GeV
    double BoundStateSpin = -1;
    double BoundStateJ = -1; 
    int BoundStateColorConfig = -1; // 1: singlet, 8: octet
};
struct IntegrateXsection {
    string name_f_xsection;       // Name of function (written in GSL convention for integrand) used to compute xsection. Either total partonic xsection or dxsection/dcosthetaCM.
    struct ParamXsection param;
    gsl_integration_workspace* workspace;
    double MinSecondKinematicVariable = NAN;    // if NAN, the default integration boundary is used
    double MaxSecondKinematicVariable = NAN;    // if NAN, the default integration boundary is used
    string name_f_change_variable = "";       // Options: "rapidity", "pT",.... Name of function used to change variable from, e.g., rapidity and pT to costhetaCM. Function must be of the form void (double pT, double s_parton, double mt, double &costhetaCM, double &Jacobian)
    double IntegrationRelativeError = NAN;    // if NAN, the default IntegrationRelativeError is used
    size_t NumberOfSubintervals = 0;    // if NAN, the default NumberOfSubintervals is used
    int DynamicScaleChoice = 0; // 0: fixed scale; 2: sum_f mT(f)/2 with mT = sqrt(m^2 + pT^2) and f is final states.
    double ScaleUncertaintyFactor = 1;  // scale *= ScaleUncertaintyFactor
    gsl_integration_workspace* workspace2 = NULL;   // in case more than 1 workspace is needed
};
struct ParamConvolute1D {
    struct IntegratePDF iPDF;
    struct IntegrateXsection iXsection;
    bool isConvolute;
    ParamConvolute1D(IntegratePDF input_iPDF, IntegrateXsection input_iXsection)    // constructor (also have consistency checks)
    {
        if (input_iPDF.param.pid1 != input_iXsection.param.pid1)
        {
            throw domain_error("Inconsistent initialization of struct ParamConvolute1D. Expect iPDF.param.pid1 == iXsection.param.pid1.");
        }
        if (input_iPDF.param.pid2 != input_iXsection.param.pid2)
        {
            throw domain_error("Inconsistent initialization of struct ParamConvolute1D. Expect iPDF.param.pid2 == iXsection.param.pid2.");
        }
        if (input_iPDF.param.s_parton != input_iXsection.param.s_parton)
        {
            throw domain_error("Inconsistent initialization of struct ParamConvolute1D. Expect iPDF.param.s_parton == iXsection.param.s_parton. Got iPDF.param.s_parton = " + to_string(input_iPDF.param.s_parton) + ", iXsection.param.s_parton = " + to_string(input_iXsection.param.s_parton));
        }
        if (input_iPDF.param.ScaleFactorization != input_iXsection.param.ScaleRenormalization)
        {
            throw domain_error("Inconsistent initialization of struct ParamConvolute1D. Expect iPDF.param.ScaleFactorization == iXsection.param.ScaleRenormalization. muF != muR has not been implemented.");
        }
        iPDF = input_iPDF;
        iXsection = input_iXsection;
    }
    void setPID(int pid1, int pid2)
    {
        iPDF.param.pid1 = pid1;
        iXsection.param.pid1 = pid1;
        iPDF.param.pid2 = pid2;
        iXsection.param.pid2 = pid2;
    }
};
struct ParamPDFProduct {
    LHAPDF::PDF* pdf;
    int pid1;
    int pid2;
    double ScaleFactorization; // GeV
    double s_parton; // GeV
    double s_hadron; // GeV
    double IntegrationRelativeError = NAN;    // if NAN, the default IntegrationRelativeError is used
    size_t NumberOfSubintervals = 0;    // if NAN, the default NumberOfSubintervals is used
};

// ============ Functions to convolute with PDF ==============================================

// Compute integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection
double PDF_times_xsection(IntegratePDF iPDF, IntegrateXsection iXsection);

// Integrate 1/x1*PDF(x1)*PDF(s_parton/(s_hadron*x1)) over x1 (from s_parton/s_hadron to 1)
double Integrated_x1_PDF_Product(gsl_integration_workspace* workspace, LHAPDF::PDF* pdf, int pid1, int pid2, double s_parton, double s_hadron, double ScaleFactorization, double &error, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);

// Convolute (differential) partonic xsection with PDF. Assume (differential) xsection depends on partonic s (s_parton) and is independent of momentum fraction, i.e. Lorentz invariant
double Convolute_PDF_dxsection_ds(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double s_parton, double costhetamin = -1, double costhetamax = 1, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);

// Convolute (differential) partonic xsection with PDF. use yt,pT,ytbar instead of x1,x2,costhetaCM
double Convolute_PDF_dxsection_dyt(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double yt, double pTmin = 0., double pTmax = 0., double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);

// Convolute (differential) partonic xsection with PDF. use yt,pT,ytbar instead of x1,x2,costhetaCM
// Assume ytbarmin = - ytbarmax
double Convolute_PDF_dxsection_dpT(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double yt, double pT, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000);

// Integrate over dhadronicxsection_dsparton over s_parton
double xsection_hadron(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double s_min = 0, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000, double inp_s_max = 0);


#endif


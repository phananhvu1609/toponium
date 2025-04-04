
/********************************************************************************
 *									  *										    *
 * 	Contains convolution procedures:                                            *
 * 	    - 1D convolution where xsection and pdf are separable					*
 * 																				*
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/

#ifndef guard_Convolution_h
#define guard_Convolution_h

#include <iostream>
#include <gsl/gsl_errno.h>
#include "Utilities.h"
#include "Computations.h"
using namespace std;


// Convolute (differential) partonic xsection with PDF. Assume (differential) xsection depends on partonic s (s_parton) and is independent of momentum fraction, i.e. Lorentz invariant
// Note: This function doesn't integrate differential partonic xsection to give total partonic xsection. s_hadron is also the maximum of s_parton.
double Convolute_PDF_1D(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double s_min = 0, double IntegrationRelativeError = 1e-7, size_t NumberOfSubintervals = 1000, double inp_s_max = 0)
{	
    double s_hadron = iPDF.param.s_hadron;
    double mt = iXsection.param.mt;
	if (s_hadron <= 4*mypow(mt,2))
	{
		return 0;
	}
	if (iXsection.name_f_change_variable == "rapidity")
	{
		throw invalid_argument("Convolute_PDF_1D can only work with Lorentz invariant kinematic varialbes. Rapidity is not Lorentz invariant.");
	}
	if (s_min < 4*mypow(mt,2))
	{
		s_min = 4*mypow(mt,2);
	}

	// Parameters
	double s_max = s_hadron;
	if (inp_s_max > 0)
	{
		s_max = min(s_hadron,inp_s_max);
		if (s_max < s_min)
		{
			return 0;
		}
	}
	
	struct ParamConvolute1D pConvolute(iPDF,iXsection);

	// Integrate
	double result;
	int run_status = gsl_integration_1D(workspace,f_PDF_xsection,s_min,s_max,&pConvolute,result,error,IntegrationRelativeError,NumberOfSubintervals);

	// Check convergence
	if (run_status != GSL_SUCCESS)
	{
		cout << "[WARNING]    Result of Convolute_PDF(pdf," << iXsection.param.SecondKinematicVariable << "," << iXsection.param.alphaS << "," << iXsection.param.mt << "," << s_hadron << "," << iXsection.param.pid1 << "," << iXsection.param.pid2 << ") = " << result << " +- " << error << " is unreliable. !!!" << endl;
	}
	return result;
}

#endif
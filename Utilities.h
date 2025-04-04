/********************************************************************************
 *									  *										    *
 * 	Contains ultility procedures:                                               *
 * 	    - Printing messanges (welcome, etc.)									*
 * 	    - New or alternative mathematical functions (pow, approx, etc.)			*
 * 																				*
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/

#ifndef guard_Utilities_h
#define guard_Utilities_h

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <complex>
using namespace std;

// ============ Some physical constants =============
#define invGeV2topb 0.389379e9

// ============ Some auxiliary functions ============

// Print a greeting message at startup
void Greetings();

// Print results. res[0] = res, res[1] = res_lower, res[2] = res_upper
template<typename Stream>
void Print_Result(Stream sout, double res[]);

// Replacement function for pow. Should be used for integer exponent < 100 (for u^n with integer n < 100)
template<typename Number>
Number mypow(Number x, unsigned int n);

template<typename Number, typename Integer>
Number mypow(Number x, Integer n);

// Approximately equal
template<typename Number, typename Real>
bool approx(Number x, Number y, Real rel_err = 1e-4);

// Round to certain digit, e.g. precision = 1e-2 if user want to round up to 2 decimal points
template<typename Number>
Number round(Number value, Number precision);


// ============ General integration procedure ============

// 1D adaptive integration using GSL
template<typename Function, typename Param>
int gsl_integration_1D(gsl_integration_workspace* w, Function f, double xmin, double xmax, Param param, double &result, double &error, double RelativeError = 1e-7, size_t NumberOfSubintervals = 1000, int routine = GSL_INTEG_GAUSS15)
{
	// Initialize integration
	gsl_function F;
	F.function = f;
	if (!(param == nullptr))
	{
		F.params = param;
	}
	

	// Do integration
	int run_status = gsl_integration_qag(&F, xmin, xmax, 0, RelativeError, NumberOfSubintervals, routine, w, &result, &error);

	// Check convergence
	if ((run_status != GSL_SUCCESS && run_status != GSL_EROUND) || w->size >= NumberOfSubintervals || (abs(error) > abs(RelativeError*result) && abs(result) > mypow(abs(min(min(1./xmin,1./xmax),min(xmin,xmax))),3) && abs(result) > 0.) || (abs(error) > abs(min(0.1,abs(20.*RelativeError))*result) && abs(result) > 0.) || (abs(error) > 1e-30 && result == 0))
	{
		cout << "!!! Warning. Integration (gsl_integration_1D) is unsuccessful (code " << run_status << "). ";
		if (w->size >= NumberOfSubintervals)
		{
			cout << "Maximum iteration reached. ";
		}
		if ((abs(error) > abs(RelativeError*result) && abs(result) > mypow(abs(min(min(1./xmin,1./xmax),min(xmin,xmax))),3) && abs(result) > 0.) || (abs(error) > abs(min(0.1,abs(20.*RelativeError))*result) && abs(result) > 0.) || (abs(error) > 1e-30 && result == 0))
		{
			cout << "Desired relative error " << RelativeError << " not reached. Actual relative error: " << abs(error/result) << ".";
		}
		
		cout << "!!!" << endl;
		if (run_status == GSL_SUCCESS || run_status == GSL_EROUND)
		{
			run_status = 1025;
		}
	}

	return run_status;
}
// 1D adaptive integration using GSL + workspace initialization (use this overload if there are parameters to pass to GSL)
template<typename Function, typename Param>
int gsl_integration_1D(Function f, double xmin, double xmax, Param param, double &result, double &error, double RelativeError = 1e-7, size_t NumberOfSubintervals = 1000, int routine = GSL_INTEG_GAUSS61);
// 1D adaptive integration using GSL (use this overload if there are NO parameters to pass to GSL)
template<typename Function>
int gsl_integration_1D(Function f, double xmin, double xmax, double &result, double &error, double RelativeError = 1e-7, size_t NumberOfSubintervals = 1000, int routine = GSL_INTEG_GAUSS61);


// ============ For testing ============

// Test function
double f_test (double, void * );

#endif

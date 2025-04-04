/********************************************************************************
 *									  *										    *
 * 	Contains ultility procedures (list can be found in Utilities.h)             *
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/


#include "Utilities.h"


/*******************************************************************************	
*									  *										   *
*								Print  messages		               	       	   *
*									  *										   *
*******************************************************************************/
void Greetings()
{
    cout << endl;
	cout << "......................................................................................................" << endl;
	cout << ".                                                                                                    ." << endl;
    cout << ".    Program \"Toponium\" starts                                                                       ." << endl;
    cout << ".    This program computes the hadronic observables by convoluting the partonic ones with the PDF    ." << endl;
	cout << ".                                                                                                    ." << endl;
	cout << "......................................................................................................" << endl;
    cout << endl;
}

/********************************************************************************	
*									  *											*
*							Mathmematical operations							*
*									  *											*
********************************************************************************/
// Replacement function for pow. Should be used for integer exponent < 100 (for u^n with integer n < 100)
template<typename Number>
Number mypow(Number x, unsigned int n) {
	/*
		Repeated squaring method. Returns 0.0^0 = 1.0, so continuous in x
	*/
  	Number value = 1.0;
	do {
		if(n & 1) value *= x;  /* for n odd */
		n >>= 1;
		x *= x;
	} while (n);

	return value;
}

template<typename Number, typename Integer>
Number mypow(Number x, Integer n) {
	unsigned int un;
	if (n == (int)n){
		if(n < 0) {
			x = 1.0/x;
			un = -n;
		} else {
			un = n;
		}
		return mypow(x,un);
	}
	else {
		return pow(x,n);
	}
}

// Approximately equal
template<typename Number, typename Real>
bool approx(Number x, Number y, Real rel_err) {
	if (abs(x-y) <= rel_err*min(abs(x),abs(y)))
	{
		return true;
	}
	else
	{
		return false;
	}
}

// Round to certain digit, e.g. precision = 1e-2 if user want to round up to 2 decimal points
template<typename Number>
Number round(Number value, Number precision)
{
    return round(value / precision) * precision;
}


/********************************************************************************	
*									  *											*
*							Integration operations								*
*									  *											*
********************************************************************************/

// 1D adaptive integration using GSL + workspace initialization
template<typename Function, typename Param>
int gsl_integration_1D(Function f, double xmin, double xmax, Param param, double &result, double &error, double RelativeError, size_t NumberOfSubintervals, int routine)
{
	// Initialize integration
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (NumberOfSubintervals);
	// Do integration
	int run_status = gsl_integration_1D(w,f,xmin,xmax,param,result,error,RelativeError,NumberOfSubintervals,routine);
	// Free up workspace
	gsl_integration_workspace_free (w);

	return run_status;
}

// 1D adaptive integration using GSL
template<typename Function>
int gsl_integration_1D(Function f, double xmin, double xmax, double &result, double &error, double RelativeError, size_t NumberOfSubintervals, int routine)
{
	return gsl_integration_1D(f,xmin,xmax,nullptr,result,error,RelativeError,NumberOfSubintervals,routine);
}


/*******************************************************************************	
*									  *										   *
*							    Test function		               	       	   *
*									  *										   *
*******************************************************************************/
double f_test (double x, void * params)
{
  double f = 1/x;
  return f;
}


/********************************************************************************	
*									  *											*
*					Instantiate all template instances							*
*									  *											*
********************************************************************************/
template double mypow<double>(double, unsigned int);
template double mypow<double,int>(double, int);
template double mypow<double,double>(double, double);
template complex<double> mypow<complex<double>,int>(complex<double>, int);
template bool approx<double,double>(double, double, double);
template bool approx<complex<double>,double>(complex<double>, complex<double>, double);
template double round<double>(double value, double precision);
template int gsl_integration_1D<double (*)(double, void*), double*>(double (*)(double, void*), double, double, double*, double&, double&, double, size_t, int);
template int gsl_integration_1D<double (*)(double, void*)>(double (*)(double, void*), double, double, double&, double&, double, size_t, int);
// template int gsl_integration_1D<double (double, void*)>(double (double, void*), double, double, double&, double&, double, size_t, int);
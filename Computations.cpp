/********************************************************************************
 *									  *										    *
 * 	Contains computational procedures (list can be found in Formulas.h)         *
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/


#include "Computations.h"


/********************************************************************************	
*									  *										    *
*		Compute xsection of g + g -> t + tbar at tree-level  (in GeV^-2)		*
*									  *									        *
********************************************************************************/
double xsection_ggttTree(gsl_integration_workspace* workspace, double alphaS, double mt, double s, double &error, double costhetamin, double costhetamax, double IntegrationRelativeError, size_t NumberOfSubintervals)
{
	if (s <= 4*mypow(mt,2))
	{
		return 0;
	}
	
	// Integrate
	double result;
	double param[] = {alphaS, mt, s};
	int run_status = gsl_integration_1D(workspace,f_dxsection_dcosthetaCM_ggttTree,costhetamin,costhetamax,param,result,error,IntegrationRelativeError,NumberOfSubintervals);

	// Check convergence
	if (run_status != GSL_SUCCESS)
	{
		cout << "[WARNING]    Result of xsection_ggttTree(" << alphaS << "," << mt << "," << s << ") = " << result << " +- " << error << " is unreliable. !!!" << endl;
	}
	return result;
}

double xsection_ggttTree_noerr(gsl_integration_workspace* workspace, double alphaS, double mt, double s, double costhetamin, double costhetamax, double IntegrationRelativeError, size_t NumberOfSubintervals)
{
	double error;
	return xsection_ggttTree(workspace, alphaS, mt, s, error, costhetamin, costhetamax, IntegrationRelativeError, NumberOfSubintervals);
}

double xsection_ggttTree_y(gsl_integration_workspace* workspace, double alphaS, double mt, double s, double &error, double ymin, double ymax, double IntegrationRelativeError, size_t NumberOfSubintervals)
{
	if (s <= 4*mypow(mt,2))
	{
		return 0;
	}
	if (isnan(ymin) || isnan(ymax) || ymin == HUGE_VALF || ymax == HUGE_VALF)
	{
		double ECM = sqrt(s)/2;
		double pCM = sqrt(s/4 - mypow(mt,2));
		ymax = abs(0.5*log((ECM + pCM) / (ECM - pCM)));
		ymin = -ymax;
	}
	
	
	// Integrate
	double result;
	double param[] = {alphaS, mt, s};
	int run_status = gsl_integration_1D(workspace,f_dxsection_dy_ggttTree,ymin,ymax,param,result,error,IntegrationRelativeError,NumberOfSubintervals);

	// Check convergence
	if (run_status != GSL_SUCCESS)
	{
		cout << "[WARNING]    Result of xsection_ggttTree_y(" << alphaS << "," << mt << "," << s << ") = " << result << " +- " << error << " is unreliable. !!!" << endl;
	}
	return result;
}


/********************************************************************************	
*									  *										    *
*				GSL conventions for functions to convolute with PDF				*
*									  *									        *
********************************************************************************/

// Compute integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection
double PDF_times_xsection(IntegratePDF iPDF, IntegrateXsection iXsection)
{
	auto name_f_xsection = iXsection.name_f_xsection;
	auto name_f_change_variable = iXsection.name_f_change_variable;
	auto paramPDF = iPDF.param;
	auto paramXsection = iXsection.param;
	double s_parton = paramXsection.s_parton;
	
	// Check consistency
	if (paramXsection.s_parton != paramPDF.s_parton)
	{
		throw invalid_argument("Inconsistent option for iPDF and iXsection (in PDF_times_xsection). Require iPDF.param.s_parton == iXsection.param.s_parton.");
	}
	if (paramXsection.pid1 != paramPDF.pid1)
	{
		throw invalid_argument("Inconsistent option for iPDF and iXsection (in PDF_times_xsection). Require iPDF.param.pid1 == iXsection.param.pid1.");
	}
	if (paramXsection.pid2 != paramPDF.pid2)
	{
		throw invalid_argument("Inconsistent option for iPDF and iXsection (in PDF_times_xsection). Require iPDF.param.pid2 == iXsection.param.pid2.");
	}

	// Check dynamic scale
	if (iXsection.DynamicScaleChoice > 0 && !paramXsection.is_diffxsection)
	{
		throw domain_error("Invalid option for iXsection.DynamicScaleChoice. Dynamic scaling requires numerically integrating over differential xsection.");
	}

	// Compute xsection or differential xsection with appropriate Jacobian
	double xsection = 0, error;
	if (!paramXsection.is_diffxsection)
	{
		bool isUseDefault = false;
		if (isnan(iXsection.MinSecondKinematicVariable) || isnan(iXsection.MaxSecondKinematicVariable) || isnan(iXsection.IntegrationRelativeError) || (iXsection.NumberOfSubintervals <= 0))
		{
			isUseDefault = true;
		}
		
		// Compute xsection
		if (name_f_xsection == "xsection_ggttTree")
		{
			if (isUseDefault)
			{
				xsection = xsection_ggttTree(iXsection.workspace,paramXsection.alphaS,paramXsection.mt,s_parton,error);
			}
			xsection = xsection_ggttTree(iXsection.workspace,paramXsection.alphaS,paramXsection.mt,s_parton,error,iXsection.MinSecondKinematicVariable,iXsection.MaxSecondKinematicVariable,iXsection.IntegrationRelativeError,iXsection.NumberOfSubintervals);
		}
		else if (name_f_xsection == "dxsect_dmToponium_partonic")
		{
			xsection = dxsect_dmToponium_partonic(paramXsection.pid1,paramXsection.pid2,s_parton,paramXsection.mtt,paramPDF.ScaleFactorization,paramXsection.ScaleRenormalization,paramXsection.mt,paramXsection.TopDecayWidth,paramXsection.BoundStateSpin,paramXsection.BoundStateJ,paramXsection.BoundStateColorConfig,iPDF.pdf);
		}
		else if (name_f_xsection == "")
		{
			xsection = 1;
		}
		else
		{
			throw domain_error("Unrecognized option \"" + name_f_xsection + "\" for name_f_xsection.");
		}
	}
	else
	{
		throw domain_error("Other differential xsection not implemented yet. Please use dsigma/dmtt instead.");
	}
	
	// Compute PDF
	double integratedPDFProduct = 1;
	if (name_f_change_variable == "rapidity")
	{
		throw domain_error("Unrecognized option \"" + name_f_change_variable + "\" for name_f_change_variable.");
	}
	else
	{
		if (isnan(iPDF.IntegrationRelativeError) || (iPDF.NumberOfSubintervals <= 0))
		{
			integratedPDFProduct = Integrated_x1_PDF_Product(iPDF.workspace,iPDF.pdf,paramPDF.pid1,paramPDF.pid2,s_parton,paramPDF.s_hadron,paramPDF.ScaleFactorization,error);
		}
		else
		{
			integratedPDFProduct = Integrated_x1_PDF_Product(iPDF.workspace,iPDF.pdf,paramPDF.pid1,paramPDF.pid2,s_parton,paramPDF.s_hadron,paramPDF.ScaleFactorization,error,iPDF.IntegrationRelativeError,iPDF.NumberOfSubintervals);
		}
	}
	
	return 1./paramPDF.s_hadron*integratedPDFProduct*xsection;
}

// GSL convention for Integrate over s_parton of integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection
// Warning: second integration variable must be independent of s_parton
double f_PDF_xsection(double s_parton, void *args)
{
	struct ParamConvolute1D* params = (ParamConvolute1D*)args;
	struct IntegratePDF iPDF = params->iPDF;
	struct IntegrateXsection iXsection = params->iXsection;
	iPDF.param.s_parton = s_parton;
	iXsection.param.s_parton = s_parton;

	return PDF_times_xsection(iPDF,iXsection);
}

// GSL convention for Integrate over costhetaCM of integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection
// Warning: second integration variable must be independent of s_parton
double f_PDF_xsection_costhetaCM(double costhetaCM, void *args)
{
	struct ParamConvolute1D* params = (ParamConvolute1D*)args;
	struct IntegratePDF iPDF = params->iPDF;
	struct IntegrateXsection iXsection = params->iXsection;
	iXsection.param.SecondKinematicVariable = costhetaCM;
	iXsection.name_f_change_variable = "";

	return PDF_times_xsection(iPDF,iXsection);
}

// GSL convention for Integrate over s_parton of integrated( integrated(1/(s*x1)*PDF(x1)*PDF(s_parton/(s_hadron*x1)),x1) * xsection, costhetaCM)
// Warning: second integration variable must be independent of s_parton
double f_dhadronicxsection_dsparton(double s_parton, void *args)
{
	struct ParamConvolute1D* params = (ParamConvolute1D*)args;
	struct IntegratePDF iPDF = params->iPDF;
	struct IntegrateXsection iXsection = params->iXsection;
	iXsection.param.SecondKinematicVariable = NAN;
	iXsection.name_f_change_variable = "";
	iXsection.param.s_parton = s_parton;
	iPDF.param.s_parton = s_parton;
	double costhetamin = iXsection.MinSecondKinematicVariable;
	double costhetamax = iXsection.MaxSecondKinematicVariable;
	double error;

	return Convolute_PDF_dxsection_ds(iPDF,iXsection,error,iXsection.workspace2,s_parton,costhetamin,costhetamax,iXsection.IntegrationRelativeError,iXsection.NumberOfSubintervals);
}

// GSL conventions for 1/x1*PDF(x1)*PDF(s_parton/(s_hadron*x1))
double f_PDF_product(double x1, void *args)
{
	struct ParamPDFProduct* params = (ParamPDFProduct*)args;
	double Q2 = mypow(params->ScaleFactorization,2);
	double x2 = (params->s_parton)/((params->s_hadron) * x1);
	LHAPDF::PDF* pdf = params->pdf;

	// Test physical conditions
	if (!( pdf->inPhysicalRangeX(x1) && pdf->inPhysicalRangeX(x2) && pdf->inPhysicalRangeQ2(Q2) ))
	{
		cout << "!!! Warning. Unphysical (x1,x2,Q2) = (" << x1 << "," << x2 << "," << Q2 << ") in f_PDF_product. Return 0. !!!" << endl;
		return 0;
	}

	// Compute
	double f1 = (pdf->xfxQ2(params->pid1,x1,Q2)) / x1;
	double f2 = (pdf->xfxQ2(params->pid2,x2,Q2)) / x2;
	if (f1 < 0)
	{
		f1 = 0;
	}
	if (f2 < 0)
	{
		f2 = 0;
	}
	
	
	return 1./x1 * f1 * f2;
}

// GSL conventions for 1/x1*PDF(x1)*PDF(s_parton/(s_hadron*x1)), u = ln(x)
double f_PDF_product_ln(double u, void *args)
{
	double x1 = exp(u), Jacobian = x1;
	return f_PDF_product(x1, args) * Jacobian;
}



/********************************************************************************	
*									  *										    *
*						Convolute a function with PDF  			    			*
*									  *									        *
********************************************************************************/

// Integrate 1/x1*PDF(x1)*PDF(s_parton/(s_hadron*x1)) over x1 (from s_parton/s_hadron to 1). Integrate over u = ln(x)
double Integrated_x1_PDF_Product(gsl_integration_workspace* workspace, LHAPDF::PDF* pdf, int pid1, int pid2, double s_parton, double s_hadron, double ScaleFactorization, double &error, double IntegrationRelativeError, size_t NumberOfSubintervals)
{
	// Parameters
	double xmin = s_parton / s_hadron;
	double xmax = pdf->xMax();
	double umin = log(xmin), umax = log(xmax);
	struct ParamPDFProduct params = {pdf,pid1,pid2,ScaleFactorization,s_parton,s_hadron};

	// Integrate
	double result;
	int run_status = gsl_integration_1D(workspace,f_PDF_product_ln,umin,umax,&params,result,error,IntegrationRelativeError,NumberOfSubintervals);

	// Check convergence
	if (run_status != GSL_SUCCESS)
	{
		cout << "[WARNING]    Result of Integrated_x1_PDF_Product(pdf," << pid1 << "," << pid2 << "," << s_parton << "," << s_hadron << "," << ScaleFactorization << ") = " << result << " +- " << error << " is unreliable. !!!" << endl;
	}
	return result;
}


// Convolute (differential) partonic xsection with PDF. Assume (differential) xsection depends on partonic s (s_parton) and is independent of momentum fraction, i.e. Lorentz invariant
double Convolute_PDF_dxsection_ds(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double s_parton, double costhetamin, double costhetamax, double IntegrationRelativeError, size_t NumberOfSubintervals)
{	
    double s_hadron = iPDF.param.s_hadron;
    double mt = iXsection.param.mt;
	if (s_hadron <= 4*mypow(mt,2))
	{
		return 0;
	}

	// parameters
	iXsection.param.s_parton = s_parton;
	iPDF.param.s_parton = s_parton;
	struct ParamConvolute1D pConvolute(iPDF,iXsection);

	// Integrate
	double result;
	int run_status = gsl_integration_1D(workspace,f_PDF_xsection_costhetaCM,costhetamin,costhetamax,&pConvolute,result,error,IntegrationRelativeError,NumberOfSubintervals);

	// Check convergence
	if (run_status != GSL_SUCCESS)
	{
		cout << "[WARNING]    Result of Convolute_PDF_dxsection_ds(pdf," << costhetamin << "," << costhetamax << "," << iXsection.param.alphaS << "," << iXsection.param.mt << "," << s_hadron << "," << iXsection.param.pid1 << "," << iXsection.param.pid2 << ") = " << result << " +- " << error << " is unreliable. !!!" << endl;
	}
	return result;
}


// Integrate over dhadronicxsection_dsparton over s_parton
double xsection_hadron(IntegratePDF iPDF, IntegrateXsection iXsection, double &error, gsl_integration_workspace* workspace, double s_min, double IntegrationRelativeError, size_t NumberOfSubintervals, double inp_s_max)
{	
    double s_hadron = iPDF.param.s_hadron;
    double mt = iXsection.param.mt;
	if (s_hadron <= 4*mypow(mt,2))
	{
		return 0;
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
	if (!iXsection.workspace2)
	{
		throw invalid_argument("Invalid iXsection.workspace 2 (in xsection_hadron). Need workspace to do integration.");
	}
	
	struct ParamConvolute1D pConvolute(iPDF,iXsection);

	// Integrate
	double result;
	int run_status = gsl_integration_1D(workspace,f_dhadronicxsection_dsparton,s_min,s_max,&pConvolute,result,error,IntegrationRelativeError,NumberOfSubintervals);

	// Check convergence
	if (run_status != GSL_SUCCESS)
	{
		cout << "[WARNING]    Result of xsection_hadron(pdf," << iXsection.param.SecondKinematicVariable << "," << iXsection.param.alphaS << "," << iXsection.param.mt << "," << s_hadron << "," << iXsection.param.pid1 << "," << iXsection.param.pid2 << ") = " << result << " +- " << error << " is unreliable. !!!" << endl;
	}
	return result;
}

/********************************************************************************	
*									  *										    *
*		Rewrite formulas (in Formulas.cpp) in a form readable by GSL  			*
*									  *									        *
********************************************************************************/
// GSL convention for dxsection_dcosthetaCM_ggttTree. params[] = {alphaS, mt, s}
double f_dxsection_dcosthetaCM_ggttTree(double costhetaCM, void * params)
{
	double* param = (double*)params;
	double alphaS = param[0], mt = param[1], s = param[2];
	return dxsection_dcosthetaCM_ggttTree(alphaS, mt, s, costhetaCM);
}
// GSL convention for dxsection_dy_ggttTree. params[] = {alphaS, mt, s}
double f_dxsection_dy_ggttTree(double y, void * params)
{
	double* param = (double*)params;
	double alphaS = param[0], mt = param[1], s = param[2];
	return dxsection_dy_ggttTree(alphaS, mt, s, y);
}
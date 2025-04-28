/********************************************************************************
 *									  *										    *
 * 	Program to convolute the partonic (differential) cross section with PDF to 	*
 * 	get the hadronic counterpart												*
 * 																				*
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/

#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <fstream>
#include <map>
#include <stdexcept>
#include <chrono>
#include <math.h>
#include <complex>
#include <filesystem>
#include <typeinfo>
#include <cstdlib>

#include "Computations.h"
#include "Convolution.h"
#include "Utilities.h"

using namespace LHAPDF;
using namespace std;


int main() {

	auto t_start = chrono::high_resolution_clock::now(), t_end = t_start;
	auto t_programstart = chrono::high_resolution_clock::now(), t_programend = t_programstart; 
	double dif,programdif;
	Greetings();

/*******************************************************************************	
*									  *										   *
*						INITIALIZE INPUT PARAMETERS		               	       *
*									  *										   *
*******************************************************************************/

	map<string, vector<string>> options;
	const int NumberOfParameters = 21; //Number of input parameters. ! Modify if options are added or deleted !
	cout << ".................. Initialization starts .................." << endl;


	//=========================== Default values for all options ===========================

	// Default values
	const string defaults[NumberOfParameters*2] = {
		// General setup (pdfset, pdfmember,order)
		"PDFSetName", "NNPDF31_nnlo_as_0118_notop",
		"PDFmember", "0",
		"Order", "NLO",
		// Output configurations
		"OutputFolder", "output",
		"mttbarmin", "335",
		"mttbarmax", "380",
		"mttbarstep", "0.1",
		// mt(GeV); TopDecayWidth(GeV); muF and muR (in units of mt), muF goes into PDF, muR goes into loop integration. Default value taken in accordance with arXiv:1701.06228
		"Mtop", "172.5",
		"TopDecayWidth", "1.326",
		"muF", "1.0",
		"muR", "1.0",
		// Other SM physics parameters (energy in GeV, alphaS at MZ)
		"EbeamCM", "6500",
		"alphaS", "0.118",
		"MZ", "91.2",
		"vev", "246.22",
		// Options for integration
		"RelativeError", "1e-5",
		"MaxNumberSubintervals", "1000",
		// Other configurations (debug = 0 (no test), 1 (test, no debug), 2 (print debug); isComputeTotXsection = 0 (don't compute), 1 (compute))
		"debug", "1",
		"InternalRelativeError", "1e-5",
		"DynamicScalingChoice", "2",
		"ScaleUncertainty", "0.5	2",
	};
	
	// Save defaults values to 'options'
	for (int i = 0; i < NumberOfParameters; i++)
	{
		string key = defaults[2*i];
		string buffer = defaults[2*i+1];	// options[key] = buffer
		stringstream line(buffer);
		string val;
		vector<string> values;
		while (line)	// read line word for word. Convert line to a vector of words.
		{
			val.clear();
			line >> val;
			if (val.empty()) continue;
			values.push_back(val);
		}
		options[key] = values;		
	}


	// ============ Override default options by those provided by users in hadronic.inp ============
	// Every option is a string at this point. 
	// Comments start with a "#"
	
	// Read parameters from input file
	ifstream input("hadronic.inp");
	while (input)
	{
		string buffer;
		getline(input, buffer);	// read 'input' line by line
		if (buffer.empty()) continue;
		
		stringstream line(buffer);
		string key, val;
		vector<string> values;
		line >> key;
		while (line)	// read line word for word. Convert line to a vector of words. Skip lines starting with '#'
		{
			val.clear();
			line >> val;
			if (val.empty()) continue;
			if (val[0] == '#') break;
			values.push_back(val);
		}
		options[key] = values;
	}
	
	// assignment input parameters to option
	string PDFSetName = options["PDFSetName"][0];
	int PDFmember = atoi(options["PDFmember"][0].c_str());
	string Order = options["Order"][0];
	string OutputFolder = options["OutputFolder"][0];
	double mttbarmin = atof(options["mttbarmin"][0].c_str());
	double mttbarmax = atof(options["mttbarmax"][0].c_str());
	double mttbarstep = atof(options["mttbarstep"][0].c_str());
	double Mtop = atof(options["Mtop"][0].c_str());
	double TopDecayWidth = atof(options["TopDecayWidth"][0].c_str());
	double muF = atof(options["muF"][0].c_str());
	double muR = atof(options["muR"][0].c_str());
	double EbeamCM = atof(options["EbeamCM"][0].c_str());
	double alphaS = atof(options["alphaS"][0].c_str());
	double MZ = atof(options["MZ"][0].c_str());
	double vev = atof(options["vev"][0].c_str());
	double RelativeError = atof(options["RelativeError"][0].c_str());
	size_t MaxNumberSubintervals = atoi(options["MaxNumberSubintervals"][0].c_str());
	int debug = atoi(options["debug"][0].c_str());
	double InternalRelativeError = atof(options["InternalRelativeError"][0].c_str());
	int DynamicScalingChoice = atoi(options["DynamicScalingChoice"][0].c_str());
	double ScaleUncertainty[2] = {atof(options["ScaleUncertainty"][0].c_str()), atof(options["ScaleUncertainty"][1].c_str())};


	// ============ Check if the parameters' values are allowed ============

	if (PDFmember < 0)
	{
		throw domain_error("Incorrect option PDFmember: must be an integer >= 0.");
	}
	if (debug < 0)
	{
		throw domain_error("Incorrect option for debug: must be an integer >= 0.");
	}
	if (DynamicScalingChoice != 0)
	{
		throw domain_error("Incorrect option for DynamicScalingChoice: only fixed scaling has been implemented.");
	}
	if (Mtop < 0)
	{
		throw domain_error("Incorrect option Mtop: mass must be >= 0.");
	}
	if (EbeamCM < 0)
	{
		throw domain_error("Incorrect option EbeamCM: beam energy must be >= 0.");
	}
	if (Order != "NLO")
	{
		throw domain_error("Incorrect option Order: Only NLO has been implemented.");
	}
	if (alphaS < 0)
	{
		throw domain_error("Incorrect option alphaS: alphaS must be >= 0.");
	}
	if (alphaS > sqrt(4*M_PI))
	{
		throw domain_error("Incorrect option alphaS: alphaS too large for perturbative computation.");
	}
	if (!approx(muF,muR,InternalRelativeError))
	{
		throw domain_error("Incorrect option for muF, muR: expect muF = muR. Separation of muF and muR has not been implemented.");
	}
	if ((PDFSetName == "NNPDF31_nnlo_as_0118_notop") && !approx(alphaS,0.118,InternalRelativeError))
	{
		cout << "[WARNING] Input PDFSet and alphaS are incompatible. NNPDF31_nnlo_as_0118_notop should work with alphaS(MZ) = 0.118." << endl;
	}
	if (mttbarmin < 0.7*2*Mtop || mttbarmax > 1.3*2*Mtop)
	{
		cout << "[WARNING] mttbar provided is far from ttbar threshold. This code assumes non-relativistic tops. Results may be unreliable." << endl;
	}
	

/*******************************************************************************	
*									  *										   *
*						INITIALIZE PDFs AND VARIABLES		               	   *
*									  *										   *
*******************************************************************************/

	// ============ Initialization ============

	// Initialize PDF	
	PDF* pdf = mkPDF(PDFSetName, PDFmember);

	// Extract info from PDF
	vector<int> pids = pdf->flavors();

	// Create save folder (if not already exist)
	filesystem::create_directory(OutputFolder);

	// char* that store today in the form YYMMDD
	char datestring [12];
	struct tm * timeinfo;
	time_t rawtime;
	time (&rawtime);
	timeinfo = localtime(&rawtime);
	strftime (datestring,12,"%y%m%d-%H%M",timeinfo);
	filesystem::create_directory(OutputFolder + "/" + datestring);

	// Prevent GSL from aborting the program
	gsl_set_error_handler_off();

	// Initialize integration workspace
	gsl_integration_workspace * w1 = gsl_integration_workspace_alloc (MaxNumberSubintervals);
	gsl_integration_workspace * w2 = gsl_integration_workspace_alloc (MaxNumberSubintervals);
	gsl_integration_workspace * w3 = gsl_integration_workspace_alloc (MaxNumberSubintervals);

	// Other physical parameters
	double ScaleFactorization = muF*Mtop;
	double ScaleRenormalization = muR*Mtop;
	double mt = Mtop;
	double s_hadron = mypow(EbeamCM*2,2);
	double s_min = mypow(mttbarmin,2);
	int pidg = 21;
	const size_t pidq_length = 5, pidquark_length = 10;
	int pidq[] = {1,2,3,4,5,6}, pidquark[pidquark_length] = {-5,-4,-3,-2,-1,1,2,3,4,5};
	double alphaS_muR = pdf->alphasQ(ScaleRenormalization);
	cout << "[INFO] alphaS(muR = " << ScaleRenormalization << " GeV) = " << alphaS_muR << " shall be used for fixed scale calculations, unless otherwise stated." << endl;
	double ScaleUncertaintyFactor[3] = {1,ScaleUncertainty[0],ScaleUncertainty[1]};


	// ============ Consistency checks ============
	if ( !approx(alphaS,pdf->alphasQ(MZ),InternalRelativeError) )
	{
		cout << "[WARNING] Input PDFSet and alphaS are incompatible (input " << alphaS << " vs PDF's " << pdf->alphasQ(MZ) << "). Please note that input alphaS shouble be taken at MZ." << endl;
	}
	if ( !(pdf->inPhysicalRangeQ(muF*Mtop)) )
	{
		throw domain_error("Incorrect option muF: muF outside of range of PDF.");
	}
	if ( !(pdf->inPhysicalRangeQ(muR*Mtop)) )
	{
		cout << "[WARNING] muR outside of range of PDF. Result can be unreliable" << endl;
	}
	

	// Initialization completed
	t_end = chrono::high_resolution_clock::now();
	dif = std::chrono::duration<double, std::milli>(t_end-t_programstart).count() / 1000.;
	cout << "[INFO] Initialization finshied (it took " << dif << " s)" << endl << endl;
	

/*******************************************************************************	
*									  *										   *
*					SOME TESTS, check if code run correctly		               *
*									  *										   *
*******************************************************************************/

	if (debug >= 1)
	{	
		// Initialize testing
		filesystem::create_directory(OutputFolder + "/test");
		filesystem::create_directory(OutputFolder + "/test/" + datestring);
		t_start = chrono::high_resolution_clock::now();
		cout << ".................. Tests start .................." << endl;
		double alphaS_muR = 0.107655, mt = 172.5, ScaleFactorization = mt, ScaleRenormalization = mt, costhetamin = -1, costhetamax = 1, EbeamCM = 6500, s_hadron = mypow(2.*EbeamCM,2), s = mypow(400.,2), vev = 246.22, MZ = 91.2;
		double TopDecayWidth = 1.326;
		int count_error = 0;
		cout << "   Default parameters: alphaS(muR) = " << alphaS_muR << ", mt = " << mt << ", muF = " << ScaleFactorization << ", muR = " << ScaleRenormalization << ", costhetamin = " << costhetamin << ", costhetamax = " << costhetamax << ", EbeamCM = " << EbeamCM << "," << endl;
		cout << "                       s_parton = (" << sqrt(s) << " GeV)^2, MZ = " << MZ << " GeV, vev = " << vev << " GeV." << endl;


		// ============ reproduce xf(x,mu) plot in fig 3.1, 1706.00428 ============
		{
			cout << " . Reproducing fig 3.1 in 1706.00428 (.dat files in " << OutputFolder << "/test)" << endl;
			const double q2_list[2] = {10,10000};
			const double DX = 0.01;
			const double MAXLOGX = 0;
			const double MINLOGX = -3;
			const int NX = (int) floor((MAXLOGX - MINLOGX)/DX) + 1;
			for (int pid : pids)
			{
				for (double q2 : q2_list)
				{
					const string spid = lexical_cast<string>(pid);
					const string sq2 = lexical_cast<string>(q2);
					const string smem = lexical_cast<string>(PDFmember);
					filesystem::create_directory(OutputFolder + "/test/q2-" + sq2);
					const string filename = OutputFolder + "/test/" + "q2-" + sq2 + "/"+ PDFSetName + "_" + smem + "_" + spid + ".dat";
					
					ofstream f(filename.c_str());
					f << "# x \t q^2 \t x*f" << endl;
					for (int ix = 0; ix < NX; ix++) {
						const double log10x = (MINLOGX + ix*DX < -1e-3) ? MINLOGX + ix*DX : 0;	// if log10x is > -1e-3, set log10x = 0
						const double x = pow(10, log10x);
						const double xf = pdf->xfxQ2(pid, x, q2);
						f << x << " " << q2 << " " << xf << endl;
					}
					f.close();
				}
			}
		}
		// Check gluon PDF
		{
			const double q2 = mypow(mt,2);
			const double DX = 0.01;
			const double MAXLOGX = 0;
			const double MINLOGX = -5;
			const int NX = (int) floor((MAXLOGX - MINLOGX)/DX) + 1;
			int pid = 21;
			const string spid = lexical_cast<string>(pid);
			const string sq2 = lexical_cast<string>(q2);
			const string smem = lexical_cast<string>(PDFmember);
			filesystem::create_directory(OutputFolder + "/test/q2-" + sq2);
			const string filename = OutputFolder + "/test/" + "q2-" + sq2 + "/"+ PDFSetName + "_" + smem + "_" + spid + ".dat";
			const string filename2 = OutputFolder + "/test/" + "q2-" + sq2 + "/integrated_PDF_"+ PDFSetName + "_" + smem + "_" + spid + ".dat";
			
			ofstream f(filename.c_str());
			f << "# x \t f" << endl;
			for (int ix = 0; ix < NX; ix++) {
				const double log10x = (MINLOGX + ix*DX < -1e-3) ? MINLOGX + ix*DX : 0;	// if log10x is > -1e-3, set log10x = 0
				const double x = pow(10, log10x);
				const double xf = pdf->xfxQ2(pid, x, q2);
				f << x << " " << xf/x << endl;
			}
			f.close();
			ofstream f2(filename2.c_str());
			f2 << "# s_parton \t integrated f" << endl;
			for (double mtt = 300.; mtt <= sqrt(s_hadron); mtt += (mtt < 400. ? 1e-2 : mtt < 1000 ? 1e-1 : 1)) {
				double s_parton = mypow(mtt,2);
				double error;
				double result = Integrated_x1_PDF_Product(w2,pdf,pid,pid,s_parton,s_hadron,sqrt(q2),error,RelativeError,MaxNumberSubintervals);
				f2 << s_parton << " " << result << endl;
			}
			f2.close();
		}


		// ============ compute tree-level total xsection and compare with top++ ============

		// Test gsl on known function
		{
			double result, error;
			double xmin = 1e-5, xmax = 1;
			cout << " . Integratiing 1/x using GSL integrator for x = " << xmin << " to " << xmax << ". " << flush;
			double expected = -log(xmin);
			gsl_integration_1D(f_test,xmin,xmax,result,error);
			cout << "Result = " << result << " +- " << error << ". Expect " << expected << "." << endl;
			if (!approx(result,expected,max(1e-4,max(InternalRelativeError,RelativeError))))
			{
				cout << "   !! Test failed !!" << endl;
				count_error++;
			}
		}

		// Compute xsection for gg > ttbar at tree level at s = (400 GeV)^2. alphaS is taken at MZ.
		{
			double expected = 3.3512039612527216e-8; // GeV^-2
			double error;
			double alphaS = 0.118;
			double result = xsection_ggttTree(w1,alphaS,mt,s,error,costhetamin,costhetamax,RelativeError,MaxNumberSubintervals);
			streamsize p = cout.precision();
			cout << " . Computing partonic xsection for tree gg > ttbar at s = (" << sqrt(s) << " GeV)^2, alphaS = " << alphaS << flush << ". Result = " << result;
			cout << " +- " << setprecision(2) << error << setprecision(p);
			cout << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
			double result_y = xsection_ggttTree_y(w1,alphaS,mt,s,error,HUGE_VALF,HUGE_VALF,RelativeError,MaxNumberSubintervals);
			cout << "   using rapidity: " << flush << "result = " << result_y << " GeV^-2." << endl;
			if (!approx(result,expected,max(1e-4,max(InternalRelativeError,RelativeError))) || !approx(result_y,expected,max(1e-4,max(InternalRelativeError,RelativeError))))
			{
				cout << "   !! Test failed !!" << endl;
				count_error++;
			}
		}
		// Convolute with PDF to get total xsection
		{
			double result, error, expect = 419.821;
			double xmin = 4*mypow(mt,2);
			int pid = 21;
			struct ParamPDF pPDF = {pid,pid,ScaleFactorization,NAN,s_hadron};
			struct IntegratePDF iPDF = {pdf,pPDF,w1,RelativeError,MaxNumberSubintervals};
			struct ParamXsection pXsection = {pid,pid,NAN,alphaS_muR,mt,ScaleRenormalization};
			struct IntegrateXsection iXsection = {"xsection_ggttTree",pXsection,w1,costhetamin,costhetamax,"",RelativeError,MaxNumberSubintervals};
			result = Convolute_PDF_1D(iPDF,iXsection,error,w2,xmin,RelativeError,MaxNumberSubintervals)*invGeV2topb;
			streamsize p = cout.precision();
			cout << " . Computing hadronic xsection for tree gg > ttbar for Ebeam = " << EbeamCM/1000 << " TeV. Result = " << result;
			cout << " +- " << setprecision(2) << error*invGeV2topb << setprecision(p);
			double expected = 411;
			cout << " pb. Expect " << expect << " (top++)." << endl;
			if (!approx(result,expected,0.1))
			{
				cout << "   !! Test failed !!" << endl;
				count_error++;
			}
		}
		// Convolute luminosity with hard function. Check with table 2 in 0812.0919
		{
			double mttbar = 2*mt, ScaleRenormalization = mt, ScaleFactorization = mt;
			double CF = 4./3.;	// Casimir operator for fundamental representation
			double CA = 3.;	// Casimir operator for adjoint representation
			double SoftScaleRenormalization = mt*CF*alphaS;
			const size_t ScaleUncertainty_count = 3;
			double ScaleUncertaintyFactor[ScaleUncertainty_count] = {1,2,4};
			double alphaS_Soft = pdf->alphasQ(SoftScaleRenormalization);
			double ImGFactor_octet = ImGreenFunction(mttbar,mt,TopDecayWidth,alphaS_Soft, CF, CA, 8) * invGeV2topb / mypow(mt,2) / mttbar;
			double ImGFactor_singlet = ImGreenFunction(mttbar,mt,TopDecayWidth,alphaS_Soft, CF, CA, 1) * invGeV2topb / mypow(mt,2) / mttbar;
			struct ParamPDF pPDF = {pidg,pidg,ScaleFactorization,NAN,s_hadron};
			struct IntegratePDF iPDF = {pdf,pPDF,w1,RelativeError,MaxNumberSubintervals};
			struct ParamXsection pXsection = {pidg,pidg,NAN,alphaS_muR,mt,ScaleRenormalization,false,NAN,MZ,vev};
			struct IntegrateXsection iXsection = {"",pXsection,w1,NAN,NAN,"",RelativeError,MaxNumberSubintervals};
			iPDF.InternalError = InternalRelativeError;
			iXsection.workspace2 = w2;
			struct ParamConvolute1D pConvolute(iPDF,iXsection);
			pConvolute.set_TopDecayWidth(TopDecayWidth);
			// gg > 1S_0^[1]
			pConvolute.set_Scale(ScaleRenormalization);
			pConvolute.setPID(pidg,pidg);
			pConvolute.set_boundstate(0,0,1);
			{
				double result[ScaleUncertainty_count];
				double expected[ScaleUncertainty_count] = {18.2e-6, 18.7e-6, 18.3e-6};
				cout << " . Computing L x F at mttbar = " << mttbar << " GeV for gg > 1S_0^[1]. (scale:result; expected) GeV^-2 = " << flush;
				for (size_t i = 0; i < ScaleUncertainty_count; i++)
				{
					double Scale = ScaleUncertaintyFactor[i]*mt;
					pConvolute.set_Scale(Scale);
					void* argConvolute = &pConvolute;
					result[i] = dsigmadmttbar(mttbar,argConvolute) / ImGFactor_singlet;
					cout << "(" << result[i] << ";" << expected[i] << ")  " << flush;
				}
				cout << endl;
				for (size_t i = 0; i < ScaleUncertainty_count; i++)
				{
					if (!approx(result[i],expected[i],0.1))
					{
						cout << "   !! Test failed !!" << endl;
						count_error++;
						break;
					}
				}
			}
			// gg > 1S_0^[8]
			pConvolute.set_Scale(ScaleRenormalization);
			pConvolute.setPID(pidg,pidg);
			pConvolute.set_boundstate(0,0,8);
			{
				void* argConvolute = &pConvolute;
				double result = dsigmadmttbar(mttbar,argConvolute) / ImGFactor_octet, expected = 55.8e-6;
				cout << " .                                         gg > 1S_0^[8]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.1))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// gg > 3S_1^[1]
			pConvolute.setPID(pidg,pidg);
			pConvolute.set_boundstate(1,1,1);
			{
				void* argConvolute = &pConvolute;
				double result = dsigmadmttbar(mttbar,argConvolute) / ImGFactor_singlet, expected = 0.175e-6;
				cout << " .                                         gg > 3S_1^[1]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.35))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// gg > 3S_1^[8]
			pConvolute.setPID(pidg,pidg);
			pConvolute.set_boundstate(1,1,8);
			{
				void* argConvolute = &pConvolute;
				double result = dsigmadmttbar(mttbar,argConvolute) / ImGFactor_octet, expected = 6.06e-6;
				cout << " .                                         gg > 3S_1^[8]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.35))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// qq > 1S_0^[1]
			{
				double result = 0., expected = 0.00664e-6;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidquark[i],-pidquark[i]);
					pConvolute.set_boundstate(0,0,1);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute) / ImGFactor_singlet;
					}
				}
				cout << " .                                         qq > 1S_0^[1]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.1))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// qq > 1S_0^[8]
			{
				double result = 0., expected = 0.0166e-6;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidquark[i],-pidquark[i]);
					pConvolute.set_boundstate(0,0,8);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute) / ImGFactor_octet;
					}
				}
				cout << " .                                         qq > 1S_0^[8]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.1))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// qq > 3S_1^[1]
			{
				double result = 0., expected = 0.;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidquark[i],-pidquark[i]);
					pConvolute.set_boundstate(1,1,1);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute) / ImGFactor_singlet;
					}
				}
				cout << " .                                         qq > 3S_1^[1]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.2))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// qq > 3S_1^[8]
			{
				double result = 0., expected = 21.7e-6;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidquark[i],-pidquark[i]);
					pConvolute.set_boundstate(1,1,8);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute) / ImGFactor_octet;
					}
				}
				cout << " .                                         qq > 3S_1^[8]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.1))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// gq > 1S_0^[1]
			{
				double result = 0., expected = -0.795e-6;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidg,pidquark[i]);
					pConvolute.set_boundstate(0,0,1);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute);
					}
					pConvolute.setPID(pidquark[i],pidg);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute);
					}
				}
				result /= ImGFactor_singlet;
				cout << " .                                         gq > 1S_0^[1]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.1))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// gq > 1S_0^[8]
			{
				double result = 0., expected = -1.99e-6;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidg,pidquark[i]);
					pConvolute.set_boundstate(0,0,8);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute);
					}
					pConvolute.setPID(pidquark[i],pidg);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute);
					}
				}
				result /= ImGFactor_octet;
				cout << " .                                         gq > 1S_0^[8]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.2))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			// gq > 3S_1^[8]
			{
				double result = 0., expected = 3.99e-6;
				for (size_t i = 0; i < pidquark_length; i++)
				{
					pConvolute.setPID(pidg,pidquark[i]);
					pConvolute.set_boundstate(1,1,8);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute);
					}
					pConvolute.setPID(pidquark[i],pidg);
					{
						void* argConvolute = &pConvolute;
						result += dsigmadmttbar(mttbar,argConvolute);
					}
				}
				result /= ImGFactor_octet;
				cout << " .                                         gq > 3S_1^[8]. Result = " << result << " GeV^-2. Expect " << expected << " GeV^-2." << endl;
				if (!approx(result,expected,0.2))
				{
					cout << "   !! Test failed !!" << endl;
					count_error++;
				}
			}
			
		}

		// Testing completed
		if (count_error > 0)
		{
			cout << "              vvvvvvvvvvvvvvvvv" << endl;
			cout << "[WARNING] >>> " << count_error << " test(s) failed. <<<" << endl;
			cout << "              ^^^^^^^^^^^^^^^^^" << endl;
		}
		else
		{
			cout << "[INFO] All test cases passed" << endl;
		}
		t_end = chrono::high_resolution_clock::now();
		dif = std::chrono::duration<double, std::milli>(t_end - t_start).count() / 1000.;
		cout << "[INFO] Tests finshied (it took " << dif << " s) <" << endl << endl;
	}

	return 0;
/*******************************************************************************	
*									  *										   *
*					COMPUTE TOPONIUM DIFF CROSS SECTION  		               *
*									  *										   *
*******************************************************************************/

	// Initialize ALP computation
	t_start = chrono::high_resolution_clock::now();
	cout << "............................................................." << endl;
	cout << "..................... Main computation ......................" << endl;
	const string result_filename = OutputFolder + "/" + datestring + "/result.dat";
	ofstream f_result(result_filename.c_str());	
	f_result << "# Parameters: alphaS(ScaleR) = " << alphaS_muR << ", m_t = " << mt << ", PDFset = " << PDFSetName << ", PDFmember = " << PDFmember << ", ScaleF = " << ScaleFactorization << ", ScaleR = " << ScaleRenormalization << ", sqrt(s_hadron) = " << sqrt(s_hadron) << ", MZ = " << MZ << ", vev = " << vev << ", RelativeError = " << RelativeError << ", ScaleUncertainty = (" << ScaleUncertainty[0] << "," << ScaleUncertainty[1] << ")" << endl;
	f_result << "# Results of toponium effect to hadronic xsection" << endl << endl;
	// double result_mtt_noCgg[3] = {0,0,0}, result_costhetaCM_noCgg[3] = {0,0,0}, result_y_noCgg[3] = {0,0,0}, result_pT_noCgg[3] = {0,0,0}, result_tot_noCgg[3] = {0,0,0};
	// double result_mtt_Cgg[3] = {0,0,0}, result_costhetaCM_Cgg[3] = {0,0,0}, result_y_Cgg[3] = {0,0,0}, result_pT_Cgg[3] = {0,0,0}, result_tot_Cgg[3] = {0,0,0};
	// double result_mtt_Cgg2[3] = {0,0,0}, result_costhetaCM_Cgg2[3] = {0,0,0}, result_y_Cgg2[3] = {0,0,0}, result_pT_Cgg2[3] = {0,0,0}, result_tot_Cgg2[3] = {0,0,0};
	// cout << " * Cggtilde(" << ScaleRenormalization << ") = " << Cggtilde_at_mu(ScaleRenormalization,Ctt,fa,mt,vev,MZ,pdf,0.) << endl;

	// ============ m_ttbar distribution ============
	double s_parton = HUGE_VALF; // this is a placeholder and can take any value. Assign NAN to check if the program runs correctly (the result shouldn't be NAN)
	struct ParamPDF pPDF = {pidg,pidg,ScaleFactorization,s_parton,s_hadron};
	struct IntegratePDF iPDF = {pdf,pPDF,w1,RelativeError,MaxNumberSubintervals};
	double integrateion_result_gg_1S01[3] = {0,0,0}, integrateion_result_gg_1S08[3] = {0,0,0}, integrateion_result_qq_3S18[3] = {0,0,0}, integrateion_result_gq_1S01[3] = {0,0,0};
	const string filename = OutputFolder + "/" + datestring + "/dist_mtt.dat", fileGreenFunction = OutputFolder + "/" + datestring + "/greenfunction.dat";
	cout << " * Exporting m_ttbar to " << filename << ". " << flush;
	ofstream f(filename.c_str());
	ofstream fGreenFunction(fileGreenFunction.c_str());
	f << "# Compute m_ttbar distribution dxsection/dm_ttbar (pb/GeV) of PP > ttbar. Energy and mass in GeV, xsection in pb." << endl;
	f << "# Parameters: alphaS(ScaleR) = " << alphaS_muR << ", m_t = " << mt << ", PDFset = " << PDFSetName << ", PDFmember = " << PDFmember << ", ScaleF = " << ScaleFactorization << ", ScaleR = " << ScaleRenormalization << ", sqrt(s_hadron) = " << sqrt(s_hadron) << ", MZ = " << MZ << ", vev = " << vev << ", RelativeError = " << RelativeError << ", DynamicScalingChoice = " << DynamicScalingChoice << ", ScaleUncertainty = (" << ScaleUncertainty[0] << "," << ScaleUncertainty[1] << ")" << endl;
	f << "# m_ttbar (GeV) \t\t\t\t gg > singlet \t [lower , \t upper] \t\t\t\t gg > octet \t [lower , \t upper] \t\t\t\t qq > octet \t [lower , \t upper]" << endl;
	fGreenFunction << "# m_ttbar (GeV) \t\t\t\t singlet \t\t\t\t octet" << endl;

	for (double mttbar = mttbarmin; mttbar <= mttbarmax+mttbarstep*0.5; mttbar += mttbarstep)
	{	
		s_parton = mypow(mttbar,2);
		double dsigmadmttbar_gg_1S01[3] = {0,0,0}, dsigmadmttbar_gg_1S08[3] = {0,0,0}, dsigmadmttbar_qq_3S18[3] = {0,0,0}, dsigmadmttbar_gq_1S01[3] = {0,0,0};
		cout << "[DEBUG] mttbar = " << mttbar; // DEBUG
		for (size_t i_scale = 0; i_scale < 3; i_scale++)
		{
			iPDF.param.s_parton = s_parton;
			struct ParamXsection pXsection = {pidg,pidg,s_parton,alphaS_muR,mt,ScaleRenormalization*ScaleUncertaintyFactor[i_scale],false,NAN,MZ,vev};
			pXsection.TopDecayWidth = TopDecayWidth;
			struct IntegrateXsection iXsection = {"",pXsection,w1,NAN,NAN,"",RelativeError,MaxNumberSubintervals};
			iPDF.param.ScaleFactorization = ScaleFactorization * ScaleUncertaintyFactor[i_scale];
			iXsection.workspace2 = w2;
			struct ParamConvolute1D pConvolute(iPDF,iXsection);
			void* argConvolute = &pConvolute;

			// gg > T(S=0,J=0,singlet)
			pConvolute.iXsection.name_f_xsection = "dxsect_dmToponium_partonic";
			pConvolute.setPID(pidg,pidg);
			pConvolute.iXsection.param.BoundStateColorConfig = 1;
			pConvolute.iXsection.param.BoundStateSpin = 0;
			pConvolute.iXsection.param.BoundStateJ = 0;
			argConvolute = &pConvolute;
			dsigmadmttbar_gg_1S01[i_scale] = dsigmadmttbar(mttbar,argConvolute);

			// gg > T(S=0,J=0,octet)
			pConvolute.setPID(pidg,pidg);
			pConvolute.iXsection.param.BoundStateColorConfig = 8;
			pConvolute.iXsection.param.BoundStateSpin = 0;
			pConvolute.iXsection.param.BoundStateJ = 0;
			argConvolute = &pConvolute;
			dsigmadmttbar_gg_1S08[i_scale] = 0; //dsigmadmttbar(mttbar,argConvolute);
			if (i_scale == 1)
			{
				cout << ", dL[" << pConvolute.iPDF.param.pid1 << "," << pConvolute.iPDF.param.pid2 << "]/dtau = " << Integrated_x1_PDF_Product(pConvolute.iPDF.workspace,pConvolute.iPDF.pdf,pConvolute.iPDF.param.pid1,pConvolute.iPDF.param.pid2,s_parton,pConvolute.iPDF.param.s_hadron,pConvolute.iPDF.param.ScaleFactorization,RelativeError,pConvolute.iPDF.IntegrationRelativeError,pConvolute.iPDF.NumberOfSubintervals); // DEBUG
			}

			// qq > T(S=0,J=0,octet)
			dsigmadmttbar_qq_3S18[i_scale] = 0;
			pConvolute.iXsection.param.BoundStateColorConfig = 8;
			pConvolute.iXsection.param.BoundStateSpin = 1;
			pConvolute.iXsection.param.BoundStateJ = 1;
			for (size_t i_quark = 0; i_quark < pidq_length; i_quark++)
			{
				pConvolute.setPID(pidq[i_quark],-pidq[i_quark]);
				argConvolute = &pConvolute;
				dsigmadmttbar_qq_3S18[i_scale] += 0;//dsigmadmttbar(mttbar,argConvolute);
				if (i_scale == 1)
				{
					cout << ", dL[" << pConvolute.iPDF.param.pid1 << "," << pConvolute.iPDF.param.pid2 << "]/dtau = " << Integrated_x1_PDF_Product(pConvolute.iPDF.workspace,pConvolute.iPDF.pdf,pConvolute.iPDF.param.pid1,pConvolute.iPDF.param.pid2,s_parton,pConvolute.iPDF.param.s_hadron,pConvolute.iPDF.param.ScaleFactorization,RelativeError,pConvolute.iPDF.IntegrationRelativeError,pConvolute.iPDF.NumberOfSubintervals); // DEBUG
				}
			}

			// gq > T(S=0,J=0,singlet)
			dsigmadmttbar_gq_1S01[i_scale] = 0;
			pConvolute.iXsection.param.BoundStateColorConfig = 1;
			pConvolute.iXsection.param.BoundStateSpin = 0;
			pConvolute.iXsection.param.BoundStateJ = 0;
			for (size_t i_quark = 0; i_quark < pidq_length; i_quark++)
			{
				pConvolute.setPID(pidg,pidq[i_quark]);
				argConvolute = &pConvolute;
				dsigmadmttbar_gq_1S01[i_scale] += dsigmadmttbar(mttbar,argConvolute);
				if (i_scale == 1)
				{
					cout << ", dL[" << pConvolute.iPDF.param.pid1 << "," << pConvolute.iPDF.param.pid2 << "]/dtau = " << Integrated_x1_PDF_Product(pConvolute.iPDF.workspace,pConvolute.iPDF.pdf,pConvolute.iPDF.param.pid1,pConvolute.iPDF.param.pid2,s_parton,pConvolute.iPDF.param.s_hadron,pConvolute.iPDF.param.ScaleFactorization,RelativeError,pConvolute.iPDF.IntegrationRelativeError,pConvolute.iPDF.NumberOfSubintervals); // DEBUG
				}
			}
			
			iPDF.param.ScaleFactorization = ScaleFactorization;
		}
		cout << endl; // DEBUG
		
		// Write to file
		f << mttbar << "\t\t\t\t" << dsigmadmttbar_gg_1S01[0] << "\t" << dsigmadmttbar_gg_1S01[1] << "\t" << dsigmadmttbar_gg_1S01[2] << "\t\t\t\t" << dsigmadmttbar_gg_1S08[0] << "\t" << dsigmadmttbar_gg_1S08[1] << "\t" << dsigmadmttbar_gg_1S08[2] << "\t\t\t\t" << dsigmadmttbar_qq_3S18[0] << "\t" << dsigmadmttbar_qq_3S18[1] << "\t" << dsigmadmttbar_qq_3S18[2] << "\t\t\t\t" << dsigmadmttbar_gq_1S01[0] << "\t" << dsigmadmttbar_gq_1S01[1] << "\t" << dsigmadmttbar_gq_1S01[2] << endl;
		double CF = 4./3., CA = 3., alphaS = pdf->alphasQ(ScaleRenormalization), SoftScaleRenormalization = mt*CF*alphaS, alphaS_Soft = pdf->alphasQ(SoftScaleRenormalization);
		fGreenFunction << mttbar << "\t\t\t\t" << ImGreenFunction(mttbar, mt, TopDecayWidth, alphaS_Soft, CF, CA, 1) << "\t\t\t\t" << ImGreenFunction(mttbar, mt, TopDecayWidth, alphaS_Soft, CF, CA, 8)  << endl;
		// Add up all bin (total should be consistent with total xsection)
		for (size_t i_scale = 0; i_scale < 3; i_scale++)
		{
			integrateion_result_gg_1S08[i_scale] += dsigmadmttbar_gg_1S08[i_scale]*mttbarstep;
			integrateion_result_qq_3S18[i_scale] += dsigmadmttbar_qq_3S18[i_scale]*mttbarstep;
			integrateion_result_gg_1S01[i_scale] += dsigmadmttbar_gg_1S01[i_scale]*mttbarstep;
		}
	}
	f_result << "# Integrated xsection from mtt distribution        (gg>1S01) \t" << integrateion_result_gg_1S01[0] << " [" << integrateion_result_gg_1S01[1] << "," << integrateion_result_gg_1S01[2] << "] pb; (gg>1S08) \t" << integrateion_result_gg_1S08[0] << " [" << integrateion_result_gg_1S08[1] << "," << integrateion_result_gg_1S08[2] << "] pb; (qq>3S18) \t" << integrateion_result_qq_3S18[0] << " [" << integrateion_result_qq_3S18[1] << "," << integrateion_result_qq_3S18[2] << "] pb." << endl;
	

	f.close();
	fGreenFunction.close();
	cout << "Done. Integrated total xsection = (gg>1S01) \t" << integrateion_result_gg_1S01[0] << " [" << integrateion_result_gg_1S01[1] << "," << integrateion_result_gg_1S01[2] << "] pb; (gg>1S08) \t" << integrateion_result_gg_1S08[0] << " [" << integrateion_result_gg_1S08[1] << "," << integrateion_result_gg_1S08[2] << "] pb; (qq>3S18) \t" << integrateion_result_qq_3S18[0] << " [" << integrateion_result_qq_3S18[1] << "," << integrateion_result_qq_3S18[2] << "] pb." << endl;

	// Computation ends
	f_result.close();
	t_end = chrono::high_resolution_clock::now();
	dif = std::chrono::duration<double, std::milli>(t_end - t_start).count() / 1000.;
	cout << "[INFO] Computation finshied (it took " << dif << " s)" << endl << endl;

	// Wrap up
	cout << "* Note: . alphaS_soft currently has no scale uncertainty band " << endl;
	cout << "        . Check that the other toponium papers are doing something similar." << endl;
	cout << "        . Go through: Peskin ch 17, 0812.0919, 1003.5827, hep-ph/9707223, 1207.2389." << endl;
	cout << "        . Can we write these in a way that can be used in MC event generator?" << endl;
	cout << "        . Has initial state radiation been already included in PDF?" << endl;
	cout << "        . What invariant mass to use?" << endl;
	cout << "        . What to do to solve: test on all scaleing available + extract pdf to check by hand." << endl;
	cout << "        . The SM prediction has no wavefunction effects, including Sommerfeld enhancement. These should be included in toponium, although they are not. See eq. 16 in 2411.18962." << endl;
	cout << ".................. Exiting .................." << endl;
	delete pdf;
	gsl_integration_workspace_free (w1);
	gsl_integration_workspace_free (w2);
	gsl_integration_workspace_free (w3);
	t_programend = chrono::high_resolution_clock::now();
	programdif = std::chrono::duration<double, std::milli>(t_programend - t_programstart).count() / 1000.;
	cout << "[INFO] Program finshied (it took " << int(programdif/60) << " mins " << fmod(programdif,60) << " s)" << endl << endl;
	return 0;
}
/********************************************************************************
 *									  *										    *
 * 	Contains physical formulas (list can be found in Formulas.h)                *
 * 	Use: LHAPDF 6.5.3, NNPDF 3.1 (NNPDF31_nnlo_as_0118_notop), GSL				*
 *									  *										    *
 ********************************************************************************/


#include "Formulas.h"


/********************************************************************************	
*									  *										    *
*	    Compute Jacobian and costhetaCM from pT and y for gg > ttbar			*
*									  *									        *
********************************************************************************/
// Warning: assume diffxsection is symmetric t <-> u
void Jacobian_and_Conversion_ggttbar_pT(double pT, double s_parton, double mt, double &costhetaCM, double &Jacobian)
{
	if (s_parton <= 4*mypow(mt,2))
	{
		throw invalid_argument("Energy not conserved (in Jacobian_and_Conversion_pT): Expect s_parton > 4*mt^2");
	}
	double mt2 = mypow(mt,2);
	double pCM = sqrt(s_parton/4. - mt2);
	
	if (pT < 0)
	{
		throw invalid_argument("Incorrect value for pT (in Jacobian_and_Conversion_pT): Expect pT >= 0");
	}
	if (s_parton < 0)
	{
		throw invalid_argument("Incorrect value for s (in Jacobian_and_Conversion_pT): Expect s >= 0");
	}
	if (mt < 0)
	{
		throw invalid_argument("Incorrect value for mt (in Jacobian_and_Conversion_pT): Expect mt >= 0");
	}
	if (pT > pCM)
	{
		Jacobian = 0;
		costhetaCM = NAN;
		cout << "Incorrect value for pT (in Jacobian_and_Conversion_pT): Expect pT <= pCM" << endl;
		throw invalid_argument("Incorrect value for pT (in Jacobian_and_Conversion_pT): Expect pT <= pCM");
		return;
	}
	
	double sinthetaCM = pT/pCM;
	costhetaCM = sqrt(1. - mypow(sinthetaCM,2));
	double tanthetaCM = sinthetaCM / costhetaCM;
	Jacobian = 2.*tanthetaCM/pCM;	// factor 2 from assuming diffxsection is symmetric
}

void Jacobian_and_Conversion_ggttbar_rapidity(double rapidity, double s_parton, double mt, double &costhetaCM, double &Jacobian)
{
	if (s_parton <= 4*mypow(mt,2))
	{
		throw invalid_argument("Energy not conserved (in Jacobian_and_Conversion_pT): Expect s_parton > 4*mt^2");
	}
	if (s_parton < 0)
	{
		throw invalid_argument("Incorrect value for s (in Jacobian_and_Conversion_pT): Expect s >= 0");
	}
	if (mt < 0)
	{
		throw invalid_argument("Incorrect value for mt (in Jacobian_and_Conversion_pT): Expect mt >= 0");
	}

	double mt2 = mypow(mt,2);
	double ECM = sqrt(s_parton)/2;
	double pCM = sqrt(s_parton/4. - mt2);
	costhetaCM = ECM/pCM * (exp(2*rapidity) - 1) / (exp(2*rapidity) + 1);
	if (costhetaCM < -1 || costhetaCM > 1)
	{
		throw invalid_argument("Incorrect value for rapidity (in Jacobian_and_Conversion_pT): rapidity out of physical range");
	}

	Jacobian = ECM/pCM * (4*exp(2*rapidity)) / mypow(1 + exp(2*rapidity),2);
}

// Kinematical cut
bool kinematical_cut(double costhetaCM, double s_parton, double mt)
{
	return true;
	// double sinthetaCM = sqrt(1. - mypow(costhetaCM,2));
	// double pCM = sqrt(s_parton/4. - mypow(mt,2));
	// double pT = pCM*sinthetaCM;

	// if (pT >= 400. && pT <= 1500.) // GeV
	// {
	// 	return true;
	// }
	// return false;
}

/********************************************************************************	
*									  *										    *
*						    Compute dynamic scale							    *
*									  *									        *
********************************************************************************/
double DynamicScale(double s, double mt, double costhetaCM, int DynamicScaleChoice)
{
	if (DynamicScaleChoice == 0)
	{
		throw domain_error("Invalid DynamicScaleChoice (in DynamicScale). DynamicScaleChoice = 0 is for static scale");
	}
	else if (DynamicScaleChoice == 2)
	{
		double pT;
		double mt2 = mypow(mt,2);
		double pCM = sqrt(s/4. - mt2), sinthetaCM = sqrt(1. - mypow(costhetaCM,2));
		pT = pCM*sinthetaCM;
		// Dynamic scale = sum_f mT(f)/2 with mT = sqrt(m^2 + pT^2) and f is final states. For 2->2, scale = mT(t)
		double scale = sqrt(mt2 + mypow(pT,2));
		return scale;
	}
	else
	{
		throw domain_error("Incorrect option for DynamicScalingChoice: must be either 0 (fixed scale) or 2 (half sum transverse mass of final particles).");
	}
	return NAN;
}


/********************************************************************************	
*									  *										    *
*  			    Compute t,u from s, masses, costhetaCM (gg > ttbar)				*
*									  *									        *
********************************************************************************/
void costhetaCM_to_t_u_ggtt(double costhetaCM, double &t, double &u, double s, double mt)
{
	if (s <= 4*mypow(mt,2))
	{
		throw invalid_argument("Energy not conserved (in costhetaCM_to_t_u_ggtt): Expect s > 4*mt^2");
	}
	if (costhetaCM < -1 or costhetaCM > 1)
	{
		throw invalid_argument("Incorrect value for costhetaCM (in costhetaCM_to_t_u_ggtt): Expect -1 <= costhetaCM <= 1");
	}
	if (s < 0)
	{
		throw invalid_argument("Incorrect value for s (in costhetaCM_to_t_u_ggtt): Expect s >= 0");
	}
	if (mt < 0)
	{
		throw invalid_argument("Incorrect value for mt (in costhetaCM_to_t_u_ggtt): Expect mt >= 0");
	}
	
	double mt2 = mypow(mt,2);
	t = mt2 + (-s + costhetaCM*sqrt(s*(-4*mt2 + s)))/2.;
	u = 2*mt2 - s - t;
}


/********************************************************************************	
*									  *										    *
*			Tree level dxsect/dy for g + g -> t + tbar (in GeV^-2)				*
*									  *									        *
********************************************************************************/

double dxsection_dy_ggttTree(double alphaS, double mt, double s, double y)
{
    if (s <= 4*mypow(mt,2))
    {
        return 0;
    }
	if (alphaS < 0 || alphaS > sqrt(4*M_PI))
	{
		throw invalid_argument("Incorrect value for alphaS (in dxsection_dcosthetaCM_ggttTree): Expect 0 <= alphaS <= 4*Pi");
	}
    double costhetaCM, Jacobian;
	Jacobian_and_Conversion_ggttbar_rapidity(y,s,mt,costhetaCM,Jacobian);

    return Jacobian*dxsection_dcosthetaCM_ggttTree(alphaS,mt,s,costhetaCM);
}


/********************************************************************************	
*									  *										    *
*		Tree level dxsection/dcosthetaCM for g + g -> t + tbar (in GeV^-2)		*
*									  *									        *
********************************************************************************/
// pCM = sqrt(s-4*mt2)
double dxsection_dcosthetaCM_ggttTree(double alphaS, double mt, double s, double costhetaCM)
{
    if (s <= 4*mypow(mt,2))
    {
        return 0;
    }
	if (alphaS < 0 || alphaS > sqrt(4*M_PI))
	{
		throw invalid_argument("Incorrect value for alphaS (in dxsection_dcosthetaCM_ggttTree): Expect 0 <= alphaS <= 4*Pi");
	}
    double t, u, mt2 = mypow(mt,2);
	costhetaCM_to_t_u_ggtt(costhetaCM,t,u,s,mt);
	double pCM = sqrt(s/4. - mt2);
	
    return -(1./96.)*2*M_PI*(mypow(alphaS,2)*2.*pCM*(2*mypow(mt2,4) - 8*mypow(mt2,3)*t + mypow(mt2,2)*(3*mypow(s,2) + 4*s*t + 12*mypow(t,2)) + 
        t*(mypow(s,3) + 3*mypow(s,2)*t + 4*s*mypow(t,2) + 2*mypow(t,3)) - mt2*(mypow(s,3) + 2*mypow(s,2)*t + 8*s*mypow(t,2) + 8*mypow(t,3)))*
      (7*mypow(mt2,2) + 4*mypow(t,2) - t*u + 4*mypow(u,2) - 7*mt2*(t + u)))/(mypow(s,3.5)*mypow(mt2 - t,2)*mypow(mt2 - u,2));
}

double dxsection_dcosthetaCM_ggttTree_components(double alphaS, double mt, double s, double costhetaCM, double ts, double tt, double tu)
{
    if (s <= 4*mypow(mt,2))
    {
        return 0;
    }
	if (alphaS < 0 || alphaS > sqrt(4*M_PI))
	{
		throw invalid_argument("Incorrect value for alphaS (in dxsection_dcosthetaCM_ggttTree): Expect 0 <= alphaS <= 4*Pi");
	}
    double t, u, mt2 = mypow(mt,2);
	costhetaCM_to_t_u_ggtt(costhetaCM,t,u,s,mt);
	double pCM = sqrt(s/4. - mt2);
	
    return -(1./96.)*2*M_PI*(mypow(alphaS,2)*2.*pCM*(2*mypow(mt2,4) - 8*mypow(mt2,3)*t + mypow(mt2,2)*(3*mypow(s,2) + 4*s*t + 12*mypow(t,2)) + t*(mypow(s,3) + 3*mypow(s,2)*t + 4*s*mypow(t,2) + 2*mypow(t,3)) - mt2*(mypow(s,3) + 2*mypow(s,2)*t + 8*s*mypow(t,2) + 8*mypow(t,3)))*(4*mypow(t,2)*mypow(tu,2) + mypow(mt2,2)*(4*mypow(tt,2) - tt*tu + 4*mypow(tu,2)) - t*tt*tu*u + 4*mypow(tt,2)*mypow(u,2) + mt2*(t*(tt - 8*tu)*tu + tt*(-8*tt + tu)*u)))/(mypow(s,3.5)*mypow(mt2 - t,2)*mypow(mt2 - u,2));
}


/********************************************************************************	
*									  *										    *
*								For toponium									*
*									  *									        *
********************************************************************************/
// See table 1 and eq. 6 in 0812.0919
double normalization_factor(int pid1, int pid2, double boundstateSpin, double boundstateJ, int boundstateColorConfig)
{
	if (boundstateSpin == 0 and boundstateJ == 0)
	{
		if (pid1 != 21 && pid2 != 21)	// Both initial partons are quarks
		{
			if (boundstateColorConfig == 1)	// Color singlet
			{
				return 3./4.;
			}
			else if (boundstateColorConfig == 8)	// Color octet
			{
				return 6;
			}
			else
			{
				throw invalid_argument("Incorrect value for boundstateColorConfig (in NormalizationFactor): Expect 1 or 8");
			}
		}
		else // At least one of the initial parton is a gluon
		{
			if (boundstateColorConfig == 1)	// Color singlet
			{
				return 1;
			}
			else if (boundstateColorConfig == 8)	// Color octet
			{
				return 5./2.;
			}
			else
			{
				throw invalid_argument("Incorrect value for boundstateColorConfig (in NormalizationFactor): Expect 1 or 8");
			}
		}
	}
	else if (boundstateSpin == 1 and boundstateJ == 1)
	{
		if (pid1 == 21 && pid2 == 21)	// Both initial partons are gluons
		{
			if (boundstateColorConfig == 1)	// Color singlet
			{
				return 9./4.;
			}
			else if (boundstateColorConfig == 8)	// Color octet
			{
				return 18;
			}
			else
			{
				throw invalid_argument("Incorrect value for boundstateColorConfig (in NormalizationFactor): Expect 1 or 8");
			}
		}
		else // At least one of the initial parton is a quark
		{
			if (boundstateColorConfig == 1)	// Color singlet
			{
				return 0;
			}
			else if (boundstateColorConfig == 8)	// Color octet
			{
				return 32./3.;
			}
			else
			{
				throw invalid_argument("Incorrect value for boundstateColorConfig (in NormalizationFactor): Expect 1 or 8");
			}
		}
	}
	else
	{
		throw invalid_argument("Incorrect value for boundstateSpin and boundstateJ (in F_hardfactor): Expect (0,0) or (1,1)");
	}
	return 0;
}

// First term, second line of eq. 6 in 0812.0919
double is_two_to_one_process(int pid1, int pid2, double boundstateSpin, double boundstateJ, int boundstateColorConfig, double z)
{
	if (z != 1)	// Check that this condition can actually be satisfied
	{
		return 0;
	}
	if (pid1 == 21 && pid2 == 21 && boundstateSpin == 0 && boundstateJ == 0)
	{
		return 1;
	}
	else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0 && boundstateSpin == 1 && boundstateJ == 1)
	{
		return 1;
	}
	return 0;
}

// See eq. 7 in 0812.0919
double Ch_coefficient(int pid1, int pid2, double boundstateSpin, double boundstateJ, int boundstateColorConfig, double mtt, double ScaleRenormalization, double beta0, double CF, double CA, int nf, double TF)
{
	if (pid1 == 21 && pid2 == 21 && boundstateSpin == 0 && boundstateJ == 0 && boundstateColorConfig == 1)
	{
		return beta0/2 * log(mypow(ScaleRenormalization/mtt,2)) + CF*(mypow(M_PI,2)/4. - 5) + CA*(1 + mypow(M_PI,2)/12.);
	}
	else if (pid1 == 21 && pid2 == 21 && boundstateSpin == 0 && boundstateJ == 0 && boundstateColorConfig == 8)
	{
		return beta0/2 * log(mypow(ScaleRenormalization/mtt,2)) + CF*(mypow(M_PI,2)/4. - 5) + CA*(3 - mypow(M_PI,2)/24.);
	}
	else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0 && boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8)	// Must the initial quarks be the same?
	{
		return beta0/2 * log(mypow(ScaleRenormalization/mtt,2)) + CF*(mypow(M_PI,2)/3. - 8) + CA*(59./9. + 2*log(2.)/3. - mypow(M_PI,2)/4.) - 10./9.*nf*TF - 16./9.*TF;
	}
	return 0;
}

// See eq. 10 in 0812.0919
double splitting_function(int pid1, int pid2, double z, double CA, double CF, double TF)
{
	if (pid1 == 21 and pid2 == 21)
	{
		return 2*CA*(1./(1-z) + 1./z  + z*(1-z) - 2.);
	}
	else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0)
	{
		return 0;
	}
	else
	{
		throw invalid_argument("Incorrect value for pid1 and pid2 (in splitting_function): Expect pid1 + pid2 = 0 or at least one of them is a gluon (21)");
	}
	return 0;
}

/**
 * @brief Hard factor for gg > T (see eq. 5 and 6 in 0812.0919), where T = (2S+1) S_J^[1,8], where S is the bound state with spin and total angular momentum J and color singlet ([1]) or octet ([8]).
 */
double F_hardfactor(int pid1, int pid2, double s_parton, double mtt, double alphaS, double ScaleFactorization, double ScaleRenormalization, double boundstateSpin, double boundstateJ, int boundstateColorConfig, double CF, double CA, double TF)
{
	double z = mypow(mtt,2)/s_parton;
	// Input check
	if (z > 1)
	{
		throw invalid_argument("Incorrect value for z (in F_hardfactor): Expect z <= 1");
	}
	if (z < 0)
	{
		throw invalid_argument("Incorrect value for z (in F_hardfactor): Expect z >= 0");
	}
	if (boundstateSpin != 0 && boundstateSpin != 1)
	{
		throw invalid_argument("Incorrect value for boundstateSpin (in F_hardfactor): Expect 0 or 1");
	}
	if (boundstateJ != 0 && boundstateJ != 1)
	{
		throw invalid_argument("Incorrect value for boundstateJ (in F_hardfactor): Expect 0 or 1");
	}
	if (boundstateColorConfig != 1 && boundstateColorConfig != 8)
	{
		throw invalid_argument("Incorrect value for boundstateColorConfig (in F_hardfactor): Expect 1 or 8");
	}
	// The initial quarks must be the each other antiparticle
	if (pid1 != 21 && pid2 != 21 && pid1 + pid2 != 0)
	{
		throw invalid_argument("Incorrect value for pid1 and pid2 (in F_hardfactor): Expect pid1 + pid2 = 0 or at least one of them is a gluon (21)");
	}

	int nf = 5;	// Number of active flavors
	double beta0 = 11./3.*CA - 4./3.*nf*TF;	// beta function at one loop
	double N = normalization_factor(pid1, pid2, boundstateSpin, boundstateJ, boundstateColorConfig);
	double delta = is_two_to_one_process(pid1, pid2, boundstateSpin, boundstateJ, boundstateColorConfig, z);
	double Ch = Ch_coefficient(pid1, pid2, boundstateSpin, boundstateJ, boundstateColorConfig, mtt, ScaleRenormalization, beta0, CF, CA, nf, TF);
	double Acollinear = 0;
	double Anoncollinear = 0;

	return N*mypow(M_PI,2)*mypow(alphaS,2)/(3*s_parton) * (1+alphaS*Ch/M_PI) * (mypow(mtt,2)*delta + alphaS/M_PI*(Acollinear + Anoncollinear));
}

double ImGreenFunction(double mtt, double mt, double TopDecayWidth, double alphaS, double CF, double CA, int boundstateColorConfig)
{
	double C = 0;
	if (boundstateColorConfig == 1)
	{
		C = CF;
	}
	else if (boundstateColorConfig == 8)
	{
		C = CF - CA/2.;
	}
	else
	{
		throw invalid_argument("Incorrect value for boundstateColorConfig (in ImGreenFunction): Expect 1 or 8");
	}

	complex<double> G = 0, I(0,1), v = sqrt((mtt + I*TopDecayWidth - 2*mt)/mt);
	// double argv = atan2(TopDecayWidth, mtt - 2*mt)/2, abs_v = sqrt(mypow((mtt-2*mt)/mt,2) + mypow(TopDecayWidth/mt,2));
	// v = abs_v * exp(I*argv);
	G = mypow(mt,2) * v/(4*M_PI) * (I + alphaS*C/v * (I*M_PI/2. - log(v)));
	// return mypow(mt,2) * ( 1./(4*M_PI)*abs_v*cos(argv) + alphaS*C/8. - alphaS*C/(4*M_PI)*argv);
	return imag(G);
}

double dxsect_dmToponium_partonic(int pid1, int pid2, double s_parton, double mtt, double ScaleFactorization, double ScaleRenormalization, double mt, double TopDecayWidth, double boundstateSpin, double boundstateJ, int boundstateColorConfig, LHAPDF::PDF* pdf)
{
	double alphaS = pdf->alphasQ(ScaleRenormalization);
	double CF = 4./3.;	// Casimir operator for fundamental representation
	double CA = 3.;	// Casimir operator for adjoint representation
	double TF = 1./2.;
	double F = F_hardfactor(pid1, pid2, s_parton, mtt, alphaS, ScaleFactorization, ScaleRenormalization, boundstateSpin, boundstateJ, boundstateColorConfig, CF, CA, TF);
	double SoftScaleRenormalization = mt*CF*alphaS;
	double alphaS_Soft = pdf->alphasQ(SoftScaleRenormalization);
	double ImG = ImGreenFunction(mtt, mt, TopDecayWidth, alphaS_Soft, CF, CA, boundstateColorConfig);
    return 1./mtt * F * ImG / mypow(mt,2);
}
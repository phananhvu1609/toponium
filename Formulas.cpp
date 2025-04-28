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
	else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0 && boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8)
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

// See eq. 10 in 0812.0919. For P_qq and P_gg, (1-z)*P_qq and (1-z)*P_gg is implemented instead, to avoid singularities. IMPORTANT: P_gq != P_qg
double splitting_function(int pid1, int pid2, double z, double CA, double CF, double TF)
{
	if (pid1 == 21 && pid2 == 21)
	{
		return 2*CA*(1 + (1-z)/z  + z*(1-z)*(1-z) - 2.*(1-z));
	}
	else if (pid1 == 21 && pid2 != 21)
	{
		return CF*(1 + mypow(1-z,2))/z;
	}
	else if (pid1 != 21 && pid2 == 21)
	{
		return TF*(mypow(z,2) + mypow(1-z,2));
	}
	else if (pid1 != 21 && pid2 != 21)
	{
		return 2*CF*(1 - (1+z)*(1-z)/2);
	}
	else
	{
		throw invalid_argument("Incorrect values for pid1 and pid2 (in splitting_function): Expect either quarks (1,2,3,4,5,6) or gluon (21)");
	}
	return 0;
}

double Acollinear_function(int pid1, int pid2, double z, double mtt, double ScaleRenormalization, double ScaleFactorization, double boundstateSpin, double boundstateJ, int boundstateColorConfig, bool is_two_to_one, double beta0, double CF, double CA, double TF, bool is_subtract = false)
{
	if (is_two_to_one)
	{
		// g+g -> 1 S_0^[1,8] + X
		if (pid1 == 21 && pid2 == 21 && boundstateSpin == 0 && boundstateJ == 0 && z == 1)
		{
			return -beta0/2.*log(mypow(ScaleFactorization/mtt,2));
		}
		else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0 && boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8 && z == 1)
		{
			return -3./2.*CF*log(mypow(ScaleFactorization/mtt,2));
		}
		else
		{
			return 0;
		}
		if (!is_two_to_one_process(pid1, pid2, boundstateSpin, boundstateJ, boundstateColorConfig, z))
		{
			return 0;
		}
	}
	// g+q -> T+X
	else if (((pid1 == 21 && pid2 != 21) || (pid1 != 21 && pid2 == 21)) && !is_subtract)
	{
		// g+q -> 1 S_0^[1,8] + X
		if (boundstateSpin == 0 && boundstateJ == 0)
		{
			return 1./2.*splitting_function(21,1,z,CA,CF,TF)*log( mypow(mtt*(1-z)/ScaleFactorization,2) / z ) + CF*z/2;
		}
		// g+q -> 3 S_1^[8] + X
		else if (boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8)
		{
			return 1./2.*splitting_function(1,21,z,CA,CF,TF)*log( mypow(mtt*(1-z)/ScaleFactorization,2) / z ) + TF*z*(1-z);
		}
		else
		{
			return 0;
		}
	}
	// g+g -> T + X
	else if (pid1 == 21 && pid2 == 21)
	{
		// g+q -> 1 S_0^[1,8] + X
		if (boundstateSpin == 0 && boundstateJ == 0)
		{
			if (!is_subtract)	// if not subtract, return normal expression
			{
				return splitting_function(pid1,pid2,z,CA,CF,TF)*(2*log(1-z)/(1-z) + 1/(1-z)*log( mypow(mtt/ScaleFactorization,2)/z )); // splitting function already include (1-z)
			}
			else	// if subtract, all terms regular at z = 1 are set to z = 1
			{	
				return splitting_function(pid1,pid2,1,CA,CF,TF)*(2*log(1-z)/(1-z) + 1/(1-z)*log( mypow(mtt/ScaleFactorization,2) )); // splitting function already include (1-z)
			}
		}
		else
		{
			return 0;
		}
	}
	// q+q -> T + X
	else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0)
	{
		// q+q -> 3 S_1^[8] + X
		if (boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8)
		{
			if (!is_subtract) // if not subtract, return normal expression
			{
				return splitting_function(pid1,pid2,z,CA,CF,TF)*(2*log(1-z)/(1-z) + 1/(1-z)*log( mypow(mtt/ScaleFactorization,2)/z )) + CF*(1-z); // splitting function already include (1-z)
			}
			else	// if subtract, all terms regular at z = 1 are set to z = 1
			{
				return splitting_function(pid1,pid2,1,CA,CF,TF)*(2*log(1-z)/(1-z) + 1/(1-z)*log( mypow(mtt/ScaleFactorization,2) )); // splitting function already include (1-z)
			}
		}
		else
		{
			return 0;
		}
	}
	return 0;
}

double Anoncollinear_function(int pid1, int pid2, double z, double mtt, double ScaleRenormalization, double ScaleFactorization, double boundstateSpin, double boundstateJ, int boundstateColorConfig, double beta0, double CF, double CA, double TF, bool is_subtract = false)
{
	double z2 = mypow(z,2), z3 = mypow(z,3), z4 = mypow(z,4), z5 = mypow(z,5), z6 = mypow(z,6), z7 = mypow(z,7), z8 = mypow(z,8);
	double Nc = 3;	// Number of colors
	double BF = (mypow(Nc,2)-4)/(4*Nc);
	// g+g -> T + X
	if (pid1 == 21 && pid2 == 21)
	{
		// g+g -> 1 S_0^[1] + X
		if (boundstateSpin == 0 && boundstateJ == 0 && boundstateColorConfig == 1 && !is_subtract)
		{
			if (z <= 0.97)
			{
				return -CA / (6*z*mypow(1-z,2)*mypow(1+z,3)) * (12 + 11*z2 + 24*z3 - 21*z4 - 24*z5 + 9*z6 - 11*z8 + 12*(-1 + 5*z2 + 2*z3 + z4 + 3*z6 + 2*z7)*log(z) );
			}
			else	// The above expression ill-behave at z = 1. The one below is an exppansion of the above expression at z = 1, to regularize the expression
			{
				return -CA / 6. * (11 + 21*mypow(-1 + z,2) - (152*mypow(-1 + z,3))/5. + (168*mypow(-1 + z,4))/5. - (1289*mypow(-1 + z,5))/35. + (1383*mypow(-1 + z,6))/35. - (1459*mypow(-1 + z,7))/35. + (652*mypow(-1 + z,8))/15. - (103877*mypow(-1 + z,9))/2310. + (3239*mypow(-1 + z,10))/70. + z);
			}
		}
		// g+g -> 1 S_0^[8] + X
		else if (boundstateSpin == 0 && boundstateJ == 0 && boundstateColorConfig == 8)
		{
			double fz = -CA; // if subtract, fz = f(z=1)
			if (z <= 0.97 && !is_subtract)
			{
				fz = CA*(-2 - (23*z2)/6. - 5*z3 + (7*z4)/2. + 4*z5 - (3*z6)/2. + z7 + (23*z8)/6. - 2*(-1 + 5*z2 + 2*z3 + 3*z4 + 5*z6 + 2*z7)*log(z))/((1 - z)*z*mypow(1 + z,3));
			}
			else if (!is_subtract)	// The above expression ill-behave at z = 1. The one below is an exppansion of the above expression at z = 1, to regularize the expression
			{
				fz = CA*(-1 - (5*mypow(-1 + z,2))/2. + (11*mypow(-1 + z,3))/6. - (26*mypow(-1 + z,4))/5. + (28*mypow(-1 + z,5))/5. - (128*mypow(-1 + z,6))/21. + (229*mypow(-1 + z,7))/35. - (2179*mypow(-1 + z,8))/315. + (4553*mypow(-1 + z,9))/630. - (1647*mypow(-1 + z,10))/220.);
			}
			return fz / (1-z);
		}
		// g+g -> 3 S_1^[1] + X
		else if (boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 1 && !is_subtract)
		{
			if (z <= 0.97)
			{
				return 256*BF/(6.*CF*mypow(Nc,2)) * z/(mypow(1-z,2)*mypow(1+z,3)) * (2 + z + 2*z2 - 4*z4 - z5 + 2*z2*(5 + 2*z + z2)*log(z));
			}
			else	// The above expression ill-behave at z = 1. The one below is an exppansion of the above expression at z = 1, to regularize the expression
			{
				return BF/(CF*mypow(Nc,2)) * ((-320*(-1 + z))/9. - (64*mypow(-1 + z,2))/9. + (448*mypow(-1 + z,3))/45. - (128*mypow(-1 + z,4))/15. + 
												(608*mypow(-1 + z,5))/105. - (352*mypow(-1 + z,6))/105. + (224*mypow(-1 + z,7))/135. - 
												(608*mypow(-1 + z,8))/945. + (1088*mypow(-1 + z,9))/10395. + (1472*mypow(-1 + z,10))/10395.);
			}
		}
		// g+g -> 3 S_1^[8] + X
		else if (boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8 && !is_subtract)
		{
			if (z <= 0.97)
			{
				return 1/(36*z*mypow(1-z,2)*mypow(1+z,3)) * (108 + 153*z + 400*z2 + 65*z3 - 356*z4 - 189*z5 - 152*z6 - 29*z7 + (108*z + 756*z2 + 432*z3 + 704*z4 + 260*z5 + 76*z6)*log(z));
			}
			else	// The above expression ill-behave at z = 1. The one below is an exppansion of the above expression at z = 1, to regularize the expression
			{
				return (-95*(-1 + z))/108. + (29*mypow(-1 + z,2))/27. - (487*mypow(-1 + z,3))/270. + (421*mypow(-1 + z,4))/180. - 
						(337*mypow(-1 + z,5))/126. + (3599*mypow(-1 + z,6))/1260. - (6665*mypow(-1 + z,7))/2268. + 
						(67247*mypow(-1 + z,8))/22680. - (13204*mypow(-1 + z,9))/4455. + (92051*mypow(-1 + z,10))/31185.;
			}
		}
		else
		{
			return 0;
		}
	}
	// g+q -> T+X
	else if (((pid1 == 21 && pid2 != 21) || (pid1 != 21 && pid2 == 21)) && !is_subtract)
	{
		// g+q -> 1 S_0^[1,8] + X
		if (boundstateSpin == 0 && boundstateJ == 0)
		{
			return -CF / z * (1-z) * (1-log(z));
		}
		// g+q -> 3 S_1^[8] + X
		else if (boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8)
		{
			return TF/4 * (1-z) * (1+3*z) + CA/(4*CF) * 1./z * ( (1-z)*(2+z+2*mypow(z,2)) + 2*z*(1+z)*log(z) ); 
		}
		else
		{
			return 0;
		}
	}
	// q+q -> T + X
	else if (pid1 != 21 && pid2 != 21 && pid1 + pid2 == 0)
	{
		// q+q -> 1 S_0^[1] + X
		if (boundstateSpin == 0 && boundstateJ == 0 && boundstateColorConfig == 1 && !is_subtract)
		{
			return 32*CF/(3*mypow(Nc,2)) * z*(1-z);
		}
		// q+q -> 1 S_0^[8] + X
		else if (boundstateSpin == 0 && boundstateJ == 0 && boundstateColorConfig == 8 && !is_subtract)
		{
			return 32*BF/(3*mypow(Nc,2)) * z*(1-z);
		}
		// q+q -> 3 S_1^[8] + X
		else if (boundstateSpin == 1 && boundstateJ == 1 && boundstateColorConfig == 8)
		{
			if (!is_subtract)	// if not subtract, return normal expression
			{
				if (z == 1)
				{
					throw invalid_argument("Incorrect value for z (in Anoncollinear_function): Expect z != 1");
				}
				return -( CF*mypow(1-z,2) + CA*(1+z+z2)/3. ) / (1-z);
			}
			else	// if subtract, all terms regular at z = 1 are set to z = 1
			{
				return - CA/(1-z);
			}
		}
		else
		{
			return 0;
		}
	}
	return 0;
}

/**
 * @brief Hard factor for gg > T (+X) (see eq. 5 and 6 in 0812.0919), where T = (2S+1) S_J^[1,8], where S is the bound state with spin and total angular momentum J and color singlet ([1]) or octet ([8]).
 * is_two_to_one is true for terms with delta(1-z) and false for terms without. subtract is for terms added to subtract 1/(1-z) divergence
 */
double F_hardfactor(int pid1, int pid2, double s_parton, double mtt, double alphaS, double ScaleRenormalization, double ScaleFactorization, double boundstateSpin, double boundstateJ, int boundstateColorConfig, double CF, double CA, double TF, bool is_two_to_one, bool is_subtract)
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
	if (pid1 != 21 && pid2 != 21 && abs(pid1) != 1 && abs(pid2) != 1 && abs(pid1) != 2 && abs(pid2) != 2 && abs(pid1) != 3 && abs(pid2) != 3 && abs(pid1) != 4 && abs(pid2) != 4 && abs(pid1) != 5 && abs(pid2) != 5 && abs(pid1) != 6 && abs(pid2) != 6)
	{
		throw invalid_argument("Incorrect value for pid1 and pid2 (in F_hardfactor): Expect either quarks (+-1, +-2, +-3, +-4, +-5, +-6) or gluons (21)");
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
	if (!is_two_to_one || is_subtract)
	{
		delta = 0;
	}
	double Ch = Ch_coefficient(pid1, pid2, boundstateSpin, boundstateJ, boundstateColorConfig, mtt, ScaleRenormalization, beta0, CF, CA, nf, TF);
	double Acollinear = Acollinear_function(pid1,pid2,z,mtt,ScaleRenormalization,ScaleFactorization,boundstateSpin,boundstateJ,boundstateColorConfig,is_two_to_one,beta0,CF,CA,TF,is_subtract);
	double Anoncollinear = 0;
	if (!is_two_to_one)
	{
		Anoncollinear = Anoncollinear_function(pid1,pid2,z,mtt,ScaleRenormalization,ScaleFactorization,boundstateSpin,boundstateJ,boundstateColorConfig,beta0,CF,CA,TF,is_subtract);
	}
	double replace_s_parton_by_mtt = 1;
	if (is_subtract)
	{
		replace_s_parton_by_mtt = 1/z; // when subtract, the 1/(3*s_parton) term should be replaced by 1/(3*mtt^2)
	}

	return N*mypow(M_PI,2)*mypow(alphaS,2)/(3*s_parton) * (1+alphaS*Ch/M_PI) * (delta + alphaS/M_PI*(Acollinear + Anoncollinear)*replace_s_parton_by_mtt);
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

// double dxsect_dmToponium_partonic(int pid1, int pid2, double s_parton, double mtt, double ScaleRenormalization, double ScaleFactorization, double mt, double TopDecayWidth, double boundstateSpin, double boundstateJ, int boundstateColorConfig, LHAPDF::PDF* pdf)
// {
// 	double alphaS = pdf->alphasQ(ScaleRenormalization);
// 	double CF = 4./3.;	// Casimir operator for fundamental representation
// 	double CA = 3.;	// Casimir operator for adjoint representation
// 	double TF = 1./2.;
// 	double F = F_hardfactor(pid1, pid2, s_parton, mtt, alphaS, ScaleRenormalization, ScaleFactorization, boundstateSpin, boundstateJ, boundstateColorConfig, CF, CA, TF);
// 	double SoftScaleRenormalization = mt*CF*alphaS;
// 	double alphaS_Soft = pdf->alphasQ(SoftScaleRenormalization);
// 	double ImG = ImGreenFunction(mtt, mt, TopDecayWidth, alphaS_Soft, CF, CA, boundstateColorConfig);
//     return 1./mtt * F * ImG / mypow(mt,2);
// }
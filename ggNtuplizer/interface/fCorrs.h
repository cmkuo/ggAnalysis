#ifndef INC_FCORRS
#define INC_FCORRS

namespace fcorrs{

//--------------- B A R R E L -----------------

  //
  //Method to correct energy loss due to the leackage of a shower from the crystals
  //
  inline double f5x5( double iEta ) {
    if ( iEta < 40.2198 ) return 1;
    return 1 - 3.03103e-6*(iEta - 40.2198)*(iEta - 40.2198);
  }

  inline double f3x3( double iEta ) {
    return 1 + 1.15581e-6*iEta*iEta - 1.35535e-10*iEta*iEta*iEta*iEta;
  }

  inline double fCorrEta( double iEta) {
    // old corrections DO NOT USE!

    Double_t fcorreta = 0;
    double p0 = 56.82;
    if ( fabs(iEta) < p0 ) {
      fcorreta = 0.99879;
    }
    if ( fabs(iEta) > p0 && fabs(iEta) <85.0 ) {
      double p1 = 1.006;
      double p2 = -2.227E-6;
      double p3 = 7.592E-11;
      fcorreta = p1 + p2 * fabs(iEta)*fabs(iEta) + p3 * fabs(iEta)*fabs(iEta)*fabs(iEta)*fabs(iEta);
    }
    if ( fcorreta == 0 )
      {
	std::cout << "Something is not right with Correction Function of Eta!!!" << std::endl;
	//return -100.000;
      }
    return fcorreta;
  }


  //
  //Method to correct energy loss due to bremsstrahlung
  //
  inline double fBrem(double brLinear, double e) {
    //
    // first parabola (for brLinear < threshold)
    // p0*x^2 + p1*x + p2
    // second parabola (for brLinear >= threshold)
    // ax^2 + bx + c, make y and y' the same in threshold
    // y = p0*threshold^2 + p1*threshold + p2
    // yprime = p1 + 2*p0*threshold
    // a = p3
    // b = yprime - 2*a*threshold
    // c = y - a*threshold^2 - b*threshold
    // final result is multiplied by cos(p5/br + p6*br + p7)

    // make NO correction if brLinear is invalid!
    if ( brLinear == 0 ) return e;
    //

    // this is for 31X
    if ( brLinear < 1.1 ) brLinear = 1.1;
    if ( brLinear > 8 ) brLinear = 8.0;

    // 3_1_X
    //double p0 = -0.05185;
    //double p1 = 0.1354;
    //double p2 = 0.9165;
    //double p3 = -0.0005626;
    //double p4 = 1.385;    

    // 3_10_0 implemented in CMSSW_4_2_X
    double p0 = -0.05289;
    double p1 = 0.1374;
    double p2 = 0.9141;
    double p3 = -0.000669;
    double p4 = 1.38;    

    double threshold = p4;
    
    double y = p0*threshold*threshold + p1*threshold + p2;
    double yprime = 2*p0*threshold + p1;
    double a = p3;
    double b = yprime - 2*a*threshold;
    double c = y - a*threshold*threshold - b*threshold;

    double fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;
    
    return e/fCorr;

  }

  inline double fullCorr(double et, double eta) {
    double fCorr = 0;

    // hybrid SC 
    // 3_10_0 implemented in CMSSW_4_2_X
    // fEtEta p0(ET)
    double c0 = 1.000;
    double c1 = -0.698;
    double c2 = 0;

    // fEtEta p1(ET)
    double c4 = 0.6605;
    double c5 = 8.825; 
    double c6 = 0.841;
 
    // final fitting
    double c7 = 1.081;  // curve point in eta distribution
    double c8 = 7.6;     // sharpness of the curve
    double c3 = -0.00181;


    double p0 = c0 + c1/(et + c2);
    double p1 = c4/(et + c5) + c6/(et*et);

    fCorr = p0 + p1*atan(c8*(c7 - fabs(eta))) + c3*fabs(eta);

    return et/fCorr;
  }

//--------------- E N D - C A P S -----------------

  inline double fBrem_ee(double brLinear, double e) {
    // first parabola (for brLinear < threshold)
    // p0*x^2 + p1*x + p2
    // second parabola (for brLinear >= threshold)
    // ax^2 + bx + c, make y and y' the same in threshold
    // y = p0*threshold^2 + p1*threshold + p2
    // yprime = p1 + 2*p0*threshold
    // a = p3
    // b = yprime - 2*a*threshold
    // c = y - a*threshold^2 - b*threshold
    // final result is multiplied by cos(p5/br + p6*br + p7)

    // make NO correction if brLinear is invalid!
    if ( brLinear == 0 ) return e;
    // make a flat correction if brLinear is too big (>9)

    // FM with preshower
    if ( brLinear > 6.5 ) brLinear = 6.5;
    if ( brLinear < 0.9 ) brLinear = 0.9;

    // ============= Fixed Matrix With Preshower SC
    // 3_10_0 implemented in CMSSW_4_2_X
    double p0 = -0.07945;
    double p1 = 0.1298;
    double p2 = 0.9147;
    double p3 = -0.001565;
    double p4 = 0.9;

    double threshold = p4;
    
    double y = p0*threshold*threshold + p1*threshold + p2;
    double yprime = 2*p0*threshold + p1;
    double a = p3;
    double b = yprime - 2*a*threshold;
    double c = y - a*threshold*threshold - b*threshold;

    double fCorr = 1;
    if ( brLinear < threshold ) 
      fCorr = p0*brLinear*brLinear + p1*brLinear + p2;
    else 
      fCorr = a*brLinear*brLinear + b*brLinear + c;
    
    return e/fCorr;

  }

  inline double fullCorr_ee(double et, double eta) {
    double fCorr = 0.;

    // 3_10_0 implemented in CMSSW_4_2_X
    double c0 = -3.516;
    double c1 = -2.362;

    double c2 = 2.151;
    double c3 = 1.572;

    double c4 = -0.336;
    double c5 = -0.2807;

    double c6 = 3.2;
    double c7 = 0.0;

    double p0 = c0 + c1/sqrt(et);
    double p1 = c2 + c3/sqrt(et);
    double p2 = c4 + c5/sqrt(et);
    double p3 = c6 + c7/sqrt(et);

    // make correction constant
    if ( fabs(eta) < 1.6 ) eta = 1.6;
    if ( fabs(eta) > 2.6 ) eta = 2.6;
    fCorr = p0 + p1*fabs(eta) + p2*eta*eta + p3/fabs(eta);
    
    if ( 1./fCorr > 1.2 ) return et/1.2;
    if ( 1./fCorr < 0.5 ) return et/0.5;
    return et/fCorr;
  }

}
#endif //INC_FCORRS

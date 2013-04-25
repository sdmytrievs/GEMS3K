//-------------------------------------------------------------------
// $Id$
//
// Stub implementation of TKinMet class for further development
//
// Copyright (C) 2013  D.Kulik, B.Thien, S.Dmitrieva
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//

#include <cmath>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "verror.h"
#include "s_kinmet.h"

//=============================================================================================
// TKinMet base class constructor (simulation of phase metastability and MWR kinetics)
// (c) DK March 2013
//
//=============================================================================================
/// generic constructor
//
TKinMet::TKinMet( const KinMetData *kmd ):
    KinProCode(kmd->KinProCod_),    KinModCode(kmd->KinModCod_),   KinSorpCode(kmd->KinSorpCod_),
    KinLinkCode(kmd->KinLinkCod_),  KinSizedCode(kmd->KinSizedCod_),  KinResCode(kmd->KinResCod_),
    NComp(kmd->NComp_), nlPh(kmd->nlPh_), nlPc(kmd->nlPc_), nPRk(kmd->nPRk_), nSkr(kmd->nSkr_),
    nrpC(kmd->nrpC_), naptC(kmd->naptC_), nAscC(kmd->nAscC_), // numpC(kmd->numpC_), iRes4(kmd->iRes4_),
    R_CONST(8.31451), T_k(kmd->T_k_), P_bar(kmd->P_bar_), kTau(kmd->kTau_), kdT(kmd->kdT_),
    IS(kmd->IS_), pH(kmd->pH_),  pe(kmd->pe_),  Eh(kmd->Eh_),  nPhi(kmd->nPh_), mPh(kmd->mPh_),
    vPh(kmd->vPh_), sAPh(kmd->sAPh_), LaPh(kmd->LaPh_), OmPh(kmd->OmPh_),
    sSAi(kmd->sSA_),  sgw(kmd->sgw_),  sgg(kmd->sgg_),  rX0(kmd->rX0_),  hX0(kmd->hX0_),
    sVp(kmd->sVp_), sGP(kmd->sGP_), nPul(kmd->nPul_), nPll(kmd->nPll_)
{
    // pointer assignments
    arfeSAr = kmd->arfeSAr_;   // [nPRk]
    arAscp = kmd->arAscp_;     // [nAscC]
    SM = kmd->SM_;             // [NComp]
    arDCC = kmd->arDCC_;       // [NComp]
    arPhXC = kmd->arPhXC_;     // [nlPh]
    arocPRk = kmd->arocPRk_;   // [nPRk]
    arxSKr = kmd->arxSKr_;     // [nSKr]
    arym = kmd->arym_;         // L
    arla = kmd->arla_;         // L
    arnx = kmd->arnx_;         // [NComp]
    arnxul = kmd->arnxul_;     // [NComp]  (DUL)
    arnxll = kmd->arnxll_;     // [NComp]  (DLL)
    arWx = kmd->arWx_;         // NComp
    arVol = kmd->arVol_;       // NComp

// 2D and 3D arrays
    arlPhc = NULL;
    arrpCon = NULL;
    arapCon = NULL;
//    arUmpCon = NULL;
    alloc_kinrtabs( );
    init_kinrtabs( kmd->arlPhc_, kmd->arrpCon_,  kmd->arapCon_ );

    /// allocation of work array of parameters and results for 'parallel reactions'
    arPRt = NULL;
    if(nPRk > 0 )
       arPRt = new TKinReact[nPRk];
//        alloc_arPRt();
    init_arPRt( );   // load data for parallel reactions

    // Work data and kinetic law calculation results
//    double spcfu[];    /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]
//    double spcfl[];    /// work array of coefficients for splitting nPul and nPll into nxul and nxll [NComp]

    double kTot = 0.;   /// Total rate constant (per m2 phase surface area)
    double rTot = 0.;   /// Current total MWR rate (mol/s)
    double vTot = 0.;   /// Total surface propagation velocity (nm/s)

    nPh = nPhi;
    sSA = sSAi;
    sSAcor = sSAi; /// Initialized corrected specific surface area (m2/g)
    sAph_c = sAPh; /// Initialized corrected surface area of the phase (m2/g)
    kdT_c = kdT;   /// Initialized corrected time step (s)

    T_k = 0.; // To trigger P-T recalculation after constructing the class instance
}

/// Destructor
TKinMet::~TKinMet()
{
  free_kinrtabs();
  if( arPRt )
      delete[] arPRt;
}

// Checks dimensions in order to re-allocate the class instance, if necessary
bool
TKinMet::testSizes( const KinMetData *kmd )
{
    bool status;
    status = (KinProCode == kmd->KinProCod_) && (KinModCode == kmd->KinModCod_) && (KinSorpCode == kmd->KinSorpCod_)
     && (KinLinkCode == kmd->KinLinkCod_) &&  (KinSizedCode == kmd->KinSizedCod_) && (KinResCode == kmd->KinResCod_)
     && (NComp == kmd->NComp_) && (nlPh == kmd->nlPh_) && (nlPc == kmd->nlPc_) && (nPRk == kmd->nPRk_)
               && (nSkr == kmd->nSkr_);
    return status;
}

/// allocates memory for TKinMet data
void
TKinMet::alloc_kinrtabs()
{
   long int j, s, lp, pr;

   if( nlPh && nlPc )
   {
     arlPhc = new double *[nlPh];
     for( lp=0; lp<nlPh; lp++)
     {
          arlPhc[lp] = new double[nlPc];
     }
   }

   if( nPRk && nrpC )
   {
     arrpCon = new double *[nPRk];
     for( pr=0; pr<nPRk; pr++)
     {
          arrpCon[pr] = new double[nrpC];
     }
   }

   if( nPRk && nSkr && naptC )
   {
        arapCon = new double **[nPRk];
        for(j=0; j<nPRk; j++)
        {
            arapCon[j]   = new double *[nSkr];
            for(s=0; s<nSkr; s++)
            {
                arapCon[j][s] = new double [naptC];
            }
        }
    }
}

/// returns 0 if o.k. or some arrays were not allocated.
/// returns -1 if error was encountered.
//
long int
TKinMet::init_kinrtabs( double *p_arlPhc, double *p_arrpCon,  double *p_arapCon  )
{
    long int j, i, s, lp, pr;

    if( arlPhc ) {

        for( lp=0; lp<nlPh; lp++)
            for( i=0; i<nlPc; i++)
                arlPhc[lp][i] = p_arlPhc[nlPc*lp+i];
    }
    if( arrpCon ) {

        for( pr=0; pr<nPRk; pr++)
            for( i=0; i<nrpC; i++)
                arrpCon[pr][i] = p_arrpCon[nrpC*pr+i];
    }
    if( arapCon ) {

        for( j=0; j<nPRk; j++)
            for( s=0; s<nSkr; s++)
                for( i=0; i<naptC; i++)
                    arapCon[j][s][i] = p_arapCon[naptC*nSkr*j+naptC*s+i];
    }
    return 0;
}

/// frees memory for TKinMet tables
//
void
TKinMet::free_kinrtabs()
{
    long int j, s, lp, pr;

    if( arlPhc )
    {
      for( lp=0; lp<nlPh; lp++)
      {
           delete[]arlPhc[lp];
      }
      delete[]arlPhc;
    }
    if( arrpCon )
    {
      for( pr=0; pr<nPRk; pr++)
      {
           delete[]arrpCon[pr];
      }
      delete[]arrpCon;
    }
    if( arapCon )
    {
         for(j=0; j<nPRk; j++)
         {
             for(s=0; s<nSkr; s++)
             {
                 delete[]arapCon[j][s];
             }
         }
         for(j=0; j<nPRk; j++)
         {
             delete[]arapCon[j];
         }
         delete[]arapCon;
    }
}

// creates the TKinReact array
//void alloc_arPRt()
//{
//    if( nPRk > 0 )
//    {
//        arPRt = new TKinReact[nPRk];
//    }
//}

/// Initializes TKinReact array and loads parameters into it
void
TKinMet::init_arPRt()
{
    long int xj;

    if( nPRk > 0 )
    {
        for(xj=0; xj<nPRk; xj++)
        {
            arPRt[xj].xPR = xj;   /// index of this parallel reaction
            arPRt[xj].nSa = nSkr; // number of species involved in parallel reactions
            arPRt[xj].ocPRk[0] = arocPRk[xj][0]; /// operation code for this kinetic parallel reaction affinity term
            arPRt[xj].ocPRk[1] = arocPRk[xj][1]; /// index of particle face (surface patch)
            arPRt[xj].xSKr = arxSKr;
            arPRt[xj].feSAr = arfeSAr[xj];
            arPRt[xj].rpCon = arrpCon[xj];
            arPRt[xj].apCon = arapCon[xj];

            // work data: unpacked rpCon[nrpC]
            if( nrpC >=4 )
            {
                arPRt[xj].kod = arPRt[xj].rpCon[0];  /// rate constant at standard temperature (mol/m2/s)
                arPRt[xj].kop = arPRt[xj].rpCon[1];  /// rate constant at standard temperature (mol/m2/s)
                arPRt[xj].Ap = arPRt[xj].rpCon[2];  /// Arrhenius parameter
                arPRt[xj].Ea = arPRt[xj].rpCon[3];  /// activation energy at st.temperature J/mol
            }
            else {
                arPRt[xj].kod = arPRt[xj].kop = arPRt[xj].Ap = arPRt[xj].Ea = 0.0;
            }
            if( nrpC >=8 )
            {
                arPRt[xj].bI = arPRt[xj].rpCon[4];
                arPRt[xj].bpH = arPRt[xj].rpCon[5];
                arPRt[xj].bpe = arPRt[xj].rpCon[6];
                arPRt[xj].bEh = arPRt[xj].rpCon[7];
            }
            else {
               arPRt[xj].bI = arPRt[xj].bpH = arPRt[xj].bpe = arPRt[xj].bEh = 0.0;
            }
            if( nrpC >=12 )
            {
                arPRt[xj].pPR = arPRt[xj].rpCon[8];
                arPRt[xj].qPR = arPRt[xj].rpCon[9];
                arPRt[xj].mPR = arPRt[xj].rpCon[10];
                arPRt[xj].uPR = arPRt[xj].rpCon[11];
            }
            else {
                arPRt[xj].pPR = arPRt[xj].qPR = arPRt[xj].mPR = arPRt[xj].uPR = 0.0;
            }
            if( nrpC > 12 )
            {
                arPRt[xj].OmEff = arPRt[xj].rpCon[12];
                arPRt[xj].nucRes = arPRt[xj].rpCon[13];
            }
            else {
                arPRt[xj].OmEff = 1.;
            }
//            arPRt[xj].Omg = OmPh; /// Input stability index non-log (d-less)

    // Results of rate term calculation
            arPRt[xj].arf = 1.;  // Arrhenius factor (temperature correction on kappa)
            arPRt[xj].cat = 1.;  // catalytic product term (f(prod(a))
            arPRt[xj].aft = 0.;  // affinity term (f(Omega))

            arPRt[xj].kd = arPRt[xj].kod;   // rate constant (involving all corrections) in mol/m2/s
            arPRt[xj].kp = arPRt[xj].kop;
            // check direction - for precipitation, arPRt[xj].kPR = arPRt[xj].kop;
            arPRt[xj].rPR = 0.;   // rate for this region (output) in mol/s
            arPRt[xj].rmol = 0.;   // rate for the whole face (output) in mol/s
//        arPRt[xj].velo,   // velocity of face growth (positive) or dissolution (negative) nm/s

        } // xj
    }
}

// frees memory for TKinReact array
//void free_arPRt()
//{
//    if( arPRt )
//        delete[]arPRt;
//}


//=============================================================================================
// TKinReact base class constructor (to keep data for parallel reaction regions)
// (c) DK March 2013
//
//=============================================================================================
// Generic constructor
//TKinReact::TKinReact( )
//{
//
//}

// Sets new specific surface area and other properties of the phase;
// also updates 'parallel reactions' area fractions
// returns false if these parameters in TKinMet instance did not change; true if they did.
//
bool
TKinMet::UpdateFSA( const double pAsk, const double pXFk, const double pFWGTk, const double pFVOLk,
                    const double pLaPh, const double pPULk, const double pPLLk, const double pYOFk,
                    const double pICa, const double ppHa, const double ppea, const double pEha )
{
    long int i;
    bool status = false;
    if( sSA != pAsk || nPh != pXFk || mPh != pFWGTk || vPh != pFVOLk || nPul != pPULk || nPll != pPLLk
        || LaPh != pLaPh || IS != pICa || pH != ppHa || pe != ppea || Eh != pEha || sGP != pYOFk*mPh )
        status = true;
    sSA = pAsk;
    nPh = pXFk;
    mPh = pFWGTk;
cout << " !!! mPh: " << mPh << endl;
    vPh = pFVOLk;
    nPul = pPULk;
    nPll = pPLLk;
    LaPh = pLaPh;
    sGP = pYOFk*mPh;
    IS = pICa;
    pH = ppHa;
    pe = ppea;
    Eh = pEha;
    sAPh = sSA*mPh;             // current surface area of this phase, m2
    OmPh = pow( 10., LaPh );    // phase stability index (activity scale) 10^LaPh_

    for( i = 0; i < nPRk; i++ )
    {
       if( arPRt[i].feSAr != arfeSAr[i] )
           status = true;
       arPRt[i].feSAr = arfeSAr[i];
    }
    return status;
}

// Returns (modified) specific surface area of the phase;
// and gets (modified) 'parallel reactions' area fractions
//
double
TKinMet::GetModFSA ( double& pPULk, double& pPLLk  )
{
    long int i;
    for( i = 0; i < nPRk; i++ )
    {
       arPRt[i].feSAr = arfeSAr[i];
    }
    pPULk = nPul;
    pPLLk = nPll;
    return sSAcor;
}

// Updates temperature to T_K and pressure to P_BAR;
// calculates Arrhenius factors and temperature-corrected rate constants in all PR regions.
//
long int
TKinMet::UpdatePT ( const double T_K, const double P_BAR )
{
    double Arf=1., kd=0., kp=0.;
    double RT = R_CONST*T_K;
    long int i;

    T_k = T_K; P_bar = P_BAR;

    for( i=0; i<nPRk; i++ )
    {
        Arf = arPRt[i].Ap * exp(-(arPRt[i].Ea)/RT);
        if( arPRt[i].Ap != 0.0 )
            arPRt[i].arf = Arf;
        else
            Arf = 1.0;
        kd = arPRt[i].kod * Arf;
        kp = arPRt[i].kop * Arf;
        arPRt[i].kd = kd;
        arPRt[i].kp = kp;
    }
    return 0;
}

// sets new time Tau and time step dTau
// returns false if neither kTau nor kdT changed; true otherwise
//
bool
TKinMet::UpdateTime( const double Tau, const double dTau )
{
    bool status = false;
    if( Tau != kTau || dTau != kdT )
        status = true;
    kTau = Tau;
    kdT = dTau;
    return status;
}

//long int
//TKinMet::PTparam()
//{
//    int iRet = 0;
//    iRet = UpdatePT(  );
//    return iRet;
//};


// calculates the rate per m2 for r-th (xPR-th) parallel reaction
double
TKinMet::PRrateCon( TKinReact &k, const long int r )
{
   long int xj, j, atopc, facex;
   double atp, ajp, aj;

cout << "kTau: " << kTau << " kd: " << k.kd << " kp: " << k.kp << " Omg: " << OmPh <<
        " nPh: " << nPh << " mPh: " << mPh << endl;

if( k.xPR != r )     // index of this parallel reaction
        cout << k.xPR << " <-> " << r << " mismatch" << endl;
   atopc = k.ocPRk[0]; // operation code for this kinetic parallel reaction affinity term
   facex = k.ocPRk[1]; // particle face index

   // activity (catalysis) product term (f(prod(a))
   k.cat = 1.;
   for( xj=0; xj < k.nSa; xj++ )
   {
       j = k.xSKr[xj];
       aj = arla[j];
       if( k.apCon[xj][0] )
       {    ajp = pow( aj, k.apCon[xj][0] );
            // may need extension in future
       }
       else
           ajp = 1.;
       k.cat *= ajp;
   }
   if( k.pPR )
       k.cat = pow( k.cat, k.pPR );
   if( k.bI )
       k.cat *= pow( IS, k.bI );
   if( k.bpH )
       k.cat *= pow( pH, k.bpH );
   if( k.bpe )
       k.cat *= pow( pe, k.bpe );
   if( k.bEh )
       k.cat *= pow( Eh, k.bEh );

   // affinity term (f(Omega))
   k.aft = 0.; atp = 0.;
   switch( atopc )    // selecting a proper affinity term
   {
    default:
    case ATOP_CLASSIC_: // = 0,     classic TST affinity term (see .../Doc/Descriptions/KinetParams.pdf)
       k.aft = 1. - pow( OmPh, k.qPR );
       if( k.mPR )
           k.aft = pow( k.aft, k.mPR );
       break;
    case ATOP_CLASSIC_REV_: // = 1, classic TST affinity term, reversed
       k.aft = pow( OmPh, k.qPR ) - 1. - k.uPR;
       if( k.mPR )
           k.aft = pow( k.aft, k.mPR );
       break;
    case ATOP_SCHOTT_: // = 2,      Schott et al. 2012 fig. 1e
       if( OmPh )
           k.aft = exp( -k.uPR/OmPh );
       else
           k.aft = 0.;
       break;
    case ATOP_HELLMANN_: // = 3,    Hellmann Tisserand 2006 eq 9
       if( OmPh )
       {
          atp = pow( k.qPR*log( OmPh ), k.uPR );
          k.aft = 1 - exp( -atp );
       }
       else {
          atp = 0.; k.aft = 0.;
       }
       break;
    case ATOP_TENG1_: // = 4,       Teng et al. 2000, eq 13
       if( OmPh )
           atp = log( OmPh );
       else
           atp = 0.;
       k.aft = k.uPR * ( OmPh - 1. ) * atp;
       break;
    case ATOP_TENG2_: // = 5,       Teng et al. 2000, Fig 6
       k.aft = pow( OmPh, k.mPR );
       break;
    case ATOP_FRITZ_: // = 6        Fritz et al. 2009, eq 6 nucleation and growth
       atp = OmPh - k.OmEff;
       k.aft = pow( atp, k.mPR );
//   nucRes
       break;
   }

   // Calculating rate for this region (output) in mol/m2/s
   if( LaPh < -0.001 ) // dissolution  (needs more flexible check based on Fa stability criterion!
   {   // kd     dissolution rate constant (corrected for T) in mol/m2/s
       k.rPR = k.kd * k.cat * k.aft;
   }
   else if( LaPh > 0.001 ) {  // precipitation rate constant (corrected for T) in mol/m2/s
       k.rPR = k.kp * k.cat * k.aft;
   }
   else {  // equilibrium
       k.rPR = 0.;
   }
   k.rPR *= k.feSAr; // Theta: fraction of surface area of the solid related to this parallel reaction

   return k.rPR;
}
//   rmol,   // rate for the whole face (output) in mol//m2/s    TBD
//   velo,   // velocity of face growth (positive) or dissolution (negative) nm/s
//       kod,  /// dissolution rate constant at standard temperature (mol/m2/s)
//       kop,  /// precipitation rate constant at standard temperature (mol/m2/s)
//       Ap,   /// Arrhenius parameter
//       Ea,   /// activation energy at st.temperature J/mol

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Implementation of TMWReaKin class (TST kinetics)
//
bool
TMWReaKin::PTparam( const double TK, const double P )
{
    bool iRet = false;
    if( TK < 273. || TK > 5273. || P < 0. || P > 1e6 )
        iRet = true;  // error
    if( fabs( TK - T_k ) > 0.1 || fabs( P - P_bar ) > 1e-5 )
    {
       iRet = UpdatePT( TK, P );
    }
    return iRet;
}


bool
TMWReaKin::RateInit( )
{   
    double RT = R_CONST*T_k;
    long int r;

    kTot = 0.;
    for( r=0; r<nPRk; r++ )
    {
       kTot += PRrateCon( arPRt[r], r ); // adds the rate constant (mol/m2/s) for r-th parallel reaction
    }
    rTot = kTot * sSAcor; // overall rate (mol/s)
    nPhi = nPh;
    sAPh = sSA * mPh;   // surface of the phase
    rTot = kTot * sAPh; // overall rate (mol/s)

cout << " init kTot: " << kTot << " rTot: " << rTot << " sSAcor: " << sSAcor << " sAPh: " << sAPh << endl;
 // Check initial rates and limit them to realistic values
    return false;
}

bool
TMWReaKin::RateMod( )
{
    double RT = R_CONST*T_k;
    long int r;

    kTot = 0.;
    for( r=0; r<nPRk; r++ )
    {
       kTot += PRrateCon( arPRt[r], r ); // adds the rate constant (mol/m2/s) for r-th parallel reaction
    }
    sSAcor = sSAi * pow( nPh/nPhi, 1./3. ); // primitive correction for specific surface area
    sAPh = sSA * mPh;   // surface of the phase
    rTot = kTot * sAPh; // overall rate (mol/s)

cout << " cur kTot: " << kTot << " rTot: " << rTot << " sSAcor: " << sSAcor << " sAPh: " << sAPh << endl;
    return false;
}

// Sets new metastability constraints based on updated kinetic rates
//
bool
TMWReaKin::SetMetCon( )
{
   double dnPh = 0., dn_dt = 0.;

   dnPh = - kdT * rTot; // change in phase amount over dt
   dn_dt = - kdT;  // minus total rate mol/s

   // First calculate phase constraints
   if( LaPh < -0.001 ) // dissolution  (needs more flexible check based on Fa stability criterion!
   {
       nPll = nPh + dnPh;
       if( nPul < nPll )
           nPul = nPll;
   }
   else if( LaPh > 0.001 ) {  // precipitation rate constant (corrected for T) in mol/m2/s
       nPul = nPh + dnPh;
       if( nPll > nPul )
           nPll = nPul;
    }
    else {  // equilibrium  - needs to be checked!
       nPul = nPh + fabs( dnPh );
       nPll = nPh - fabs( dnPh );
    }

    if( NComp > 1 )
       return false;  // DC constraints will be set in SplitMod() or SplitInit()
    // setting DC metastability constraints
    arnxul[0] = nPul;
    arnxll[0] = nPll;
    return false;
}


bool
TMWReaKin::SplitInit( )
{
    return false;
}

bool
TMWReaKin::SplitMod( )
{
    return false;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Implementation of TUptakeKin class (uptake kinetics)
//
TUptakeKin::TUptakeKin( KinMetData *kmd, long int p_numpC, double *p_arUmpCon ):TKinMet( kmd )
{
    numpC = p_numpC;
    arUmpCon = NULL;

    alloc_upttabs();
    init_upttabs( p_arUmpCon );
};

/// Destructor
TUptakeKin::~TUptakeKin()
{
  free_upttabs();
}

/// allocates memory for TUptakeKin data
void
TUptakeKin::alloc_upttabs()
{
   long int j;

   if( NComp && numpC )
   {
     arUmpCon = new double *[NComp];
     for( j=0; j<NComp; j++)
     {
          arUmpCon[j] = new double[numpC];
     }
   }

}

/// returns 0 if o.k. or some arrays were not allocated.
/// returns -1 if error was encountered.
//
long int
TUptakeKin::init_upttabs( double *p_arUmpCon  )
{
    long int j, i;

    if( arUmpCon ) {

        for( j=0; j<NComp; j++)
            for( i=0; i<numpC; i++)
                arUmpCon[j][i] = p_arUmpCon[numpC*j+i];
    }
    return 0;
}

/// frees memory for TKinMet tables
//
void
TUptakeKin::free_upttabs()
{
    long int j;

    if( arUmpCon )
    {
      for( j=0; j<NComp; j++)
      {
           delete[]arUmpCon[j];
      }
      delete[]arUmpCon;
    }
}


bool TUptakeKin::PTparam( const double TK, const double P )
{
    int iRet = 0;
    iRet = UpdatePT( TK, P );
    return iRet;
};


// Initializes uptake rates
bool
TUptakeKin::UptakeInit()
{
    return false;
}


// Calculates uptake rates
bool
TUptakeKin::UptakeMod()
{

    return false;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



/*
// -----------------------------------------------------------------------------
// Implementation of derived class main functionality
//
bool PTparam()
{
    return 0;
};

bool RateMod()
{
    return 0;
};

bool SplitMod()
{
    return 0;
};

bool SplitInit()
{
    return 0;
};

bool UptakeMod()
{
    return 0;
};

bool UptakeInit()
{
    return 0;
};

bool SetMetCon()
{
    return 0;
};


*/


//--------------------- End of s_kinmet.cpp ---------------------------


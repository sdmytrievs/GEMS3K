//-------------------------------------------------------------------
/// \file solmodengine.h
///
/// Declaration of the SolModEngine class - c++ API for phase models
/// Decorator for TSolMod and derived classes implementing built-in models
/// of mixing in fluid, liquid, aqueous, and solid-solution phases
//
// Copyright (c) 2023 S.Dmytriyeva
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
#ifndef SOLMODENGINE_H
#define SOLMODENGINE_H

#include <map>
#include "s_solmod.h"

/// Addition parameters for specific models,
/// could be add to main SolutionData
struct AddSolutionData {
    double *arZ;
    double *arM;
    double *ardenW;
    double *arepsW;
    double *arG0;
    double *arFWGT;
    double *arX;
};

/// Decorator for TSolMod and derived classes implementing built-in models of mixing
/// in fluid, liquid, aqueous, and solid-solution phases
class SolModEngine
{
   friend class SolModFactory;

public:

    /// The constructor
    SolModEngine(long int k, long int jb, SolutionData& sd, const AddSolutionData& addsd);
    /// The empty constructor
    explicit SolModEngine(long int k=0, long int jb=0, const std::string& phase_name="undefined");

    /// Code of the mixing model
    char model_code() const {
        return mod_code;
    }


    /// Call for calculation of temperature and pressure correction for the phase
    void SolModPTParams();
    // What is the difference with UpdatePT() ?
    // UpdatePT() only set new value of T an P, real changes only into SolModFactory UpdateThermoData
    // I done this function private
    /// Call for calculation of activity coefficients of species endmembers
    void SolModActivityCoeffs();
    /// Call for calculation of bulk phase excess mixing properties
    std::map<std::string, double> SolModExcessProps();
    /// Call for calculation of bulk phase ideal mixing properties
    std::map<std::string, double> SolModIdealProps();
    /// Call for retrieving bulk phase Darken quadratic terms  (not yet implemented?)
    std::map<std::string, double> SolModDarkenProps();
    /// Call for retrieving bulk phase standard state terms
    std::map<std::string, double> SolModStandProps();

    /// high-level method to retrieve pure fluid fugacities
    /// Implemented for TPRSVcalc, TCGFcalc, TSRKcalc, TPR78calc, TCORKcalc, TSTPcalc
    // never call in gems, do we need here?
    // PureSpecies();

    // Get functions ("getters")

    /// Get number of species (endmembers) in this solution phase
    long int Get_SpeciesNumber() {
        return dc_num;
    }
    /// Get names of species (endmembers) in this phase as a list of strings
    std::vector<std::string> Get_SpeciesNames() {
       return dc_names;
    }

    /// Copy ln of activity coefficients of species (end members) into provided array lngamma
    void Get_lnActivityCoeffs(double* lngamma);
    /// Get ln of activity coefficients of chemical species (end members) in a component map
    std::map<std::string, double> GetlnActivityCoeffs();

    /// NB: for an endmember (species) in solution phase:
    ///     ln_ActivityCoeff = lnConfTerm + lnRecipTerm + lnExcessTerm + lnDQFTerm
    ///       Some or all of terms on r.h.s. can be zeros
    /// Access methods to all components of ln_ActivityCoeff below

    /// Copy configurational terms adding to overall activity into provided array
    void Get_lnConfTerms(double* lnGamConf);
    /// Get configurational terms adding to overall activity into component map
    std::map<std::string, double> GetlnConfTerms();

    // lnRecipTerm
    // Used in models TBerman, TCEFmod, TMBWmod
    // ( SM_BERMAN-'B', SM_CEF-'$', SM_MBW-'#' )
    /// Copy reciprocal terms adding to overall activity coefficients into provided array
    /// Implemented for the mixing model 'B', '$', '#'
    void Get_lnRecipTerms(double* lnGamRecip);
    /// Get reciprocal terms adding to overall activity coefficients into component map
    /// Implemented for the mixing model 'B', '$', '#'
    std::map<std::string, double> GetlnRecipTerms();

    // lnExcessTerm
    // Used in models TBerman, TCEFmod, TMBWmod
    // ( SM_BERMAN-'B', SM_CEF-'$', SM_MBW-'#' )
    /// Copy excess energy terms adding to overall activity coefficients into provided array
    /// Implemented for the mixing model 'B', '$', '#'
    void Get_lnExcessTerms(double* lnGamEx);
    /// Get excess energy terms adding to overall activity coefficients into component map
    /// Implemented for the mixing model 'B', '$', '#'
    std::map<std::string, double> GetlnExcessTerms();

    // lnDQFTerm
    // Used in model TSubregular
    // ( SM_MARGB-'M' )
    /// Copy DQF terms adding to overall activity coefficients into provided array
    /// Implemented for the mixing model 'M'
    void Get_lnDQFTerms(double* lnGamDQF);
    /// Get DQF terms adding to overall activity coefficients into component map
    /// Implemented for the mixing model 'M'
    std::map<std::string, double> GetlnDQFTerms();

    /// Copy increments to molar G0 values of DCs from pure gas fugacities or DQF terms
    /// into provided array
    void Get_G0Increments(double* aGEX);
    /// Get increments to molar G0 values of DCs from pure gas fugacities or DQF terms
    /// into component map
    std::map<std::string, double> GetG0Increments();

    // aVol
    // Used in models TPRSVcalc, TCGFcalc, TSRKcalc, TPR78calc, TCORKcalc, TSTPcalc
    // ( SM_PRFLUID-'P', SM_CGFLUID-'F', SM_SRFLUID-'E', SM_PR78FL-'7', SM_CORKFL-'8', SM_STFLUID -'6' )
    /// Copy molar volumes of species into provided array
    /// Implemented for the mixing model 'P', 'F', 'E', '7', '8', '6'
    void Get_MolarVolumes(double* aVol);
    /// Get molar volumes of species into component map
    /// Implemented for the mixing model 'P', 'F', 'E', '7', '8', '6'
    std::map<std::string, double> GetMolarVolumes();

    /// Get phase volume, cm3/mol
    /// Implemented for the mixing model 'P', 'F', 'E', '7'
    double GetPhaseVolume();

    // Pparc
    // Used in models TPRSVcalc, TCGFcalc, TSRKcalc, TPR78calc, TCORKcalc, TSTPcalc
    // ( SM_PRFLUID-'P', SM_CGFLUID-'F', SM_SRFLUID-'E', SM_PR78FL-'7', SM_CORKFL-'8', SM_STFLUID -'6' )
    /// Copy partial pressures or fugacities of pure DC into provided array aPparc
    /// Implemented for the mixing model 'P', 'F', 'E', '7', '8', '6'
    void Get_PartialPressures(double* aPparc);
    /// Get partial pressures or fugacities of pure DC into component map
    /// Implemented for the mixing model 'P', 'F', 'E', '7', '8', '6'
    std::map<std::string, double> GetPartialPressures();

    // Set functions

    /// Set species (end member) mole fractions from a provided array aWx -> dc_num
    void Set_MoleFractions(double* aWx);
    /// Set species (end member) mole fractions from a provided map (dictionary) val_map
    void SetMoleFractions(const std::map<std::string, double>& val_map, double def_val=0.);

    // Y_m
    // Used in models TPitzer, TSIT, TEUNIQUAC, TKarpov, TDebyeHueckel, TLimitingLaw,
    // TShvarov, THelgeson, TDavies, TELVIS
    // ( SM_AQPITZ - 'Z', SM_AQSIT - 'S', SM_AQEXUQ - 'Q', SM_AQDH3 - '3', SM_AQDH2 - '2', SM_AQDH1 - '1',
    //   SM_AQDHS - 'Y', SM_AQDHH - 'H', SM_AQDAV - 'D', SM_AQELVIS 'J')
    /// Set molalities of aqueous species and sorbates from provided array  aM -> dc_num
    /// Used for the mixing models 'Z','S','Q','3','2','1','Y','H','D'
    void Set_SpeciesMolality(double* aM);
    /// Set molalities of aqueous species and sorbates from a provided map
    /// Used for the mixing models 'Z','S','Q','3','2','1','Y','H','D'
    void SetSpeciesMolality(const std::map<std::string, double>& val_map, double def_val=0.);

    // X (mole amounts), mass (FWGT in grams)
    // Used in CG fluid model TCGFcalc (SM_CGFLUID - 'F')
    /// Set chemical species (DC) quantities (mole amounts) at eqstate from provided array  aX -> dc_num
    /// Used for the mixing model 'F'
    void Set_SpeciesAmounts(double* aX);
    /// Set DC quantities at eqstate from component map
    /// Used for the mixing model 'F'
    void SetSpeciesAmounts(const std::map<std::string, double>& val_map, double def_val=0.);

    /// Set phase (carrier) mass, g
    /// Used for the mixing model 'F'
    void Set_PhaseMass(double aFWGT);

    /// Writing the input structure of the current phase to JSON format file
    void to_json_file(const std::string& path) const
    {
        if(solmod_task.get()) {
            solmod_task->to_json_file(path);
        }
    }
    /// Writing the input structure for the current phase to JSON stream
    void to_json(std::iostream& ff) const
    {
        if(solmod_task.get()) {
            solmod_task->to_json_stream(ff);
        }
    }

    /// Trace writing of the current phase arrays to a key-value format file
    void to_text_file(const std::string& path, bool append=false)
    {
        if(solmod_task) {
            solmod_task->to_text_file(path, append);
        }
    }

protected:

    /// Model name (posible add to TSolMod structure for test output)
    std::string model_name;
    /// Code of the mixing model
    char mod_code;

    /// Name of phase
    std::string phase_name;
    /// Index of phase in IPM problem
    long int phase_ndx;
    /// Index first of DCs included into phase
    long int dc_ndx;
    /// Number of DCs included into phase
    long int dc_num;
    /// Names of DCs included into phase
    std::vector<std::string> dc_names;

    /// Mole fractions Wx of DC in multi-component phases -> dc_num
    double *arWx = nullptr;
    /// Molalities of aqueous species and sorbates
    double *arM = nullptr;
    /// phase (carrier) masses, g
    double *arFWGT = nullptr;
    ///  DC quantities at eqstate x_j, moles
    double *arX = nullptr;

    /// Output: activity coefficients of species (end members)
    double *arlnGam = nullptr;
    /// DQF terms adding to overall activity coefficients
    double *arlnDQFt = nullptr;
    /// reciprocal terms adding to overall activity coefficients
    double *arlnRcpt = nullptr;
    /// excess energy terms adding to overall activity coefficients
    double *arlnExet = nullptr;
    /// configurational terms adding to overall activity
    double *arlnCnft = nullptr;
    /// Pure-species fugacities, G0 increment terms
    double *arGEX = nullptr;
    /// molar volumes of end-members (species) cm3/mol
    double *arVol = nullptr;
    /// phase volumes, cm3/mol
    double *aphVOL = nullptr;
    /// Partial pressures
    double *arPparc = nullptr;

    /// TSolMod description
    std::shared_ptr<TSolMod> solmod_task;

    void SolMod_create(SolutionData &sd, const AddSolutionData &addsd);
    bool check_mode(char ModCode);
    std::map<std::string, double> property2map(double* dcs_size_array);
    void map2property(const std::map<std::string, double>& dsc_name_map,
                      double* dcs_size_array, double def_value);

    /// Updates P and T in the mixing and activity model for this phase
    ///  (if those have changed)
    void UpdatePT( double T_k, double P_bar)
    {
        if(solmod_task) {
            solmod_task->UpdatePT(T_k, P_bar);
        }
    }
};

#endif // SOLMODENGINE_H

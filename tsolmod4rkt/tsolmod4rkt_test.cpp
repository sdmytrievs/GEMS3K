//--------------------------------------------------------------------
//
/// Demo test of usage of the standalone SolModEngine.
//
// Copyright (C) S.Dmytriyeva, D.Kulik
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>
//-------------------------------------------------------------------

// Can this ifdef be moved somewhere so that it is not needed in examples?
#ifdef OVERFLOW_EXCEPT
#ifdef __linux__
#include <cfenv>
#elif _MSC_VER
#include <float.h>
#else
#include <cfenv>
#endif
#endif

#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <string>
#include <iomanip>
#include <memory>
#include "jsonconfig.h"
#include "v_service.h"
#include "solmodfactory.h"

// Thermo-time-in/series1-dat.lst
//The simplest case: data exchange using disk files only
int main( int argc, char* argv[] )
{

    // Can this ifdef be moved somewhere so that it is not needed in examples?
    // It can be removed, but now better left check for different model cases or we make addition
    // application to test other models posible with  check threading and memory leaks
#if  defined(OVERFLOW_EXCEPT)
#ifdef __linux__
    feenableexcept (FE_DIVBYZERO|FE_OVERFLOW|FE_UNDERFLOW);
#elif _MSC_VER
    _clearfp();
    _controlfp(_controlfp(0, 0) & ~(_EM_INVALID | _EM_ZERODIVIDE | _EM_OVERFLOW),
               _MCW_EM);
#else

#endif
#endif

    // Read gems3k logger config file                   (is this needed in the example?)
    // this is example to config loggers
    gemsSettings().gems3k_update_loggers(true, "", 3);

    try{
        // Analyzing command line arguments
        // Defaults
        std::string after_reading = "_AfterReading.txt";
        std::string input_system_file_list_name = "Thermo-time-in/series1-dat.lst";

        // Get the list file name (of DCH, IPM and DBR input files) for a GEMS3K fileset
        if (argc >= 2 )
            input_system_file_list_name = argv[1];

        // Initialize SolModFactory from the GEMS3K file set
        SolModFactory task(input_system_file_list_name);

        // Optional: Check the data read for SolModFactory initialization
        std::string task_data_file_name = input_system_file_list_name + after_reading;
        task.to_text_file( task_data_file_name );

        // Getting the number of solution phases
        auto PhSolNumber = task.Get_SolPhasesNumber();

        // Getting the list of names of solution phases
        auto PhSolNames = task.Get_SolPhasesNames();

        std::cout << "\nTask: " << task_data_file_name << "\nPhSolNumber: " << PhSolNumber
                  << "\nPhSolNames: ";
        for_each(PhSolNames.begin(), PhSolNames.end(), [](const auto& element)
        { //printing using " " separator
            std::cout << "'" << element << "' ";
        });
        std::cout << std::endl;

        // Example of model calculations for the task (Thermo-time-in/series1-dat.lst)
        // Two feldspar phases co-existing on binodal solvus
        // In future this could be a test input from dbr file or directly from json string

        // Getting SolModEngine for a feldspar phase 1 by name (and print in C style)
        std::vector<double> x1v = { 0.674, 7e-08, 0.326 };
        std::vector<double> lnGamma1v(3, 0.);
        auto phase1 = task.SolPhase("Alkali feldspar");
        
        // Checking the mixing model in this phase
        auto phase_name_p1 = phase1.Get_SolPhaseName();
        auto mod_code_p1 = phase1.Get_MixModelCode();
        auto mod_type_p1 = phase1.Get_MixModelType();
        // Getting the number of endmembers
        auto n_species_p1 = phase1.Get_SpeciesNumber();
        std::cout << "\nPhase1: name '" << phase_name_p1 << "', mixing/activity model type '"
                  << mod_type_p1 << "', mixing/activity model code '" <<  mod_code_p1
                  << "', number of endmembers " << n_species_p1 << std::endl;

        // Getting the list of endmember names
        auto em_names_p1 = phase1.Get_SpeciesNames();

        // Setting phase composition
        phase1.Set_MoleFractions(x1v.data());
        
        // Getting phase composition in C++ style (e.g. for checking)
        double* x_ph1 = new double[n_species_p1];
        phase1.Get_MoleFractions(x_ph1);
        // Printing input phase composition in C style
        for(size_t j=0; j<n_species_p1; j++)
        {
            std::cout << "   '" << em_names_p1[j] << "': x= " << x_ph1[j] << std::endl;
        }
        delete[] x_ph1;
        
        // Calculating activity coefficients of end members
        phase1.SolModActivityCoeffs();
        
        // Writing results to a text file
        phase1.to_text_file("solmod_act_coef.txt", true);

        // Get activity coefficients and print them in C style
        phase1.Get_lnActivityCoeffs(lnGamma1v.data());
        std::cout << "Calculated activity coefficients: " << std::endl;
        for(size_t j=0; j<n_species_p1; j++)
        {
            std::cout << "   '" << em_names_p1[j] << "': ln(gamma)= " << lnGamma1v[j] <<
                         ", gamma= " << exp(lnGamma1v[j]) << std::endl;
        }
        
        // Getting SolModEngine for a felsdpar phase 2 by index
        //  (and get and print results in dict style in similar order as for phase1
        std::map<std::string, double> x2m = {
            {"Albite", 0.187},
            {"Anorthite", 3.5e-09},
            {"Sanidine", 0.813}};
        std::map<std::string, double> lnGamma2m;

        size_t p2index = 2;
        auto phase2 = task.Sol_Phase( p2index );

        // Checking the mixing model in this phase
        auto phase_name_p2 = phase2.Get_SolPhaseName();
        auto mod_code_p2 = phase2.Get_MixModelCode();
        auto mod_type_p2 = phase2.Get_MixModelType();
        // Getting the number of endmembers
        auto n_species_p2 = phase2.Get_SpeciesNumber();
        std::cout << "\nPhase2: name '" << phase_name_p2 << "', mixing/activity model type '"
                  << mod_type_p2 << "', mixing/activity model code '" <<  mod_code_p2
                  << "', number of endmembers " << n_species_p2 << std::endl;

        // Setting phase composition
        phase2.SetMoleFractions(x2m);

        auto x_ph2 = phase2.GetMoleFractions();
        // Printing input phase composition in dict style
        for(const auto& item: x_ph2 ) {
            std::cout << "   '" << item.first << "': x= " << item.second << std::endl;
        }

        // Calculating activity coefficients of end members
        phase2.SolModActivityCoeffs();

        // Writing results to a text file
        phase2.to_text_file("solmod_act_coef.txt", true);

        // Get activity coefficients and print them in dict style
        lnGamma2m = phase2.GetlnActivityCoeffs();
        std::cout << "Calculated activity coefficients: " << std::endl;
        for(const auto& item: lnGamma2m ) {
            std::cout << "   '" << item.first << "': ln(gamma)= " << item.second <<
                         ", gamma= " << exp(item.second) << std::endl;
        }


        // Calculate (a dict) of ideal properties of mixing in the phase2
        auto map_ideal = phase2.SolModIdealProps();
        phase2.to_text_file("solmod_act_coef.txt", true);
        std::cout << "\nIdeal properties of mixing:\n";
        for(const auto& item: map_ideal ) {
            std::cout << "   '" << item.first << "': " << item.second << std::endl;
        }

        // Calculate (a dict) of excess properties of mixing in the phase2
        auto map_excess = phase2.SolModExcessProps();
        phase2.to_text_file("solmod_act_coef.txt", true);
        std::cout << "\nExcess properties of mixing:\n";
        for(const auto& item: map_excess ) {
            std::cout << "   '" << item.first << "': " << item.second << std::endl;
        }


        return 0;
    }
    catch(TError& err)
    {
        std::cout  << err.title << err.mess << std::endl;
    }
    catch(std::exception& e)
    {
        std::cout  << "std::exception: " << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout  << "unknown exception" << std::endl;
    }
    return -1;
}


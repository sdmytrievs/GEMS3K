//-------------------------------------------------------------------
// $Id$
//
// Debugging version of a finite-difference 1D advection-diffusion
// mass transport model supplied by Dr. Frieder Enzmann (Uni Mainz)
// coupled with GEMIPM2K module for calculation of chemical equilibria
//
// Direct access to the TNodeArray class for storing all data for nodes
//
// Copyright (C) 2005,2007 S.Dmytriyeva, F.Enzmann, D.Kulik
//
//-------------------------------------------------------------------

#include "m_gem2mt.h"
#include "GEMS3K/jsonconfig.h"

static int task_from_file(const std::string& gem2mt_file, const std::string& ipm_lst, const std::string& dbr_lst);
static int task_A(const std::string& ipm_lst, const std::string& dbr_lst);


//---------------------------------------------------------------------------
// Test of 1D advection (finite difference method provided by Dr. F.Enzmann,
// Uni Mainz) coupled with GEMIPM2K kernel (PSI) using the TNodeArray class.
// Finite difference calculations split over independent components
// (through bulk composition of aqueous phase).
// Experiments with smoothing terms on assigning differences to bulk composition
// of nodes

// "TestVTK/GEM2MT-task.json" "TestVTK/CaWBoundC-dat.lst" "TestVTK/CaWBoundC-dbr.lst" "TestVTK/VTK"

int main( int argc, char* argv[] )
{

#ifndef USE_NLOHMANNJSON
    std::string gem2mt_in1 = "TestAD/CalColumnAD.dat";
    std::string ipm_lst = "TestAD/CalcColumn-dat.lst";
    std::string dbr_lst = "TestAD/CalcColumn-dbr.lst";
#else
    std::string gem2mt_in1 = "TestAD1/CalcColumnAD.json";
    std::string ipm_lst = "TestAD1/CalcColumn-dat.lst";
    std::string dbr_lst = "TestAD1/CalcColumn-dbr.lst";
#endif


    // from argv
    if (argc >= 2 )
        gem2mt_in1 = argv[1];
    if (argc >= 3 )
        ipm_lst = argv[2];
    if (argc >= 4 )
        dbr_lst = argv[3];

    GemsSettings().gems3k_update_loggers( true, "gems3k_logger.log", spdlog::level::info);


    try{
        //return task_from_file(gem2mt_in1, ipm_lst, dbr_lst);
        return task_A(ipm_lst, dbr_lst);
    }
    catch(TError& err) {
        TNode::ipmlog_file->error("Error {} : {}", err.title, err.mess);
        return 1;
    }

    return 0;
}

//---------------------------------------------------------------------------
// Test of 1D advection (finite difference method provided by Dr. F.Enzmann,
// Uni Mainz) coupled with GEMIPM2K kernel (PSI) using the TNodeArray class.
// Finite difference calculations split over independent components
// (through bulk composition of aqueous phase).
// Experiments with smoothing terms on assigning differences to bulk composition
// of nodes
// "TestVTK/GEM2MT-task.json" "TestVTK/CaWBoundC-dat.lst" "TestVTK/CaWBoundC-dbr.lst"
int task_from_file(const std::string& gem2mt_file, const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(gem2mt_file.empty() || ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT(0) );
    TGEM2MT::pm = mt_task.get();

    // Here we read the GEM2MT structure, prepared from GEMS or by hand
    if(TGEM2MT::pm->ReadTask(gem2mt_file, "")) {
        return 1;  // error reading files
    }

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    // TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}


int task_A(const std::string& ipm_lst, const std::string& dbr_lst)
{
    if(ipm_lst.empty() || dbr_lst.empty()) {
        Error( "Start task", "No inital files");
    }

    // The NodeArray must be allocated here
    std::shared_ptr<TGEM2MT> mt_task( new TGEM2MT('A', 201) );
    TGEM2MT::pm = mt_task.get();

    // Set up sizes, flags and values different from default
    mt_task->setName("Test of 1D coupled advection problem (dissolved Ca, Mg)");
    mt_task->setComment("@");

    // Use smart initial approximation in GEM IPM (+); SIA internal (*); AIA (-)
    mt_task->setSIA('+');
    // Set type flux Phase ( 0 undef, 1 - aq; 2 - gas; 3 - aq+gas, 4 - solids ) (default 1)
    mt_task->setTypeFluxPhase('1');
    // Use non stop debug output for nodes (+ -) (default +)
    mt_task->setOutput(true);

    // Set physical time iterator (start,end,step)
    mt_task->setTau(0, 500000, 5000);
    // Set spatial dimensions of the medium defines topology of nodes ( x y z )
    mt_task->setSpatialDimensions(0.2, 0, 0);

    // Set  M(H2O) (mass of water-solvent for molalities)
    mt_task->setMassofWaterSolvent(1.);
    // Set Maq (mass of aqueous solution for ppm etc.)
    mt_task->setMassofAqueousSolution(1.);
    // Set Vaq (volume of aqueous solution for molarities)
    mt_task->setVolumeofAqueousSolution(1.);
    // Set advection/diffusion mass transport: time step reduction factor (usually 1)
    mt_task->setTimeStepReductionFactor(1.);
    // Set initial total node volume (m^3)
    mt_task->setInitialTotalNodeVolume(0.001);
    // Set fluid advection velocity (m/sec)
    mt_task->setFluidAdvectionVelocity(1e-7);
    // Set initial node effective porosity (0 < eps < 1), usually 1
    mt_task->setInitialNodeEffectivePorosity(1.);
    // Set initial effective permeability, m2, usually 1
    mt_task->setInitialEffectivePermeability(1.);
    // Set initial value of specific longitudinal dispersivity (m), usually 1e-3
    mt_task->setInitialDispersivity(0.001);
    // Set initial general aqueous medium diffusivity (m2/sec), usually 1e-9
    mt_task->setInitialDiffusivity(0.);
    // Set initial tortuosity factor, usually 1
    mt_task->setInitialTortuosityFactor(1.);
    // Set cutoff for IC amount differences in the node between time steps (mol), usually 1e-9
    mt_task->setCutofffICamount(1e-9);
    // Set  cutoff for minimal amounts of IC in node bulk compositions (mol), usually 1e-12
    mt_task->setCutoffMinimalAmountsIC(1e-11);


    // Here we define the GEM2MT structure
    mt_task->alloc_gem2mt_arrays();

    // Update input arrays

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(TGEM2MT::pm->MassTransInit(ipm_lst, dbr_lst)) {
        return 1;  // error reading files
    }

    TGEM2MT::pm->WriteTask("gem2mt_out.dat");

    // here we call the mass-transport finite-difference coupled routine
    TGEM2MT::pm->RecCalc();

    return 0;
}

//---------------------------------------------------------------------------


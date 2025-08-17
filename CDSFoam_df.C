/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/
#include "CDS_Solver.H"   //  
#include "heRhoThermo.H"

#ifdef USE_PYTORCH
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h> //used to convert
#endif

#ifdef USE_LIBTORCH
#include <torch/script.h>
#include "DNNInferencer.H"
#endif 

//#include "fvCFD.H"
#include "dynamicFvMesh.H"
// #include "psiThermo.H"
//#include "rhoThermo.H"
//#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
//#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "PstreamGlobals.H"
#include "CombustionModel.H"

//#include "multivariateScheme.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#ifdef USE_PYTORCH
    pybind11::scoped_interpreter guard{};//start python interpreter
#endif
    #define NO_CONTROL

    #include "postProcess.H"

    // #include "setRootCaseLists.H"
    #include "listOptions.H"
    #include "setRootCase2.H"
    #include "listOutput.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"  
    
    #include "createRDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        //used for AMR
        solver.AMR_Criteria();
        mesh.update();  

        #include "centralCourantNo.H"

        if (!LTS)
        {
            #include "setDeltaT.H"
            runTime++;
        }

        solver.update_maxp();  // for smoot foil

        if (solver.ddtSchemes() == "RK2SSP" || solver.ddtSchemes() == "RK3SSP")
        {
            int rk = (solver.ddtSchemes() == "RK2SSP") ? 2 : 3;
            
            for (int nrk = 0; nrk < rk; nrk++)
            {
                Info<< nl << solver.ddtSchemes() << ": step " << nrk+1 << endl;
		solver.update();	

                if (nrk == 0)
                {
                    if (LTS)
                    {
                        #include "setRDeltaT.H"
                        runTime++;
                    }
                    Info<< "Time = " << runTime.timeName() << nl << endl;
                }
                solver.RKsolve(nrk);
            }
        }
        else
        {
            solver.update();
            if (LTS)
            {
                #include "setRDeltaT.H"
                runTime++;
            }
            Info<< "Time = " << runTime.timeName() << nl << endl;
            
            solver.Implicit_solver();
        }

        turbulence.correct();

        runTime.write();

        Info<< "============================================"<<nl<< endl;
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //

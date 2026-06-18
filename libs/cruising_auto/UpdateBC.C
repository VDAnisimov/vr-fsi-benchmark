/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "UpdateBC.H"
#include "dictionary.H"
#include "Time.H"
#include "Pstream.H"
#include "IOmanip.H"
#include "fvMesh.H"
#include "dimensionedTypes.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "wordReList.H"

#include "fixedValueFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "inletOutletFvPatchFields.H"
#include "volFields.H"
#include "PstreamReduceOps.H"
//#include "DataEntry.H" 
#include "fvsPatchFields.H"

 #include "surfaceInterpolate.H"
 #include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(UpdateBC, 0);
    addToRunTimeSelectionTable(functionObject, UpdateBC, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::UpdateBC::createFiles()
{
    if (writeToFile() && !coeffFilePtr_)
    {
        coeffFilePtr_ = createFile("coefficient");
        writeIntegratedHeader("Coefficients", coeffFilePtr_());
    }
}


void Foam::functionObjects::UpdateBC::writeIntegratedHeader
(
    const word& header,
    Ostream& os
) const
{

    writeHeader(os, "Force coefficients");
    writeCommented(os, "Time");
    writeTabbed(os, "U");
    writeTabbed(os, "F");
    writeTabbed(os, "A");
    os  << endl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::UpdateBC::UpdateBC
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool readFields
)
:
    forces(name, runTime, dict, false),
    omega_(0),
    v_st_(0.1),
    dv_start_(0.01),
    per_max_(5),
    per_first_(7),
    corr_time_(1),
    epsilon_(0.001),
    inletPatchName_(word::null),
    outletPatchName_(word::null),
    v_st1(-20),
    v_st0(0),
    F_av(0),
    F_av_last(0),
    F_av_old(0),
    timer1(1),
    per(0),
    per_m(7),
    da(0),
    coeffFilePtr_()
{
if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict, "liftDir", "dragDir");
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::UpdateBC::read(const dictionary& dict)
{
    forces::read(dict);

    rhoName_ = "rhoInf";

    dict.readEntry("inletPatchName", inletPatchName_);
    dict.readEntry("outletPatchName", outletPatchName_);
    dict.readEntry("omega", omega_);
    dict.readEntry("v_st", v_st_);

    dict.readEntry("dv_start", dv_start_);
    dict.readEntry("per_max", per_max_);
    dict.readEntry("per_first", per_first_);
    dict.readEntry("corr_time", corr_time_);
    dict.readEntry("epsilon", epsilon_);


    return true;
}


bool Foam::functionObjects::UpdateBC::execute()
{
    forces::calcForcesMoment();

    createFiles();

                     //------------------------------------------------//
	//----------// Update boundary conditions at inlet and outlet //-------------//
	           //------------------------------------------------//
        
         const fvMesh& mesh = refCast<const fvMesh>(obr_); //get mesh pointer
        
        label inletPatchId = mesh.boundaryMesh().findPatchID(inletPatchName_); //find inlet boundary by name
        label outletPatchId = mesh.boundaryMesh().findPatchID(outletPatchName_); //find outlet boundary by name
        
        const volVectorField& U = mesh.lookupObject<volVectorField>(UName_); //find velosity (U) field

	//------------------- mixed BC velocity outlet--------------------------------------
        const inletOutletFvPatchVectorField& constVv1=
    		refCast<const inletOutletFvPatchVectorField>(U.boundaryField()[outletPatchId]); // values at the boundary
        inletOutletFvPatchVectorField& vv1 = 
    		const_cast<inletOutletFvPatchVectorField&>(constVv1); //const -> no const
        vectorField& rVvv1 = vv1.refValue(); //- Value field

	//------------------- mixed BC velocity inlet--------------------------------------
        const inletOutletFvPatchVectorField& constVv2=
    		refCast<const inletOutletFvPatchVectorField>(U.boundaryField()[inletPatchId]); // values at the boundary
        inletOutletFvPatchVectorField& vv2 = 
    		const_cast<inletOutletFvPatchVectorField&>(constVv2); //const -> no const
        vectorField& rVvv2 = vv2.refValue(); //- Value field

	
        //-------------------------- average velocity--------------------------------------
        
        scalar tstep = mesh.time().deltaT().value();
        scalar frequency_ = omega_/2.0/M_PI;
        
         const auto& coordSys = coordSysPtr_();
        const Field<vector> localForce(coordSys.localVector(force_[0] + force_[1] + force_[2]));

    	//Field<vector> totForce(force_[0] + force_[1] + force_[2]); //full force (pressure + viscous parts)
        scalar Fx = localForce[0][0];
        
        if (Pstream::master()) Info << nl << Fx <<nl;

        scalar time_p = obr_.time().value()+tstep; //real time
	    time_p -= floor((obr_.time().value()+tstep)*frequency_)/frequency_; //time in a period

        if (v_st1<=-19) {
        v_st1 = v_st0 = v_st_;
        per_m = per_first_;
        per=0;
        } // start velocity 

        if (time_p<tstep)
        {
            F_av = F_av*tstep*frequency_;

             if (Pstream::master()) Info << nl << (fabs(F_av)>epsilon_) <<tab << ( per>per_m) << tab
                    << (fabs(F_av-F_av_last)<epsilon_/10) <<nl;

            //if (Pstream::master()) Info << nl << per << tab <<per_m << nl;

            if (fabs(F_av)>epsilon_ && per>per_m && fabs(F_av-F_av_last)<epsilon_/50)
            {
                timer1 = corr_time_; //correct velocity for corr_time periods (T*corr_time_)
            
                scalar dv;

                if (per_m == per_first_) //only for the first correction 
                {
                    dv = dv_start_;
                    per_m = per_max_;
                }
                else
                {
                  
                    dv = -(F_av*(v_st1-v_st0)/(F_av-F_av_old));//secant method
                }
 
                da = dv*frequency_/timer1;
                per = 0;
                F_av_old = F_av;
                v_st0 = v_st1;
            }

            if (timer1>0.99999) timer1-=1.0;
            else da = 0.0; // corrector

            F_av_last = F_av;
            F_av = 0;
            per++;
        }

        if (Pstream::master()) Info << nl << "Av. force: " << F_av_last << " Av. velocity: " <<v_st1<< nl;
        F_av = F_av + Fx;
        v_st1 = v_st1 + da*tstep;

        forAll(vv1, iFace)
    	{   

	        rVvv1[iFace] = vector(v_st1, 0, 0); //val
	    }
	    forAll(vv2, iFace)
    	{   
            rVvv2[iFace] = vector(v_st1, 0, 0); //val
	    }

    	
	    vv1.updateCoeffs(); //update boundary conditions
        vv2.updateCoeffs(); //update boundary conditions


    if (writeToFile())
    {
         writeCurrentTime(coeffFilePtr_());
        coeffFilePtr_() << tab << F_av_last << tab<< v_st1 <<endl;
    }

    return true;
}


bool Foam::functionObjects::UpdateBC::write()
{

    return true;
}


// ************************************************************************* //

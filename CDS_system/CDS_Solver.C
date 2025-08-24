
#include "CDS_Solver.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static member functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CDS_Solver, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CDS_Solver::CDS_Solver
(
    const Time& runTime,
    const fvMesh& mesh
)
:
    runTime_(runTime),
    mesh_(mesh),
    
    fluxScheme_(fluxScheme::New(mesh)),
    
    thermoPtr_(new heRhoThermo<rhoThermo, CanteraMixture>(mesh_, word::null)),
    thermo_(*thermoPtr_),
    rho_
    (
        IOobject
        (
            "rho",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.rho()
    ),
    U_
    (
        IOobject
        (
            "U",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    p_(thermo_.p()),
    T_(thermo_.T()),
    ea_(thermo_.he()),
    ha_(ea_ + p_/rho_),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*(ea_ + 0.5*magSqr(U_))
    ), 
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimVelocity*dimArea, 0.0)
    ), 
    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimDensity*dimVelocity*dimArea, 0.0)
    ),
    rhoUPhi_
    (
        IOobject
        (
            "rhoUPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimDensity*sqr(dimVelocity)*dimArea, Zero)
    ),
    rhoEPhi_
    (
        IOobject
        (
            "rhoEPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimDensity*pow3(dimVelocity)*dimArea, 0.0)
    ),

   /* turbulence_
    (
        compressible::turbulenceModel::New
        (
            rho_,
            U_,
            rhoPhi_,
            thermo_
        )
    ),

    combustion_
    (
        CombustionModel<basicThermo>::New(thermo_, turbulence_())
    ),
    chemistry_(combustion_->chemistry()),
    
    Y_(chemistry_->Y()),*/
       
    // 扩散相关场
    diffAlphaD_
    (
        IOobject
        (
            "diffAlphaD",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimTime/dimVolume, 0)
    ),
    hDiffCorrFlux_
    (
        IOobject
        (
            "hDiffCorrFlux",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimensionSet(1,0,-3,0,0,0,0), Zero)
    ),
    sumYDiffError_
    (
        IOobject
        (
            "sumYDiffError",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("sumYDiffError", dimDynamicViscosity/dimLength, Zero)
    ),   
    gradP_
    (
        IOobject
        (
            "gradP",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(fvc::grad(thermo_.p()))
    ),
    maxp_
    (
        IOobject
        (
            "maxp",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo_.p()
    ),
    normalisedGradrho
    (
        IOobject
        (
            "normalisedGradrho",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimless, 0.0)
    ),
    error  // AMR need
    (
        IOobject
        (
            "error",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        0.0,
        wordList(p_.boundaryField().size(), "zeroGradient")
    ),
    c_
    (
        IOobject
        (
            "c",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(thermo_.Cp()/(thermo_.Cv()*thermo_.psi()))
    ),
    rhoSave_
    (
        IOobject
        (
            "rhoSave",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo_.rho() 
    ),
    rhoUSave_
    (
        IOobject
        (
            "rhoUSave",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U_
    ), 
    rhoESave_
    (
        IOobject
        (
            "rhoESave",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*(ea_ + 0.5*magSqr(U_))
    )
{
    // createFields
    ddtSchemes_ = "Euler";
    if (mesh_.schemesDict().readIfPresent("timeScheme", ddtSchemes_))
    {
        if(ddtSchemes_ == "RK2SSP" || ddtSchemes_ == "RK3SSP" || ddtSchemes_ == "Euler")
        {
            Info<< "ddtSchemes: " << ddtSchemes_ << endl;
            if(ddtSchemes_ == "RK2SSP" || ddtSchemes_ == "RK3SSP")
            {
                Info<< "!! Note: RK2SSP and RK3SSP are not available for two-phase flow simulaiton. "
                    << "If you want to simulate two-phase flows, please change time scheme to 'Euler' ."
                    << endl;
            }
        }
        else
        {
            FatalErrorInFunction
                << "This timeScheme is not a valid choice. "
                << "Please use Euler, RK2SSP or RK3SSP scheme."
                << abort(FatalError);
        }
    }
    
    chemScheme_ = "ode";
    mesh_.schemesDict().readIfPresent("chemScheme", chemScheme_);
    if ((chemScheme_ == "direct") || (chemScheme_ == "ode"))
    {
        Info<< "chemScheme: " << chemScheme_ << endl;
    }
    else
    {
        FatalErrorInFunction
            << "chemScheme: " << chemScheme_
            << " is not a valid choice. "
            << "Options are: 'ode' | 'direct'"
            << abort(FatalError);
    }

    if((ddtSchemes_ != "RK2SSP") && (ddtSchemes_ != "RK3SSP"))
    {
        if(chemScheme_ == "direct")
        {
            FatalErrorInFunction
                << "This combination is not a valid choice. "
                << "If you want to use direct integrate for chemistry, please use RK2SSP or RK3SSP scheme."
                << abort(FatalError);
        }
    }
    
    Info<< "Reading thermophysical properties\n" << endl;
    CanteraMixture::setEnergyName("ea");
    
    dictionary thermoDict
    (
        IOdictionary
        (
            IOobject
            (
                "thermophysicalProperties",
                runTime_.constant(),
                mesh_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    inviscid_ = thermoDict.lookupOrDefault("inviscid",false);
    rho_ = thermo_.rho();   
    Info<< "Creating turbulence model" << endl;
    turbulence_.reset
    (
        compressible::turbulenceModel::New
        (
            rho_,
            U_,
            rhoPhi_,
            thermo_
        ).ptr()
    );
    const word turbName = mesh_.objectRegistry::lookupObject<IOdictionary>("turbulenceProperties").lookup("simulationType");
    turbName_ = turbName;

    Info<< "Creating reaction model" << endl;
    combustion_.reset
    (
        CombustionModel<basicThermo>::New(thermo_, turbulence_()).ptr()
    );
    //chemistry_ = dynamic_cast<dfChemistryModel<basicThermo>*>(combustion_->chemistry());
    chemistry_ = combustion_->chemistry();
 
   // 获取惰性组分索引
    word inertSpecie(chemistry_->lookup("inertSpecie"));
    inertIndex_ = chemistry_->species()[inertSpecie];
    chemistry_->updateEnergy();
    // 更新能量场
    //ea_ = thermo_.he();
    ha_ = ea_ + p_/rho_;
    // 更新守恒量
    rho_ = thermo_.rho();
    rhoU_ = rho_*U_;
    rhoE_ = rho_*(ea_ + 0.5*magSqr(U_));
    chemistry_->correctThermo();
    
    rhoSave_ = rho_;
    rhoUSave_ = rhoU_;
    rhoESave_ = rhoE_;

    Info<< "At initial time, min/max(T) = " 
        << min(T_).value() << ", " << max(T_).value() << endl;

    // 添加组分到场表
    PtrList<volScalarField>& Y_ = chemistry_->Y();
    forAll(Y_, i)
    {
        fields_.add(Y_[i]);
    }
    fields_.add(thermo_.he());

    // 创建组分质量分数场
    nspecies_ = chemistry_->species().size();
    rhoYi_.resize(nspecies_);
    rhoPhiYi_.resize(nspecies_);
    rhoYiSave_.resize(nspecies_);
    
    forAll(rhoYi_, i)
    {
        rhoYi_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoYi" + Y_[i].name(),
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rho_*Y_[i]
            )
        );
    }
    
    forAll(rhoPhiYi_, i)
    {
        rhoPhiYi_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "rhoPhi" + Y_[i].name(),
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("0", dimDensity*dimVelocity*dimArea, 0.0)
            )
        );
    }

    forAll(rhoYiSave_, i)
    {
        rhoYiSave_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rhoYiSave" + Y_[i].name(),
                    runTime_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                rho_*Y_[i]
            )
        );
    }    
    
    
    // 设置 Schmidt 数
    Sct_ = chemistry_->lookupOrDefault("Sct", 1.);

    // 初始化RK系数（根据时间格式）
    if (ddtSchemes_ == "RK2SSP")
    {
        rk_ = 2;
        rkcoe1_.resize(rk_);
        rkcoe2_.resize(rk_);
        rkcoe3_.resize(rk_);
        rkcoe1_[0] = 1.0; rkcoe2_[0] = 0.0; rkcoe3_[0] = 1.0;
        rkcoe1_[1] = 0.5; rkcoe2_[1] = 0.5; rkcoe3_[1] = 0.5;
    }
    else if (ddtSchemes_ == "RK3SSP")
    {
        rk_=3;
        rkcoe1_.resize(rk_);
        rkcoe2_.resize(rk_);
        rkcoe3_.resize(rk_);        
        rkcoe1_[0]=1.0; rkcoe2_[0]=0.0; rkcoe3_[0]=1.0;
        rkcoe1_[1]=0.75; rkcoe2_[1]=0.25; rkcoe3_[1]=0.25;
        rkcoe1_[2]=1.0/3.0; rkcoe2_[2]=2.0/3.0; rkcoe3_[2]=2.0/3.0;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CDS_Solver::~CDS_Solver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*void Foam::CDS_Solver::updateEnergy()
{
    // 更新能量场
    ha_ = ea_ + p_/rho_;
    chemistry_->correctThermo();
}*/

void Foam::CDS_Solver::saveFields()
{
    rhoSave_ = rho_;

    rhoUSave_ = rho_*U_;

    rhoESave_ = rho_*(ea_ + 0.5*magSqr(U_));

    PtrList<volScalarField>& Y_ = chemistry_->Y();
    rhoYiSave_.resize(nspecies_);
    forAll(rhoYiSave_,i)
    {
        rhoYiSave_[i] = rho_*Y_[i];
    }
}

void Foam::CDS_Solver::RKsolve(int nrk)
{
    volScalarField muEff("muEff", turbulence_->muEff());
    volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U_))));
    
    if (nrk == 0)
    {
        saveFields();
    }                
    
    // --- Solve density
    volScalarField rho_rhs_ = -fvc::div(rhoPhi_);
    rho_ = rkcoe1_[nrk]*rhoSave_
         + rkcoe2_[nrk]*rho_ 
         + rkcoe3_[nrk]*rho_rhs_*runTime_.deltaT();
    //Info <<"in rk"<< nrk+1 << " finish calculate rho" << nl << endl;

    // --- Solve momentum
    volVectorField rhoU_rhs_ = -fvc::div(rhoUPhi_);
    rhoU_ = rkcoe1_[nrk]*rhoUSave_
          + rkcoe2_[nrk]*rhoU_
          + rkcoe3_[nrk]*rhoU_rhs_*runTime_.deltaT();

    if (!inviscid_)
    {
        rhoU_rhs_ = fvc::laplacian(turbulence_->muEff(),U_) + fvc::div(tauMC);
        rhoU_ += rkcoe3_[nrk]*rhoU_rhs_*runTime_.deltaT();
    }

    U_ = rhoU_ / rho_;
    U_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() == rho_.boundaryField() * U_.boundaryField();             
    //Info << "in rk" << nrk+1 << ", min / max mag(U) ; " << gMin(mag(U_)()) << " / " << gMax(mag(U_)()) << nl << endl;
    Info << "\nmin / max mag(U) ; " << gMin(mag(U_)()()) << " / " << gMax(mag(U_)()()) << nl << endl;

    // --- Solve species
    PtrList<volScalarField>& Y_ = chemistry_->Y();
    if (!inviscid_)
    {
        hDiffCorrFlux_ = Zero;
        diffAlphaD_ = Zero;
        sumYDiffError_ = Zero;
        forAll(Y_, i)
        {
            sumYDiffError_ += chemistry_->rhoD(i)*fvc::grad(Y_[i]);
        }
    }
    
    if(chemScheme_ == "direct")
    {
        chemistry_->calculateW();
    }

    volScalarField Yt(0.0*Y_[0]);
    
    forAll(Y_, i)
    {
        volScalarField& Yi = Y_[i];

        volScalarField rhoYi_rhs_ = -fvc::div(rhoPhiYi_[i]);
        rhoYi_[i] = rkcoe1_[nrk]*rhoYiSave_[i]
                  + rkcoe2_[nrk]*rhoYi_[i]
                  + rkcoe3_[nrk]*rhoYi_rhs_*runTime_.deltaT();
        Yi = rhoYi_[i] / rho_;
        Yi.max(0.0);
                          
        if (!inviscid_)
        {      
            hDiffCorrFlux_ += chemistry_->hei(i)*(chemistry_->rhoD(i)*fvc::grad(Yi) - Yi*sumYDiffError_);
            diffAlphaD_ += fvc::laplacian(thermo_.alpha()*chemistry_->hei(i), Yi);
            const surfaceScalarField phiUc = linearInterpolate(sumYDiffError_) & mesh_.Sf();
            tmp<volScalarField> DEff = chemistry_->rhoD(i) + turbulence_->mut()/Sct_;
            rhoYi_rhs_ =(
                   turbName_ == "laminar"?(fvc::laplacian(DEff(), Yi) - fvc::div(phiUc,Yi,"div(phi,Yi_h)")):fvc::laplacian(DEff(), Yi));

            rhoYi_[i] += rkcoe3_[nrk]*rhoYi_rhs_*runTime_.deltaT();
            Yi = rhoYi_[i]/rho_;
            Yi.max(0.0);
        }
        Yi.correctBoundaryConditions();
        rhoYi_[i] = rho_*Yi;
        Yt += Yi;
    }
    
    forAll(Y_, i)
    {
        Y_[i] = Y_[i]/Yt;
        Y_[i].max(0.0);
        rhoYi_[i] = rho_*Y_[i];
    }
 
    // --- Solve energy 
    surfaceScalarField sigmaDotU
    (
        "sigmaDotU",
        (
            fvc::interpolate(turbulence_->muEff()) * mesh_.magSf() * fvc::snGrad(U_)
          + (mesh_.Sf() & fvc::interpolate(tauMC))
        )
        & (a_pos() * U_pos() + a_neg() * U_neg())
    );

    volScalarField rhoE_rhs_ = -fvc::div(rhoEPhi_) + fvc::div(sigmaDotU);
    rhoE_ = rkcoe1_[nrk]*rhoESave_
          + rkcoe2_[nrk]*rhoE_
          + rkcoe3_[nrk]*rhoE_rhs_*runTime_.deltaT();

    ea_ = rhoE_ / rho_ - 0.5 * magSqr(U_);
    ea_.correctBoundaryConditions();

    chemistry_->correctThermo();
    rhoE_.boundaryFieldRef() == rho_.boundaryField() * (ea_.boundaryField() + 0.5 * magSqr(U_.boundaryField()));

    if (!inviscid_)
    {
        rhoE_rhs_ = fvc::laplacian(turbulence_->alphaEff() * thermo_.gamma(), ea_)
                 - diffAlphaD_
                 + fvc::div(hDiffCorrFlux_);

        rhoE_ += rkcoe3_[nrk] * rhoE_rhs_ * runTime_.deltaT();
        ea_ = rhoE_ / rho_ - 0.5 * magSqr(U_);
        ea_.correctBoundaryConditions();

        chemistry_->correctThermo();
    }

    Info << "min/max(T) = " << min(T_).value() << ", " << max(T_).value() << endl;

    p_.ref() = rho_() / thermo_.psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() == thermo_.psi().boundaryField() * p_.boundaryField();
    
    if ((nrk == rk_-1) && (chemScheme_ == "ode"))
    {
        calculateR();
    }
}

void Foam::CDS_Solver::calculateR()
{
    combustion_->correct();
    PtrList<volScalarField>& Y_ = chemistry_->Y();
    volScalarField Yt(0.0*Y_[0]);
    forAll(Y_, i)
    {
        volScalarField& Yi = Y_[i];
                 
        rhoYi_[i].ref() += chemistry_->RR(i) * runTime_.deltaT();
        Info << "max reaction rate " << Yi.name() 
                 << " is " << max(chemistry_->RR(i)).value() << endl;

        Yi = rhoYi_[i] / rho_;
        Yi.max(0.0);
        Yi.correctBoundaryConditions();
        rhoYi_[i] = rho_ * Yi;
        Yt += Yi;
    }
    
    forAll(Y_, i)
    {
        Y_[i] = Y_[i]/Yt;
        Y_[i].max(0.0);
        rhoYi_[i] = rho_*Y_[i];
    }

    chemistry_->correctThermo();
    Info << "min/max(T) = "
         << min(T_).value() << ", " << max(T_).value() << endl;
    
    p_.ref() = rho_() / thermo_.psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() == thermo_.psi().boundaryField() * p_.boundaryField();
    
    Info << "Finish calculate Reaction" << nl << endl;
}

void Foam::CDS_Solver::fluxPredictor()
{
    if (!pos.valid())
    {
        pos = surfaceScalarField::New
        (
            "pos",
            mesh_,
            dimensionedScalar(dimless, 1)
        );

        neg = surfaceScalarField::New
        (
            "neg",
            mesh_,
            dimensionedScalar(dimless, -1.0)
        );
    }

    const volScalarField p = p_;
    const volScalarField rho = rho_;
    const volVectorField U = U_;

    rho_pos = interpolate(rho, pos());
    rho_neg = interpolate(rho, neg());

    const volVectorField rhoU(rho*U);
    rhoU_pos = interpolate(rhoU, pos(), U.name());
    rhoU_neg = interpolate(rhoU, neg(), U.name());

    U_pos = surfaceVectorField::New("U_pos", rhoU_pos()/rho_pos());
    U_neg = surfaceVectorField::New("U_neg", rhoU_neg()/rho_neg());

    const volScalarField& T = T_;

    const volScalarField rPsi("rPsi", 1.0/thermo_.psi());
    const surfaceScalarField rPsi_pos(interpolate(rPsi, pos(), T.name()));
    const surfaceScalarField rPsi_neg(interpolate(rPsi, neg(), T.name()));

    p_pos = surfaceScalarField::New("p_pos", rho_pos()*rPsi_pos);
    p_neg = surfaceScalarField::New("p_neg", rho_neg()*rPsi_neg);

    surfaceScalarField phiv_pos("phiv_pos", U_pos() & mesh_.Sf());
    surfaceScalarField phiv_neg("phiv_neg", U_neg() & mesh_.Sf());

    // Make fluxes relative to mesh-motion
    if (mesh_.moving())
    {
        phiv_pos -= mesh_.phi();
        phiv_neg -= mesh_.phi();
    }

    const volScalarField c("c", sqrt(thermo_.Cp()/thermo_.Cv()*rPsi));
    const surfaceScalarField cSf_pos
    (
        "cSf_pos",
        interpolate(c, pos(), T.name())*mesh_.magSf()
    );
    const surfaceScalarField cSf_neg
    (
        "cSf_neg",
        interpolate(c, neg(), T.name())*mesh_.magSf()
    );

    const dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0);

    const surfaceScalarField ap
    (
        "ap",
        max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
    );
    const surfaceScalarField am
    (
        "am",
        min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
    );

    word fluxScheme("Kurganov");
    mesh_.schemesDict().readIfPresent("fluxScheme", fluxScheme);

    a_pos = surfaceScalarField::New
    (
        "a_pos",
        fluxScheme == "Tadmor"
          ? surfaceScalarField::New("a_pos", mesh_, 0.5)
          : ap/(ap - am)
    );

    a_neg = surfaceScalarField::New("a_neg", 1.0 - a_pos());

    phiv_pos *= a_pos();
    phiv_neg *= a_neg();

    aSf = surfaceScalarField::New
    (
        "aSf",
        fluxScheme == "Tadmor"
          ? -0.5*max(mag(am), mag(ap))
          : am*a_pos()
    );

    aphiv_pos = surfaceScalarField::New("aphiv_pos", phiv_pos - aSf());
    aphiv_neg = surfaceScalarField::New("aphiv_neg", phiv_neg + aSf());
}

void Foam::CDS_Solver::Implicit_solver()
{

    volScalarField muEff("muEff", turbulence_->muEff());
    volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U_))));

    solve
    (
        fvm::ddt(rho_)
      + fvc::div(rhoPhi_)
    );

    solve
    (
          fvm::ddt(rhoU_)
        + fvc::div(rhoUPhi_)
    );
    U_.ref() = rhoU_()/rho_();
    U_.correctBoundaryConditions();
    rhoU_.boundaryFieldRef() == rho_.boundaryField()*U_.boundaryField();

    if (!inviscid_)
    {
        solve
        (
              fvm::ddt(rho_, U_) - fvc::ddt(rho_, U_)
           -  fvm::laplacian(turbulence_->muEff(), U_)
           -  fvc::div(tauMC)
        );
        rhoU_ = rho_*U_;
    }

    Info << "\nmin / max mag(U) ; " << gMin(mag(U_)()()) << " / " << gMax(mag(U_)()()) << nl << endl;
    
PtrList<volScalarField>& Y_ = chemistry_->Y();
    if (!inviscid_)
    {
        hDiffCorrFlux_ = Zero;
        diffAlphaD_ = Zero;
        sumYDiffError_ = Zero;
        forAll(Y_, i)
        {
            sumYDiffError_ += chemistry_->rhoD(i)*fvc::grad(Y_[i]);
        }
    }

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh_,
            fields_,
            rhoPhi_,
            mesh_.divScheme("div(phi,Yi_h)")
        )
    );
    
    combustion_->correct();
    
    volScalarField Yt(0.0*Y_[0]);

    forAll(Y_, i)
    {
        volScalarField& Yi = Y_[i];

        solve
        (
            fvm::ddt(rhoYi_[i])
          + fvc::div(rhoPhiYi_[i])
         ==
            chemistry_->RR(i)
        );
        Yi = rhoYi_[i] / rho_;
        Yi.max(0.0);
                          
        if (!inviscid_)
        {      
            hDiffCorrFlux_ += chemistry_->hei(i)*(chemistry_->rhoD(i)*fvc::grad(Yi) - Yi*sumYDiffError_);
            diffAlphaD_ += fvc::laplacian(thermo_.alpha()*chemistry_->hei(i), Yi);
            const surfaceScalarField phiUc = linearInterpolate(sumYDiffError_) & mesh_.Sf();
            tmp<volScalarField> DEff = chemistry_->rhoD(i) + turbulence_->mut()/Sct_;

            fvScalarMatrix YiEqn
            (
                    fvm::ddt(rho_, Yi) - fvc::ddt(rho_, Yi)
                 -
                    (
                        turbName_ == "laminar"
                        ?  (fvm::laplacian(DEff(), Yi) - mvConvection->fvmDiv(phiUc, Yi))
                        :  fvm::laplacian(DEff(), Yi)
                    )
            );

            YiEqn.relax();
            YiEqn.solve("Yi");
            Yi.max(0.0);
        }
        Yi.correctBoundaryConditions();
        rhoYi_[i] = rho_*Yi;
        Yt += Yi;
    }
    
    forAll(Y_, i)
    {
        Y_[i] = Y_[i]/Yt;
        Y_[i].max(0.0);
        rhoYi_[i] = rho_*Y_[i];
    }
 
    // --- Solve energy
    
    surfaceScalarField sigmaDotU
    (
        "sigmaDotU",
        (
            fvc::interpolate(turbulence_->muEff()) * mesh_.magSf() * fvc::snGrad(U_)
          + (mesh_.Sf() & fvc::interpolate(tauMC))
        )
        & (a_pos() * U_pos() + a_neg() * U_neg())
    );

    solve
    (
            fvm::ddt(rhoE_)
          + fvc::div(rhoEPhi_)
          - fvc::div(sigmaDotU)
    );

    ea_ = rhoE_ / rho_ - 0.5 * magSqr(U_);
    ea_.correctBoundaryConditions();

    chemistry_->correctThermo();
    rhoE_.boundaryFieldRef() == rho_.boundaryField() * (ea_.boundaryField() + 0.5 * magSqr(U_.boundaryField()));

    if (!inviscid_)
    {
        fvScalarMatrix eEqn
        (
              fvm::ddt(rho_, ea_) - fvc::ddt(rho_, ea_)
            ==
                (
                    turbName_ == "laminar"
                    ?
                    (
                        // alpha in deepflame is considered to calculate h by default (kappa/Cp), so multiply gamma to correct alpha
                        fvm::laplacian(turbulence_->alphaEff()*thermo_.gamma(), ea_)
                    -   diffAlphaD_
                    +   fvc::div(hDiffCorrFlux_)
                    )
                    :
                    (
                        fvm::laplacian(turbulence_->alphaEff()*thermo_.gamma(), ea_)
                    )
                )
        );

        eEqn.solve("ea_");
        chemistry_->correctThermo();
        rhoE_ = rho_*(ea_ + 0.5*magSqr(U_));
        rhoE_.boundaryFieldRef() == rho_.boundaryField()*(ea_.boundaryField() + 0.5*magSqr(U_.boundaryField()));
    }

    Info << "min/max(T) = " << min(T_).value() << ", " << max(T_).value() << endl;

    p_.ref() = rho_() / thermo_.psi();
    p_.correctBoundaryConditions();
    rho_.boundaryFieldRef() == thermo_.psi().boundaryField() * p_.boundaryField();
}

void Foam::CDS_Solver::update()
{
    fluxScheme_->update
    (
        rho_,
        rhoYi_,
        nspecies_,
        U_,
        ea_,
        p_,
        c_,
        phi_,
        rhoPhi_,
        rhoPhiYi_,
        rhoUPhi_,
        rhoEPhi_
    );
}

void Foam::CDS_Solver::update_maxp()
{
    gradP_ = mag(fvc::grad(p_));
    c_ = sqrt(thermo_.Cp()/(thermo_.Cv()*thermo_.psi()));
    forAll(maxp_, celli)
    {
        if (p_[celli] > maxp_[celli])
        {
            maxp_[celli] = p_[celli];
        }
    }
    gradP_.correctBoundaryConditions();
    c_.correctBoundaryConditions();
    maxp_.correctBoundaryConditions(); 
}

void Foam::CDS_Solver::AMR_Criteria()
{
    tmp<volScalarField> tmagGradrho = mag(fvc::grad(rho_));
    dimensionedScalar refCri
    (
        "refCri",
        dimensionSet(1, -4, 0, 0, 0),
        max(tmagGradrho()).value()
    );
    normalisedGradrho = tmagGradrho() / refCri;
    tmagGradrho.clear();

    // Calculate error across faces based on density gradient
    if (!isA<staticFvMesh>(mesh_))
    {
        vectorField gradRho(fvc::grad(rho_)().primitiveField());
        scalarField dL(mesh_.V()/fvc::surfaceSum(mesh_.magSf()));

        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();
        const label nInternalFaces = mesh_.nInternalFaces();
        error = 0.0;

        vector solutionD((vector(mesh_.solutionD()) + vector::one)/2.0);

        for (label facei = 0; facei < nInternalFaces; facei++)
        {
            label own = owner[facei];
            label nei = neighbour[facei];
            vector dr = mesh_.C()[nei] - mesh_.C()[own];
            scalar magdr = mag(dr);

            // Ignore error in empty directions
            if (mag(solutionD & (dr/magdr)) > 0.1)
            {
                scalar dRhodr = (rho_[nei] - rho_[own])/magdr;
                scalar rhoc = (rho_[nei] + rho_[own])*0.5;
                scalar dl = (dL[own] + dL[nei])*0.5;
                scalar dRhoDotOwn = gradRho[own] & (dr/magdr);
                scalar dRhoDotNei = gradRho[nei] & (-dr/magdr);
                scalar eT =
                    max
                    (
                        mag(dRhodr - dRhoDotNei)/(0.3*rhoc/dl + mag(dRhoDotNei)),
                        mag(dRhodr - dRhoDotOwn)/(0.3*rhoc/dl + mag(dRhoDotOwn))
                    );
                error[own] = max(error[own], eT);
                error[nei] = max(error[nei], eT);
            }
        }
        error.correctBoundaryConditions();
    }

}

Foam::tmp<Foam::surfaceScalarField>
Foam::CDS_Solver::amaxSf() const
{
    return max(mag(aphiv_pos()), mag(aphiv_neg()));
}

// ************************************************************************* //

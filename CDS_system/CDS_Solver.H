
#ifndef CDS_Solver_H
#define CDS_Solver_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dfChemistryModel.H"
#include "CanteraMixture.H"

#include "fvCFD.H"
#include "fluxScheme.H"
#include "rhoThermo.H"
#include "heRhoThermo.H"
#include "CombustionModel.H"
#include "turbulentFluidThermoModel.H"
#include "multivariateScheme.H"

#include "staticFvMesh.H" //


namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class CDS_Solver Declaration
\*---------------------------------------------------------------------------*/

class CDS_Solver
{
protected:
// Protected data
    const Time& runTime_;
    const fvMesh& mesh_;
    
    autoPtr<fluxScheme> fluxScheme_;

    // 时间积分方案
    word ddtSchemes_;
    
    // 化学求解方案
    word chemScheme_;
    
    // 热物理模型
    autoPtr<rhoThermo> thermoPtr_;
    rhoThermo& thermo_;

    // 基本场
    volScalarField rho_;
    volVectorField U_;
    volScalarField& p_;
    const volScalarField& T_;
    volScalarField& ea_;
    volScalarField ha_;
    
    // 守恒量场
    volVectorField rhoU_;
    volScalarField rhoE_;
    
    // 通量场
    surfaceScalarField phi_;
    surfaceScalarField rhoPhi_;
    surfaceVectorField rhoUPhi_;
    surfaceScalarField rhoEPhi_;
    
    // 模型
    //autoPtr<compressible::turbulenceModel> turbulence_;
    //autoPtr<CombustionModel<basicThermo>> combustion_;
    //autoPtr<dfChemistryModel<basicThermo>> chemistry_;
    autoPtr<compressible::turbulenceModel> turbulence_;
    autoPtr<CombustionModel<basicThermo>> combustion_;    
    dfChemistryModel<basicThermo>* chemistry_;
    
    
    //PtrList<volScalarField>& Y_;

    // 多变量插值表
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;
    
    // 参数
    bool inviscid_;
    word turbName_;
    scalar Sct_;
    label inertIndex_;
    label nspecies_;
    
    // 扩散相关场
    volScalarField diffAlphaD_;
    volVectorField hDiffCorrFlux_;
    volVectorField sumYDiffError_;
    
    // Post-process
    volScalarField gradP_;
    volScalarField maxp_;
    volScalarField normalisedGradrho;
    volScalarField error;
    volScalarField c_;
    
    // 组分场
    PtrList<volScalarField> rhoYi_;
    PtrList<surfaceScalarField> rhoPhiYi_;
    

    // 新增成员变量：RK时间格式相关
    scalar rk_;                       // RK步数
    List<scalar> rkcoe1_;             // RK系数1
    List<scalar> rkcoe2_;             // RK系数2
    List<scalar> rkcoe3_;             // RK系数3

    volScalarField rhoSave_;
    volVectorField rhoUSave_;
    volScalarField rhoESave_;
    PtrList<volScalarField> rhoYiSave_;


    // Cached temporary fields

        tmp<surfaceScalarField> pos;
        tmp<surfaceScalarField> neg;

        tmp<surfaceScalarField> rho_pos;
        tmp<surfaceScalarField> rho_neg;

        tmp<surfaceVectorField> rhoU_pos;
        tmp<surfaceVectorField> rhoU_neg;

        tmp<surfaceVectorField> U_pos;
        tmp<surfaceVectorField> U_neg;

        tmp<surfaceScalarField> p_pos;
        tmp<surfaceScalarField> p_neg;

        tmp<surfaceScalarField> a_pos;
        tmp<surfaceScalarField> a_neg;

        tmp<surfaceScalarField> aSf;

        tmp<surfaceScalarField> aphiv_pos;
        tmp<surfaceScalarField> aphiv_neg;

private:

    // Private Member Functions
        
        //- Interpolate field vf according to direction dir
        template<class Type>
        inline tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const surfaceScalarField& dir,
            const word& reconFieldName = word::null
        )
        {
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf
            (
                fvc::interpolate
                (
                    vf,
            	    dir,
            	    "reconstruct("
          	    + (reconFieldName != word::null ? reconFieldName : vf.name())
          	    + ')'
        	)
    	    );

    	    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    	    sf.rename(vf.name() + '_' + dir.name());

    	    return tsf;
        }


public:

    TypeName("CDS_Solver");

    // Constructor
    CDS_Solver(const Time& runTime, const fvMesh& mesh);

    //- Destructor
    virtual ~CDS_Solver();
    


    // Member functions

        //- Save old value
        virtual void saveFields();
        
        //- For CFL, from rhocentralFoam
        virtual void fluxPredictor();

        //- Runge-Kutta solver
        virtual void RKsolve(int nrk=0);
        
        //- Implicit solver
        virtual void Implicit_solver();

        //- Calculate reaction
        virtual void calculateR();

        //- Update fluxes
        virtual void update();
        
        //- Detonation cell
        virtual void update_maxp();
        
        //- Detonation cell
        virtual void AMR_Criteria();

    // Member Access Functions
    // 访问函数
    //const fvMesh& mesh() const { return mesh_; }
    //const Time& runTime() const { return runTime_; }

    // 获取基本场引用
    const volScalarField& rho() const { return rho_; }
    const volVectorField& U() const { return U_; }
    const volScalarField& p() const { return p_; }
    const volScalarField& T() const { return T_; }
    
    // 获取通量引用
    const surfaceScalarField& phi() const { return phi_; }
    const surfaceScalarField& rhoPhi() const { return rhoPhi_; }
    const surfaceVectorField& rhoUPhi() const { return rhoUPhi_; }
    const surfaceScalarField& rhoEPhi() const { return rhoEPhi_; }
    
    // 获取模型引用
    rhoThermo& thermo() { return thermo_; }
    compressible::turbulenceModel& turbulence() { return turbulence_(); }
    CombustionModel<basicThermo>& combustion() { return combustion_(); }
    dfChemistryModel<basicThermo>& chemistry() { return *chemistry_; }
    
    // 其他场引用
    //volScalarField& gradP() { return gradP_; }
    //volScalarField& maxp() { return maxp_; }   
         
    //- Return speed of sound
    tmp<surfaceScalarField> amaxSf() const;
    tmp<volScalarField> speedOfSound() const;
    
    // 添加访问函数，返回ddtSchemes_的值
    const word& ddtSchemes() const
    {
        return ddtSchemes_;
    }
    
    /*const word& chemScheme() const
    {
        return chemScheme_;
    }*/
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


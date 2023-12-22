/**
 *  @file ChungTransport.h
 *  Interface for class ChungTransport
 */

#ifndef CT_TRANSCRITICALTRAN_H
#define CT_TRANSCRITICALTRAN_H

// // Cantera includes
// #include "MultiTransport.h"
// Cantera includes
#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/transport/MultiTransport.h"


namespace Cantera
{

//! Class MultiTransport implements transport properties for
//! high pressure gas mixtures.
/*!
 * Chung's methods are used for viscosity and thermal conductivity, due to the
 * fact that i) we don't have critical properties for radicals, ii) the
 * computational cost using the mixing rule withing Chung's method is expansive,
 * iii) we have singularities in viscosity when we use the mixing rule by Chung,
 * we apply a simple mole-fraction averaged mixing rule for both viscosity and
 * thermal conductivity;
 *
 * Takahashi's correction is applied to binary mixture coefficients, before the
 * evaluation of diffusion coeffs through function getMultiDiffCoeffs() and
 * getThermalDiffCoeffs().
 *
 * This is a derived class of MultiTransport, although for we are not
 * calculating thermal conductivity as in class MuultiTransport.
 *
 * @ingroup tranprops
 */
class ChungTransport : public MultiTransport
{
public:
    //! default constructor
    /*!
     *   @param thermo  Optional parameter for the pointer to the ThermoPhase object
     */
    //ChungTransport(thermo_t* thermo = 0){};
    ChungTransport(ThermoPhase* thermo=0);

    virtual std::string transportType() const { 
        return "Chung"; 
    }

    // Chung's high pressure viscosity, mole-fraction averaged
    virtual doublereal viscosity();

    // Chung's high pressure thermal conductivity, mole-fraction averaged
    // virtual doublereal thermalConductivity();
    virtual double thermalConductivity();
    // this is where we modify binary diff coeffs with Takahashi
    // and this will be used when we evaluate multidiff and thermaldiff
    virtual void updateThermal_T();

    virtual void init(ThermoPhase* thermo, int mode = 0, int log_level = 0);
    virtual void getMixDiffCoeffsMass(doublereal* const d);

    friend class TransportFactory;
protected:
    // hard-coded, same as in class PengRobinsonGasPhase
    // TODO: maybe better share this with class PengRobinsonGasPhase
    virtual void ReadCriticalProperties();

    // this is for Takahashi copied from class HighPressureGasTransport
    // Set value of parameter values for Takahashi correlation, by interpolating
    // table of constants vs. Pr
    virtual doublereal setPcorr(doublereal Pr, doublereal Tr);

    virtual doublereal Tcrit_i(size_t i);

    virtual doublereal Pcrit_i(size_t i);

    virtual doublereal Vcrit_i(size_t i);

    virtual doublereal Zcrit_i(size_t i);


    vector_fp store(size_t i, size_t nsp);

    // critical properties
    vector_int IsCrit;
    vector_fp Tcrit;
    vector_fp Pcrit;
    vector_fp Vcrit;
    vector_fp Zcrit;
    vector_fp omega;
    vector_fp dipole;
    vector_fp kappa;
};

}

#endif

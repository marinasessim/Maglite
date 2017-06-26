/****************************************************************/
/*            Elastic Recoil Cross Section                      */
/****************************************************************/

#ifndef ELASTICRECOILCROSSSECTIONUSEROBJECT_H
#define ELASTICRECOILCROSSSECTIONUSEROBJECT_H

// MOOSE includes
#include "GeneralUserObject.h"

// Forward Declarations
class ElasticRecoilCrossSectionUserObject;

template <>
InputParameters validParams<ElasticRecoilCrossSectionUserObject>();

class ElasticRecoilCrossSectionUserObject : public GeneralUserObject
{
public:
  ElasticRecoilCrossSectionUserObject(const InputParameters & parameters);

  virtual ~ElasticRecoilCrossSectionUserObject();

  virtual void initialize() override;

  virtual void execute() override;

  virtual void finalize() override;

protected:
  /// Function representing the neutron spectrum
  Function & _neutron_spectrum;

  /// Function representing the scattering law
  Function & _scattering_law;

  /// Function representing the elastic cross section
  Function & _elastic_xs;

  /// Method that implements Legendre polynomails up to n = 5
  Real legendreP(unsigned int n, Real x);

  /// Method that finds the neutron energy group given a neutron energy
  unsigned int findNeutronEnergyGroup(Real energy);

  QBase * _quadrature;

  ///@{ Extract 1 dim quadrature rule from libmesh and store in points and weights
  std::vector<Real> _quad_points;
  std::vector<Real> _quad_weights;
  ///@}

  /// Order of Legendre polynomials, defaults to 5
  unsigned int _L;

  /// Number of recoil energy bins
  unsigned int _T;

  /// Number of neutron energy groups
  unsigned int _G;

  /// Atomic mass of the isotope being analyzed
  const Real _atomic_mass;

  /// Gamma is a mass fraction that limits the energy transfer from the neutron to the recoil
  Real _gamma;

  /// Limits that characterize the neutron energy groups, must be in descending order
  const std::vector<Real> _neutron_energy_limits;

  /// Limits that characterize the recoil energy bins, must be in descending order
  const std::vector<Real> _recoil_energy_limits;

  /// Neutron spectrum in our environment, usually 1/E
  std::vector<Real> _xi_g;

  /// alpha = 4 * A / (A + 1)^2
  Real _alpha;

  /**
   * Elastic recoil cross section (g -> t)
   * L (order of Legendre polynomials expansion)
   * T (recoil energy bins)
   * G (neutron energy groups)
   */
  std::vector<std::vector<std::vector<Real>>> _erxs_coeff;
  std::vector<std::vector<std::vector<Real>>> _erxs_sum;

  /// Name of the output file
  std::string _file_name;
};

#endif

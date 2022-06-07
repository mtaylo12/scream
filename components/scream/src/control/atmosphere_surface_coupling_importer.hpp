#ifndef SCREAM_IMPORTER_HPP
#define SCREAM_IMPORTER_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/atm_process/SCDataManager.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class SurfaceCouplingImporter : public AtmosphereProcess
{
public:

  template<typename DevT, typename DataT>
  using view_1d = typename KokkosTypes<DevT>::template view_1d<DataT>;
  template<typename DevT, typename DataT>
  using view_2d = typename KokkosTypes<DevT>::template view_2d<DataT>;

  template<typename DevT, typename ScalarT>
  using uview_1d = Unmanaged<view_1d<DevT, ScalarT>>;
  template<typename DevT, typename ScalarT>
  using uview_2d = Unmanaged<view_2d<DevT, ScalarT>>;

  // Constructors
  SurfaceCouplingImporter (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const {
    return AtmosphereProcessType::SurfaceCouplingImporter;
  }

  // The name of the subcomponent
  std::string name () const { return "SurfaceCouplingImporter"; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const int dt);
  void finalize_impl   ();

  void setup_surface_coupling_data(const SCDataManager &sc_data_manager);

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;

  // Number of imports and exports
  Int m_num_imports;

  // Views storing a 2d array with dims (num_cols,num_fields) for import/export data.
  // The field idx strides faster, since that's what mct does (so we can "view" the
  // pointer to the whole x2a and a2x arrays from Fortran)
  view_2d <DefaultDevice, Real> m_cpl_imports_view_d;
  uview_2d<HostDevice,    Real> m_cpl_imports_view_h;


}; // class SurfaceCouplingImporter

} // namespace scream

#endif // SCREAM_CLD_FRACTION_HPP

#ifndef SCREAM_EXPORTER_HPP
#define SCREAM_EXPORTER_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
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

class SurfaceCouplingExporter : public AtmosphereProcess
{
public:

  using PF      = scream::PhysicsFunctions<DefaultDevice>;
  using KT      = ekat::KokkosTypes<DefaultDevice>;
  using Spack   = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;

  template<typename DevT, typename DataT>
  using view_1d = typename KokkosTypes<DevT>::template view_1d<DataT>;
  template<typename DevT, typename DataT>
  using view_2d = typename KokkosTypes<DevT>::template view_2d<DataT>;

  template<typename DevT, typename ScalarT>
  using uview_1d = Unmanaged<view_1d<DevT, ScalarT>>;
  template<typename DevT, typename ScalarT>
  using uview_2d = Unmanaged<view_2d<DevT, ScalarT>>;

  // Constructors
  SurfaceCouplingExporter (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const {
    return AtmosphereProcessType::SurfaceCouplingExporter;
  }

  // The name of the subcomponent
  std::string name () const { return "SurfaceCouplingExporter"; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_1d_scalar     = 7;
    static constexpr int num_2d_vector_mid = 2;
    static constexpr int num_2d_vector_int = 1;

    uview_1d<DefaultDevice, Real> Sa_z;
    uview_1d<DefaultDevice, Real> Sa_ptem;
    uview_1d<DefaultDevice, Real> Sa_dens;
    uview_1d<DefaultDevice, Real> Sa_pslv;
    uview_1d<DefaultDevice, Real> Faxa_rainl;
    uview_1d<DefaultDevice, Real> Faxa_snowl;
    uview_1d<DefaultDevice, Real> zero_view;

    uview_2d<DefaultDevice, Spack> dz;
    uview_2d<DefaultDevice, Spack> z_mid;
    uview_2d<DefaultDevice, Spack> z_int;
  };

  // A device-friendly helper struct, storing column information about the import/export.
  struct ColumnInfo {
    // Set to invalid, for ease of checking
    KOKKOS_INLINE_FUNCTION
    ColumnInfo () : data(nullptr) {}

    KOKKOS_INLINE_FUNCTION
    ColumnInfo& operator= (const ColumnInfo&) = default;

    // Stride between the 1st entry of two consecutive columns to be imported/exported.
    // Note: this is >= that number of scalars in a column. E.g., for a vector field layout like
    //       (ncols,2,nlevs), where we import/export only the 1st vector component, the stride
    //       is 2*nlevs
    int col_stride;

    // Offset to surface field from the column start. Should be 0 for scalar fields, but
    // may be non-zero for vector quantities for which we import/export the 2nd (or larger)
    // component (the layout would be something like (num_cols,2,num_levs), so the 1st
    // entry to import export would be at index num_levs.
    int col_offset;

    // For Computed fields, we do not want to export during initialization as these field
    // are not guarenteed to have valid entries.
    bool export_during_initialization;

    // Pointer to the scream field device memory
    Real* data;
  };

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const int dt);
  void finalize_impl   ();

  // Creates an helper field, not to be shared with the AD's FieldManager
  void create_helper_field (const std::string& name,
                            const FieldLayout& layout,
                            const std::string& grid);

  // Query if a local field exists
  bool has_helper_field (const std::string& name) const { return m_helper_fields.find(name)!=m_helper_fields.end(); }

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  void setup_surface_coupling_data(const SCDataManager &sc_data_manager);

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;

  // Some helper fields.
  std::map<std::string,Field> m_helper_fields;

  // Struct which contains local variables
  Buffer m_buffer;

  // Number of imports and exports
  Int m_num_exports;

  // Views storing a 2d array with dims (num_cols,num_fields) for cpl export data.
  // The field idx strides faster, since that's what mct does (so we can "view" the
  // pointer to the whole a2x array from Fortran)
  view_2d <DefaultDevice, Real> m_cpl_exports_view_d;
  uview_2d<HostDevice,    Real> m_cpl_exports_view_h;

  // Column info used during export
  view_1d<DefaultDevice, ColumnInfo> m_info;

}; // class SurfaceCouplingExporter

} // namespace scream

#endif // SCREAM_CLD_FRACTION_HPP

#include "atmosphere_surface_coupling_exporter.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
// =========================================================================================
SurfaceCouplingExporter::SurfaceCouplingExporter (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{

}
// =========================================================================================
void SurfaceCouplingExporter::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  const auto& grid_name = m_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);

  m_num_cols = grid->get_num_local_dofs();      // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels(); // Number of levels per column

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Qunit = kg/kg;
  Qunit.set_string("kg/kg");
  auto Wm2 = W / m / m;
  Wm2.set_string("W/m2)");

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  FieldLayout scalar2d_layout      { {COL   },      {m_num_cols                 } };
  FieldLayout horiz_wind_layout    { {COL,CMP,LEV}, {m_num_cols, 2, m_num_levs  } };
  FieldLayout scalar3d_layout_mid  { {COL,LEV},     {m_num_cols,    m_num_levs  } };
  FieldLayout scalar3d_layout_int  { {COL,ILEV},    {m_num_cols,    m_num_levs+1} };

  add_field<Required>("p_mid",            scalar3d_layout_mid,  Pa,    grid_name);
  add_field<Required>("p_int",            scalar3d_layout_int,  Pa,    grid_name);
  add_field<Required>("pseudo_density",   scalar3d_layout_mid,  Pa,    grid_name);
  add_field<Required>("horiz_winds",     horiz_wind_layout,    m/s,   grid_name);
  add_field<Required>("T_mid",           scalar3d_layout_mid,  K,     grid_name);
  add_field<Required>("qv",              scalar3d_layout_mid,  Qunit, grid_name, "tracers");
  add_field<Computed>("precip_liq_surf",  scalar2d_layout,      m/s,   grid_name);
  add_field<Computed>("sfc_flux_dir_nir", scalar2d_layout,      Wm2,   grid_name);
  add_field<Computed>("sfc_flux_dir_vis", scalar2d_layout,      Wm2,   grid_name);
  add_field<Computed>("sfc_flux_dif_nir", scalar2d_layout,      Wm2,   grid_name);
  add_field<Computed>("sfc_flux_dif_vis", scalar2d_layout,      Wm2,   grid_name);
  add_field<Computed>("sfc_flux_sw_net" , scalar2d_layout,      Wm2,   grid_name);
  add_field<Computed>("sfc_flux_lw_dn"  , scalar2d_layout,      Wm2,   grid_name);

  create_helper_field("Sa_z",       scalar2d_layout, grid_name);
  create_helper_field("Sa_ptem",    scalar2d_layout, grid_name);
  create_helper_field("Sa_dens",    scalar2d_layout, grid_name);
  create_helper_field("Sa_pslv",    scalar2d_layout, grid_name);
  create_helper_field("Faxa_rainl", scalar2d_layout, grid_name);
  create_helper_field("Faxa_snowl", scalar2d_layout, grid_name);
  create_helper_field("zero_view",  scalar2d_layout, grid_name);
}
// =========================================================================================
void SurfaceCouplingExporter::create_helper_field (const std::string& name,
                                                   const FieldLayout& layout,
                                                   const std::string& grid)
{
  using namespace ekat::units;
  FieldIdentifier id(name,layout,Units::nondimensional(),grid);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation();
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}
// =========================================================================================
size_t SurfaceCouplingExporter::requested_buffer_size_in_bytes() const
{
  // Number of Reals needed by local views in the interface
  return /*Buffer::num_1d_scalar*m_num_cols*sizeof(Real) +*/
         Buffer::num_2d_vector_mid*m_num_cols*ekat::npack<Spack>(m_num_levs)*sizeof(Spack) +
         Buffer::num_2d_vector_int*m_num_cols*ekat::npack<Spack>(m_num_levs+1)*sizeof(Spack);
}
// =========================================================================================
void SurfaceCouplingExporter::init_buffers(const ATMBufferManager &buffer_manager)
{
  const int nlev_packs       = ekat::npack<Spack>(m_num_levs);
  const int nlevi_packs      = ekat::npack<Spack>(m_num_levs+1);

  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
//  m_buffer.Sa_z = decltype(m_buffer.Sa_z)(mem, m_num_cols);
//  mem += m_buffer.Sa_z.size();
//  m_buffer.Sa_ptem = decltype(m_buffer.Sa_ptem)(mem, m_num_cols);
//  mem += m_buffer.Sa_ptem.size();
//  m_buffer.Sa_dens = decltype(m_buffer.Sa_dens)(mem, m_num_cols);
//  mem += m_buffer.Sa_dens.size();
//  m_buffer.Sa_pslv = decltype(m_buffer.Sa_pslv)(mem, m_num_cols);
//  mem += m_buffer.Sa_pslv.size();
//  m_buffer.Faxa_rainl = decltype(m_buffer.Faxa_rainl)(mem, m_num_cols);
//  mem += m_buffer.Faxa_rainl.size();
//  m_buffer.Faxa_snowl = decltype(m_buffer.Faxa_snowl)(mem, m_num_cols);
//  mem += m_buffer.Faxa_snowl.size();
//  m_buffer.zero_view = decltype(m_buffer.zero_view)(mem, m_num_cols);
//  mem += m_buffer.zero_view.size();

  // 2d views packed views
  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  m_buffer.dz = decltype(m_buffer.dz)(s_mem, m_num_cols, nlev_packs);
  s_mem += m_buffer.dz.size();
  m_buffer.z_mid = decltype(m_buffer.z_mid)(s_mem, m_num_cols, nlev_packs);
  s_mem += m_buffer.z_mid.size();
  m_buffer.z_int = decltype(m_buffer.z_int)(s_mem, m_num_cols, nlevi_packs);
  s_mem += m_buffer.z_int.size();

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);

  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SurfaceCouplingExporter.");
}
// =========================================================================================
void SurfaceCouplingExporter::setup_surface_coupling_data(const SCDataManager &sc_data_manager)
{
  m_num_exports = sc_data_manager.get_num_fields();

  EKAT_ASSERT_MSG(m_num_cols == sc_data_manager.get_field_size(), "Error! Surface Coupling exports need to have size ncols.");

  m_cpl_exports_view_h = decltype(m_cpl_exports_view_h) (sc_data_manager.get_field_data_ptr(),
                                                         m_num_cols, m_num_exports);
  m_cpl_exports_view_d = Kokkos::create_mirror_view(DefaultDevice(), m_cpl_exports_view_h);
}
// =========================================================================================
void SurfaceCouplingExporter::initialize_impl (const RunType /* run_type */)
{
  // Set data

  // Compute column offset and stride

  // Can it be exported

}
// =========================================================================================
void SurfaceCouplingExporter::run_impl (const int /* dt */)
{

}
// =========================================================================================
void SurfaceCouplingExporter::finalize_impl()
{

}
// =========================================================================================
} // namespace scream

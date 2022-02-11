#include "atmosphere_surface_coupling_importer.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
// =========================================================================================
SurfaceCouplingImporter::SurfaceCouplingImporter (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{

}
// =========================================================================================
void SurfaceCouplingImporter::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  Units nondim(0,0,0,0,0,0,0);
  auto Wm2 = W / m / m;
  Wm2.set_string("W/m2)");
  const auto m2 = m*m;

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  FieldLayout scalar2d_layout      { {COL     }, {m_num_cols   } };
  FieldLayout surf_mom_flux_layout { {COL, CMP}, {m_num_cols, 2} };

  add_field<Updated>("sfc_alb_dir_vis",  scalar2d_layout,      nondim, grid_name);
  add_field<Updated>("sfc_alb_dir_nir",  scalar2d_layout,      nondim, grid_name);
  add_field<Updated>("sfc_alb_dif_vis",  scalar2d_layout,      nondim, grid_name);
  add_field<Updated>("sfc_alb_dif_nir",  scalar2d_layout,      nondim, grid_name);
  add_field<Updated>("surf_sens_flux",   scalar2d_layout,      W/m2,   grid_name);
  add_field<Updated>("surf_latent_flux", scalar2d_layout,      W/m2,   grid_name);
  add_field<Updated>("surf_mom_flux",    surf_mom_flux_layout, N/m2,   grid_name);
}
// =========================================================================================
void SurfaceCouplingImporter::setup_surface_coupling_data(const SCDataManager &sc_data_manager)
{
  m_num_imports = sc_data_manager.get_num_fields();

  EKAT_ASSERT_MSG(m_num_cols == sc_data_manager.get_field_size(),
                  "Error! Surface Coupling imports need to have size ncols.");

  m_cpl_imports_view_h = decltype(m_cpl_imports_view_h) (sc_data_manager.get_field_data_ptr(),
                                                         m_num_cols, m_num_imports);
  m_cpl_imports_view_d = Kokkos::create_mirror_view_and_copy(DefaultDevice(),
                                                             m_cpl_imports_view_h);


  std::cout << "IMPORTER: " << m_num_cols << "," << m_num_imports << std::endl;
}
// =========================================================================================
void SurfaceCouplingImporter::initialize_impl (const RunType /* run_type */)
{
  // Loop through fields and compute column information

}
// =========================================================================================
void SurfaceCouplingImporter::run_impl (const int /* dt */)
{

}
// =========================================================================================
void SurfaceCouplingImporter::finalize_impl()
{

}
// =========================================================================================
} // namespace scream

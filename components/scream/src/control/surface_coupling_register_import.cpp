#include "control/surface_coupling.hpp"

#include "share/field/field_utils.hpp"

namespace scream {
namespace control {

void SurfaceCoupling::
register_import(const std::string& fname,
                const int cpl_idx,
                const int vecComp)
{
  // Two separate checks rather than state==Open, so we can print more specific error messages
  EKAT_REQUIRE_MSG (m_state!=RepoState::Clean,
                      "Error! Registration phase has not started yet.\n"
                      "       Did you forget to call set_num_fields(..) ?\n");
  EKAT_REQUIRE_MSG (m_state!=RepoState::Closed, "Error! Registration phase has already ended.\n");

  if (fname == "unused")
  {
    // Do nothing
  } else {
    // Check that we still have room
    EKAT_REQUIRE_MSG(m_num_scream_imports<m_scream_imports_host.extent_int(0),
                       "Error! Imports view is already full. Did you call 'set_num_fields' with the wrong arguments?\n");

    // Get the field, and check that is valid
    Field field = m_field_mgr->get_field(fname);

    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Import field view has not been allocated yet.\n");

    EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

    // Check that this cpl_idx wasn't already registered
    for (int i=0; i<m_num_scream_imports; ++i) {
      EKAT_REQUIRE_MSG(cpl_idx!=m_scream_imports_host(i).cpl_idx,
                       "Error! An import with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
    }

    auto& info = m_scream_imports_host(m_num_scream_imports);

    // Set view data ptr
    info.data = field.get_internal_view_data<Real>();

    // Set cpl index
    info.cpl_idx = cpl_idx;

    // For import fluxes, we must change the sign as cpl and atm interprete the direction differently.
    if (fname == "surf_mom_flux"    || fname == "surf_sens_flux" ||
        fname == "surf_latent_flux" || fname == "surf_lw_flux_up") {
      m_cpl_scream_sign_change_host(m_num_scream_imports) = -1;
    } else {
      m_cpl_scream_sign_change_host(m_num_scream_imports) = 1;
    }

    // Get column offset and stride
    get_col_info (field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

    // Store the identifier of this field, for debug purposes
    m_imports_fids.insert(field.get_header().get_identifier());

    // Update number of imports stored
    ++m_num_scream_imports;
  }
}

} // namespace control
} // namespace scream

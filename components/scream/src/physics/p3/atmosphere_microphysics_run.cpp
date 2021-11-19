#include "physics/p3/atmosphere_microphysics.hpp"

namespace scream {

void P3Microphysics::run_impl (const int dt)
{
  timer.StartTimer();
  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_preproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();
  timer.StopTimer();
  pre_proc_times.push_back(timer.report_time("      pre-proc-time:",get_comm()));

  // Update the variables in the p3 input structures with local values.

  infrastructure.dt = dt;
  infrastructure.it++;

  timer.StartTimer();
  // Reset internal WSM variables.
  workspace_mgr.reset_internals();
  timer.StopTimer();
  wsm_reset_times.push_back(timer.report_time("wsm_reset",get_comm()));

  timer.StartTimer();
  // Run p3 main
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, lookup_tables, workspace_mgr, m_num_cols, m_num_levs);
  timer.StopTimer();
  p3_main_times.push_back(timer.report_time("      p3_main-time:",get_comm()));

  timer.StartTimer();
  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_postproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();
  timer.StopTimer();
  post_proc_times.push_back(timer.report_time("      post-proc-time:",get_comm()));
}

} // namespace scream

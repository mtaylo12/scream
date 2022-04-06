// This is a tiny program that calls p3_init() to generate tables used by p3

#include "physics/p3/p3_f90.hpp"

#include <ekat/ekat_assert.hpp>

int main(int argc, char** argv) {
  EKAT_REQUIRE_MSG (argc==2,
      "Error! p3_generate_tables *requires* to provide the output directory.\n");

  bool write_tables = true;
  auto output_dir = argv[1];
  scream::p3::p3_init(write_tables, output_dir);

  return 0;
}

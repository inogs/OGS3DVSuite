#include <sys/stat.h>

#include <filesystem>

#include "OGS/OGS3DVSuite.h"
#include "globals.h"
#include "gtest/gtest.h"

TEST(Simulation, ComplainIfOGSFileIsMissing) {
  std::string source_dir(SOURCE_DIR);
  EXPECT_THROW(std::make_unique<OGS::Simulation>(
                   source_dir + "/tests/data/test_grid/not_found.ogs");
               , std::ios_base::failure);
}



TEST(Simulation, ReadMainRectilinearFile) {
  std::string source_dir(SOURCE_DIR);
  const auto simulation = std::make_unique<const OGS::Simulation>(
      source_dir + "/tests/data/test_grid/rectilinear.ogs");
  ASSERT_EQ(simulation->n_meshes(), 2);
}



TEST(Simulation, ReadMainStructuredFile) {
  std::string source_dir(SOURCE_DIR);
  const auto simulation = std::make_unique<const OGS::Simulation>(
      source_dir + "/tests/data/test_grid/struct2d.ogs");
  ASSERT_EQ(simulation->n_meshes(), 2);
}



TEST(Simulation, get_mesh_index) {
  std::string source_dir(SOURCE_DIR);
  const auto simulation = std::make_unique<const OGS::Simulation>(
      source_dir + "/tests/data/test_grid/rectilinear.ogs");
  ASSERT_EQ(simulation->get_mesh_index("Mercator"), 0);
  ASSERT_EQ(simulation->get_mesh_index("Cylindrical"), 1);

  EXPECT_THROW(int err_value = simulation->get_mesh_index("missing"),
               std::invalid_argument);
}



TEST(Simulation, check_loading_mesh) {
  std::string source_dir(SOURCE_DIR);
  const auto simulation = std::make_unique<OGS::Simulation>(
      source_dir + "/tests/data/test_grid/rectilinear.ogs");
  EXPECT_THROW(simulation->load_mesh(2), std::invalid_argument);

  simulation->load_mesh(1);
  simulation->load_mesh(0);
}



TEST(Simulation, var_path) {
  std::string source_dir(SOURCE_DIR);
  const auto simulation = std::make_unique<OGS::Simulation>(
      source_dir + "/tests/data/test_grid/rectilinear.ogs");
  EXPECT_THROW(simulation->load_mesh(2), std::invalid_argument);

  struct stat sb {};

  const unsigned int n_timesteps = simulation->ntsteps();
  for (const auto var_type : OGS::AllVarTypes) {
    const int n_vars = simulation->var_n(var_type);
    for (unsigned int i = 0; i < n_vars; ++i)
      for (unsigned int t = 0; t < n_timesteps; ++t) {
        const std::filesystem::path var_path =
            simulation->var_path(var_type, i, t);
        EXPECT_TRUE(stat(var_path.c_str(), &sb) == 0 &&
                    !(sb.st_mode & S_IFDIR));
      }
  }
}



TEST(Simulation, check_loading_structured_mesh) {
  std::string source_dir(SOURCE_DIR);
  const auto simulation = std::make_unique<OGS::Simulation>(
      source_dir + "/tests/data/test_grid/struct2d.ogs");
  EXPECT_THROW(simulation->load_mesh(2), std::invalid_argument);

  simulation->load_mesh(1);
  simulation->load_mesh(0);
}

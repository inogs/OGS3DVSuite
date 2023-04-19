#include "OGS/OGS3DVSuite.h"
#include "globals.h"
#include "gtest/gtest.h"

TEST(RectilinearMeshVertexCoords, ReadFile) {
  std::string source_dir(SOURCE_DIR);
  const auto mesh_vertex_coords =
      OGS::RectilinearMeshVertexCoords::create_rectilinear_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/rectilinear.ogsmsh.Mercator");
  EXPECT_EQ(mesh_vertex_coords->nlon(), 395);
  EXPECT_EQ(mesh_vertex_coords->nlat(), 161);
  EXPECT_EQ(mesh_vertex_coords->nlev(), 44);

  for (const auto mask_type : OGS::AllMaskTypes)
    EXPECT_EQ(mesh_vertex_coords->ncells(),
              mesh_vertex_coords->get_mask(mask_type).get_n());
}



TEST(RectilinearMeshVertexCoords, CoordsOnLines) {
  std::string source_dir(SOURCE_DIR);
  const auto mesh_vertex_coords =
      OGS::RectilinearMeshVertexCoords::create_rectilinear_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/rectilinear.ogsmsh.Cylindrical");

  const auto hline = mesh_vertex_coords->get_horizontal_grid_line(2, 5);
  const auto vline = mesh_vertex_coords->get_vertical_grid_line(12, 9);
  const auto column = mesh_vertex_coords->get_column_grid_line(100, 100);

  EXPECT_EQ(std::get<0>(*hline).size(), mesh_vertex_coords->nlon());
  EXPECT_EQ(std::get<1>(*hline).size(), mesh_vertex_coords->nlon());
  EXPECT_EQ(std::get<2>(*hline).size(), mesh_vertex_coords->nlon());

  EXPECT_EQ(std::get<0>(*vline).size(), mesh_vertex_coords->nlat());
  EXPECT_EQ(std::get<1>(*vline).size(), mesh_vertex_coords->nlat());
  EXPECT_EQ(std::get<2>(*vline).size(), mesh_vertex_coords->nlat());

  EXPECT_EQ(std::get<0>(*column).size(), mesh_vertex_coords->nlev());
  EXPECT_EQ(std::get<1>(*column).size(), mesh_vertex_coords->nlev());
  EXPECT_EQ(std::get<2>(*column).size(), mesh_vertex_coords->nlev());

  for (unsigned int i = 1; i < mesh_vertex_coords->nlev(); ++i)
    EXPECT_GT(std::get<2>(*column).at(i), std::get<2>(*column).at(i - 1));
}



TEST(Struct2DMeshVertexCoords, ReadFile) {
  std::string source_dir(SOURCE_DIR);
  const auto mesh_vertex_coords =
      OGS::Struct2DMeshVertexCoords::create_struct2d_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/struct2d.ogsmsh.Cylindrical.nc");
  EXPECT_EQ(mesh_vertex_coords->nlon(), 395);
  EXPECT_EQ(mesh_vertex_coords->nlat(), 161);
  EXPECT_EQ(mesh_vertex_coords->nlev(), 44);

  for (const auto mask_type : OGS::AllMaskTypes)
    EXPECT_EQ(mesh_vertex_coords->ncells(),
              mesh_vertex_coords->get_mask(mask_type).get_n());
}



TEST(Struct2DMeshVertexCoords, CoordsOnLines) {
  std::string source_dir(SOURCE_DIR);
  const auto mesh_vertex_coords =
      OGS::Struct2DMeshVertexCoords::create_struct2d_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/struct2d.ogsmsh.Mercator.nc");

  const auto hline = mesh_vertex_coords->get_horizontal_grid_line(2, 5);
  const auto vline = mesh_vertex_coords->get_vertical_grid_line(12, 9);
  const auto column = mesh_vertex_coords->get_column_grid_line(100, 100);

  EXPECT_EQ(std::get<0>(*hline).size(), mesh_vertex_coords->nlon());
  EXPECT_EQ(std::get<1>(*hline).size(), mesh_vertex_coords->nlon());
  EXPECT_EQ(std::get<2>(*hline).size(), mesh_vertex_coords->nlon());

  EXPECT_EQ(std::get<0>(*vline).size(), mesh_vertex_coords->nlat());
  EXPECT_EQ(std::get<1>(*vline).size(), mesh_vertex_coords->nlat());
  EXPECT_EQ(std::get<2>(*vline).size(), mesh_vertex_coords->nlat());

  EXPECT_EQ(std::get<0>(*column).size(), mesh_vertex_coords->nlev());
  EXPECT_EQ(std::get<1>(*column).size(), mesh_vertex_coords->nlev());
  EXPECT_EQ(std::get<2>(*column).size(), mesh_vertex_coords->nlev());

  for (unsigned int i = 1; i < mesh_vertex_coords->nlev(); ++i)
    EXPECT_GT(std::get<2>(*column).at(i), std::get<2>(*column).at(i - 1));
}



TEST(Struct2DVsRectilinear, CompareLevels) {
  std::string source_dir(SOURCE_DIR);
  const auto rectilinear =
      OGS::RectilinearMeshVertexCoords::create_rectilinear_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/rectilinear.ogsmsh.Mercator");
  const auto struct2d =
      OGS::Struct2DMeshVertexCoords::create_struct2d_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/struct2d.ogsmsh.Mercator.nc");

  const unsigned int nlev = rectilinear->nlev();
  ASSERT_EQ(nlev, struct2d->nlev());

  std::vector c1_rectilinear =
      std::get<2>(*rectilinear->get_column_grid_line(42, 43));
  std::vector c1_struct2d =
      std::get<2>(*struct2d->get_column_grid_line(42, 43));

  std::vector c2_rectilinear =
      std::get<2>(*rectilinear->get_column_grid_line(194, 61));
  std::vector c2_struct2d =
      std::get<2>(*struct2d->get_column_grid_line(194, 61));

  for (unsigned int i = 0; i < nlev; ++i) {
    EXPECT_FLOAT_EQ(c1_rectilinear[i], c1_struct2d[i]);
    EXPECT_FLOAT_EQ(c2_rectilinear[i], c2_struct2d[i]);
  }
}



TEST(Struct2DVsRectilinear, CompareMasks) {
  std::string source_dir(SOURCE_DIR);
  const auto rectilinear =
      OGS::RectilinearMeshVertexCoords::create_rectilinear_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/rectilinear.ogsmsh.Cylindrical");
  const auto struct2d =
      OGS::Struct2DMeshVertexCoords::create_struct2d_mesh_vertex_coords(
          source_dir + "/tests/data/test_grid/struct2d.ogsmsh.Cylindrical.nc");

  const unsigned int ncells = rectilinear->ncells();
  ASSERT_EQ(ncells, struct2d->ncells());

  std::vector<unsigned int> n_of_components(3);
  n_of_components[OGS::BASINS] = 16;
  n_of_components[OGS::COASTS] = 1;
  n_of_components[OGS::LAND] = 1;

  for (const auto mask_type : OGS::AllMaskTypes)
    for (unsigned int i = 0; i < ncells; ++i)
      for (unsigned int j = 0; j < n_of_components[mask_type]; ++j)
        EXPECT_EQ(rectilinear->get_mask(mask_type)[i][j],
                  struct2d->get_mask(mask_type)[i][j])
            << "Error with mask " << OGS::masktype2string(mask_type)
            << " on index (" << i << ", " << j << ");";
}

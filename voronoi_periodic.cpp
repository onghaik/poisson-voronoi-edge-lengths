#include "voro++/voro++.hh"
#include <cmath>
#include <fstream>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

typedef std::vector<double> Point;

double distance(const Point &p1, const Point &p2) {
  double dx = p1[0] - p2[0];
  double dy = p1[1] - p2[1];
  double dz = p1[2] - p2[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to generate random points in 3D space using Poisson point process
std::vector<Point> generate_poisson_points(double intensity, double x_min,
                                           double x_max, double y_min,
                                           double y_max, double z_min,
                                           double z_max, unsigned int seed) {
  int num_points;
  std::mt19937 gen;
  gen.seed(seed);

  double vol = (x_max - x_min) * (y_max - y_min) * (z_max - z_min);

  // Poisson distribution with mean intensity * box_size^3
  std::poisson_distribution<int> poisson(intensity * vol);
  num_points = poisson(gen);

  std::vector<Point> points;
  std::uniform_real_distribution<double> x_dis(x_min, x_max);
  std::uniform_real_distribution<double> y_dis(y_min, y_max);
  std::uniform_real_distribution<double> z_dis(z_min, z_max);

  for (int i = 0; i < num_points; ++i) {
    Point point = {x_dis(gen), y_dis(gen), z_dis(gen)};
    points.push_back(point);
  }

  return points;
}

int main() {
  // Simulation parameters
  const double intensity = 1.0; // Intensity parameter for Poisson point process
  // const unsigned int seed = 2;

  // Set up constants for the container geometry
  const double x_min = -1, x_max = 1;
  const double y_min = -1, y_max = 1;
  const double z_min = -1, z_max = 1;
  const double cvol = (x_max - x_min) * (y_max - y_min) * (x_max - x_min);

  // Set up the number of blocks that the container is divided into
  const int n_x = 6, n_y = 6, n_z = 6;

  // Output file to save edge lengths
  std::string edge_lengths_file = "edge_lengths.txt";
  std::string edge_statistics_file = "edge_statistics.txt";

  std::vector<double> edge_lengths;
  int nEdges = 0;

  for (unsigned int seed = 0; seed < 50; seed++) {

    // Generate Poisson-distributed random points
    std::vector<Point> points = generate_poisson_points(
        intensity, x_min, x_max, y_min, y_max, z_min, z_max, seed);

    // Compute the periodic Voronoi diagram and save edge lengths
    voro::container con(x_min, x_max, y_min, y_max, z_min, z_max, n_x, n_y, n_z,
                        true, true, true, 8);

    // Add points to the container with unique IDs
    for (int i = 0; i < points.size(); i++) {
      con.put(i, points[i][0], points[i][1], points[i][2]);
    }

    std::vector<double> distances;
    voro::voronoicell_neighbor con_cell;
    /* double total_edge_distance = 0.0; */

    voro::c_loop_all vl(con);

    if (vl.start()) {
      do {
        if (con.compute_cell(con_cell, vl)) {
          nEdges += con_cell.number_of_edges();
          /* total_edge_distance += con_cell.total_edge_distance(); */
          for (int i = 0; i < con_cell.p; i++) {
            for (int j = 0; j < con_cell.nu[i]; j++) {
              int k = con_cell.ed[i][j];
              if (k > i) {
                double dx = con_cell.pts[k << 2] - con_cell.pts[i << 2];
                double dy =
                    con_cell.pts[(k << 2) + 1] - con_cell.pts[(i << 2) + 1];
                double dz =
                    con_cell.pts[(k << 2) + 2] - con_cell.pts[(i << 2) + 2];
                double dst = 0.5 * sqrt(dx * dx + dy * dy + dz * dz);
                edge_lengths.push_back(dst);
              }
            }
          }
          /*   std::cout << "Sum of edge distances vector: " */
          /*             << std::accumulate(distances.begin(), distances.end(),
           * 0.0) */
          /*             << std::endl; */
          /*   std::cout << "Total edge distances: " << total_edge_distance */
          /*             << std::endl; */
          /*   std::cout << "Difference: " */
          /*             << std::accumulate(distances.begin(), distances.end(),
           * 0.0) - */
          /*                    total_edge_distance */
          /*             << std::endl; */
        }
      } while (vl.inc());
    }
  }

  // Save the edge lengths to a .txt file
  std::ofstream lengths_outfile(edge_lengths_file);
  if (lengths_outfile) {
    for (const auto &edge : edge_lengths) {
      lengths_outfile << edge << std::endl;
    }
    std::cout << "Poisson-Voronoi tessellation with periodic boundary "
                 "conditions completed.\n"
              << "Lengths of " << nEdges << " edges saved to "
              << edge_lengths_file << std::endl;
  } else {
    std::cerr << "Error: Unable to open the output file." << std::endl;
  }

  std::ofstream stats_outfile(edge_statistics_file);
  if (stats_outfile) {
    stats_outfile << "Mean: "
                  << gsl_stats_mean(edge_lengths.data(), 1, edge_lengths.size())
                  << '\n'
                  << "Variance: "
                  << gsl_stats_variance(edge_lengths.data(), 1,
                                        edge_lengths.size())
                  << std::endl;
  }

  // con.draw_particles("voronoi_periodic_p.gnu");

  // con.draw_cells_gnuplot("voronoi_periodic_v.gnu");

  return 0;
}

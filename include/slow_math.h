#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

/* calculate raw rmse */
double root_mean_squared_error( std::vector< double > const &v0, 
           std::vector< double > const &v1 );

/* rmse comparison function */
bool rmse_based_comparison( std::vector< double > const &v0, 
                      std::vector< double > const &v1,
                      double const tolerance );

/* generate a row major rotation matrix that rotates "angle" about the y-axis */
std::vector< double > generate_3d_rotation_matrix( double radians_around_y );

/* dot product function for row major matrix */
std::vector< double > slow_dot_product( std::vector< double > const &matrix,
                                        std::vector< double > const &vector,
                                        int rows,
                                        int cols );

std::vector< double > slow_cross_product( std::vector< double > const &v0,
                                          std::vector< double > const &v1 );

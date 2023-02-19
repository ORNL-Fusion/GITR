#include "data.hpp"
#include "tests_general.hpp"

/* test the data component with relevant floating point data types */
TEMPLATE_TEST_CASE("data component", "container tests", double, float )
{
  /* The test works by calling "get_vector()" on each addressable 1d vector in order and
     concatenating the result. If "get_vector()" works right, the final result will equal
     the original vector */
  SECTION("tensor::get_vector() test")
  {
    int const dim_0 = 2;
    int const dim_1 = 3;
    int const dim_2 = 4;
    int const dim_3 = 5;
    int const dim_4 = 3;

    std::vector< TestType > numbers( dim_0 * dim_1 * dim_2 * dim_3 * dim_4 );

    for( int i = 0; i < static_cast<int>(numbers.size()); ++i ) numbers[ i ] = i;

    data::tensor< TestType, dim_0, dim_1, dim_2, dim_3, dim_4 > t( numbers.data() );

    std::vector< TestType > v_gold;

    for( int i = 0; i < dim_0; ++i )
    {
      for( int ii = 0; ii < dim_1; ++ii )
      {
        for( int iii = 0; iii < dim_2; ++iii )
        {
          for( int iiii = 0; iiii < dim_3; ++iiii )
          {
            std::vector< TestType > v = t.get_vector( i, ii, iii, iiii ).get_std_vector();

            v_gold.insert( v_gold.end(), v.begin(), v.end() );
          }
        }
      }
    }

    REQUIRE( v_gold == numbers );
  }
}

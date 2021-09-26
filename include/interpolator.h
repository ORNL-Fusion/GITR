#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

double interpolate_basis_cell( std::vector< double > const &basis_cell,
                               std::vector< double > const &coordinates,
                               std::vector< int > const &d_len,
                               std::vector< double > const &max_values )
{
  std::cout << "Ahoy, Captain! Inside interpolate_basis_cell" << std::endl;

  assert( coordinates.size() == d_len.size() );

  assert( coordinates.size() == max_values.size() );

  std::vector< double > spacing( max_values.size() );

  for( int i = 0; i < spacing.size(); i++ )
  {
    spacing[ i ] = max_values[ i ] / double(d_len[ i ]);
    std::cout << "spacing[ i ] = " << spacing[ i ] << std::endl;
  }
  std::cout << std::endl;

  /* each basis component will have 2 numbers associated with it */
  std::vector< double > normalized_fractions( coordinates.size() * 2 );

  std::cout << "fractions:" << std::endl;

  for( int i = 0; i < coordinates.size(); i++ )
  {
    /* high fraction, matches with a 1 bit */
    normalized_fractions[ i * 2 + 1 ] =
    ( coordinates[ i ] - ( std::floor( coordinates[ i ] / spacing[ i ] ) * spacing[ i ] ) ) /
    spacing[ i ];

    /* show the calculation of the fractions */
    std::cout << "coordinates[ " << i << " ] = " << coordinates[ i ]
              << " base_point[ " << i << " ] = "  
              << ( std::floor( coordinates[ i ] / spacing[ i ] ) * spacing[ i ] )
              << " difference = " 
              << ( coordinates[ i ] - ( std::floor( coordinates[ i ] / spacing[ i ] ) 
                   * spacing[ i ] ) )
              << std::endl;

    /* low fraction, matches with a 0 bit */
    normalized_fractions[ i * 2 ] = 1 - normalized_fractions[ i * 2 + 1 ];

    std::cout << "low " << i << " = " << normalized_fractions[ i * 2 ]
              << " high " << i << " = " << normalized_fractions[ i * 2 + 1 ] << std::endl;
  }
  std::cout << std::endl;

  /* we have the basis cell, what else is needed? */
  double sum = 0;

  for( int i = 0; i < basis_cell.size(); i++ )
  {
    double product = 1;

    std::cout << "product on term " << i << ":" << std::endl;

    /* iterate the bits in "i" */
    for( int j = 0; j < coordinates.size(); j++ )
    {
      int index = j * 2 + ( ( i & ( 0x1 << j ) ) >> j );

      std::string which = ( index % 2 ? "high" : "low" );

      std::cout << which << "*";

      /* Captain! Should I calculate instead of store? */
      /* Captain! Do i >> j & 0x1 instead for a shorter expression, save a clock cycle */
      product *= normalized_fractions[ j * 2 + ( ( i & ( 0x1 << j ) ) >> j ) ];
    }
    std::cout << std::endl;
    std::cout << "product = " << product << std::endl;
    std::cout << "basis_cell[ " << i << " ] = " << basis_cell[ i ] << std::endl;

    sum += product * basis_cell[ i ];

    std::cout << "sum = " << sum << std::endl << std::endl;
  }

  return sum;
}

/*

What do you need to do next? Next, we need to figure out how to make more tests.

*/
/* can you turn this into a variadic template? If all of them are you can create a header only
   library */
/* d_len is the length of each dimension */
/* max values is the biggest x-value there is */
/* data is the whole data field */
/* coordinates are single dimensional float values */
std::vector< double > fetch_basis_cell( std::vector< double > const &data,
                                        std::vector< double > const &coordinates,
                                        std::vector< int > const &d_len,
                                        std::vector< double > const &max_values )
{
  assert( coordinates.size() == d_len.size() );

  assert( coordinates.size() == max_values.size() );

  std::vector< double > spacing( max_values.size() );

  for( int i = 0; i < spacing.size(); i++ )
  {
    spacing[ i ] = max_values[ i ] / double(d_len[ i ]);
  }

  std::vector< int > d_multipliers( d_len.size(), 1 );

  for( int i = 0; i < d_multipliers.size() - 1; i++ )
  {
    d_multipliers[ i + 1 ] = d_multipliers[ i ] * d_len[ i ];
  }

  int base_index = 0;

  for( int i = 0; i < coordinates.size(); i++ )
  {
    std::cout << "coordinates[ " << i << " ] = " << coordinates[ i ]
              << " spacing[ " << i << " ] = " << spacing[ i ]
              << " d_len[ " << i << " ] = " << d_multipliers[ i ] << std::endl;

    std::cout << "base_index before = " << base_index << std::endl;

    base_index += ( std::floor( coordinates[ i ] / spacing[ i ] ) * d_multipliers[ i ] );

    std::cout << "base_index after = " << base_index << std::endl;
  }

  std::cout << "fetch base_index: " << base_index << std::endl;

  std::vector< double > basis_cell( 1 << d_len.size() );

  for( int i = 0; i < basis_cell.size(); i++ )
  {
    int flat_index = base_index;

    for( int j = 0; j < coordinates.size(); j++ )
    {
      flat_index += ( ( ( i & ( 0x1 << j ) ) >> j ) * d_multipliers[ j ] );
    }

    basis_cell[ i ] = data[ flat_index ];

    std::cout << "basis_cell[ " << i << " ] = "
              << "haystack[ " << flat_index << " ]" << std::endl;
  }

  return basis_cell;
}

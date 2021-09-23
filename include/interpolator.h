#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
/* can you turn this into a variadic template? If all of them are you can create a header only
   library */
/* d_len is the length of each dimension */
/* max values is the biggest x-value there is */
/* data is the whole data field */
/* coordinates are single dimensional float values */
std::vector< double > fetch_basis_cell( std::vector< double > data,
                                        std::vector< double > coordinates,
                                        std::vector< int > d_len,
                                        std::vector< double > max_values )
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

  /* Captain! There's an off-by-1 error in this base_index calculation code. Where? */
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

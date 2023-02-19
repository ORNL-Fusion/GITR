/*

This component defines the 2 basic data structures in the sample library.

vector_1d:

  An ordered, linear array of data

tensor:

  An ordered, linear array of data indexed with an ordered n-tuple.

*/

#pragma once
#include <iostream>
#include <cassert>
#include <array>
#include <vector>

/* Ahoy, Captain! Create the time class here */

/* all components belong in an appropriately-titled namespace */
namespace data
{

/*
  An ordered, linear array of data.
  Optionally manages memory.

T = data type
N = length
*/
/* Improvement: use size_t instead of "int" for container sizes */
template< typename T, int N >
class vector_1d
{
  public:

    /* if constructed with empty parameter list, allocate from the heap and delete when done */
    vector_1d( T *data = nullptr );
    ~vector_1d();

    /* access and assignment operators */
    T operator []( int const i ) const { assert( i < size() && i >= 0 ); return data[ i ]; }
    T& operator []( int const i ) { assert( i < size() && i >= 0 ); return data[ i ]; }

    /* useful function to copy data into a std::vector and return it */
    std::vector< T > get_std_vector(); 

    int size() const { return N; }

  private:

    T *data;
    bool manage;
};

/* constructor */
template< typename T, int N >
vector_1d< T, N >::vector_1d( T *data )
  :
  data( data )
{
  /* conditionally allocate memory if the user does not supply a pointer */
  if( data == nullptr )
  {
    manage = true;
    data = new T[ N ];
  }
  else manage = false;

  return;
}

/* destructor */
template< typename T, int N >
vector_1d< T, N >::~vector_1d()
{
  /* conditionally free memory if the user indicated that this object manage the memory */
  if( manage == true )
  {
    delete[] data;
  }
  return;
}

template< typename T, int N >
std::vector< T > vector_1d< T, N >::get_std_vector()
{
  return std::vector< T >( data, data + N );
}

/* Ahoy, Captain! vector_1d trash above */



/* for variadic templates, it is of interest to obtain the last value in a non-type parameter
   pack. The following 2 classes implement this */
template< int T >
class last_parameter
{
public:
  static int const val = T;
};

template< int ... dimlist >
class get_last_parameter
{
public:

  /* unary right fold expression: the parameter pack "dimlist" is expanded across the fold
     expression. The i-th element of parameter pack "p" is denoted as "p[ i ]". Fold expression
     expands like:

     static int const val =
     ( 
       last_parameter< dimlist[ 0 ] >{}, 
       last_parameter< dimlist[ 1 ] >{}, 
       ...
       last_parameter< dimlist[ n ] >{}

     ).val

     The comma operator "," evaluates the expression and throws away the result. So, the above
     reduces to:

     static int const val = last_parameter< dimlist[ n ] >{}.val;
  */
  static int const val = ( last_parameter< dimlist >{}, ... ).val;
};

/*

An ordered, linear array of data indexed with an ordered tuple. The minimum addressable units
for this object are vector_1d objects. This class does not manage memory.

T = data type
dimlist = ordered index tuple 

  Given a dimension specification:

  dimlist "d" = { d_i } = { d_0, d_1, ... , d_n }

  vector_1d objects are accessed via an index tuple:

  index tuple "t" = { t_i } = { t_1, t_2, ... , t_n } 

  0 <= t_i < d_i for i = { 1, ... , n }

  and is evaluated as: 

  linear position =
  t_0 * product( { d_1, ... , d_n } ) +
  t_1 * product( { d_2, ... , d_n } ) +
  t_2 * product( { d_3, ... , d_n } ) +
  ... +
  t_m * d_n

*/
/* Future improvement: add data management as an option for this class */
template< typename T, int ... dimlist >
class tensor
{
  /* functions */
  public:

    /* constructor takes pointer to data managed elsewhere */
    tensor( T *data );

    /* accesses a vector_1d according to process outlined in class description. The ordered index
       tuple is specified element-by-element as a variadic non-type parameter pack template
       argument */
    template< typename ... pack >
    vector_1d< T, get_last_parameter< dimlist ... >::val > get_vector( pack ... ) const;

  /* functions */
  private:

    /* The next 2 function overloads implement core functionality of the public data access
       routine above */
    template< typename ... indices_remaining >
    int get_vector_recurse( int dim, int index, indices_remaining ... ) const;

    /* for the final element */
    template< typename ... indices_remaining >
    int get_vector_recurse( int dim, int index ) const;

  /* data members */
  private:

    /* dimension specification { d_i } from description */
    std::array< int, sizeof...( dimlist ) > dims_array;

    /* Linear array of data values */
    T *data;

    /* used by "get_vector" functions to address a vector_1d */
    std::array< int, sizeof...( dimlist ) - 1 > offset_factors;
};


/* this function overload defines a base case for the compiler as it recursively builds the
   public "get_vector" routine according to the metaprogrammed template definition */
template< typename T, int ... dimlist > /* class template's parameter list */
template< typename ... indices_remaining > /* member function template's paramater list */
int /* function return value */
tensor< T, dimlist ... >:: /* scope specification */
get_vector_recurse( int dim, int index ) /* function name and parameter list */
const /* cv-qualifier */
{
  return index * offset_factors[ dim ];
}

/* The compiler will unroll the recursion. Note that this function accepts 3 arguments, yet it
   invokes itself with only 2. Each invokation shifts the leading argument of the parameter pack
   to the position called "index" in the parameter list, and puts the remaining arguments as a
   parameter pack in the next position in the parameter list. For finite-length lists, the
   remainder will eventually be empty and the base-case overload defined above will be called.
   This is compile-time recursion, not runtime recursion, so no recursive calls will occur during
   runtime. */
template< typename T, int ... dimlist > /* class template's parameter list */
template< typename ... indices_remaining > /* member function template's parameter list */
int /* function return value */
tensor< T, dimlist ... >:: /* scope specification */
get_vector_recurse( int dim, int index, indices_remaining ... indices ) /* function name/params */
const /* cv-qualifier */
{
  /* */
  return index * offset_factors[ dim ] + get_vector_recurse( dim + 1, indices ... );
}

/* given an ordered tuple, return the requested 1d vector */
/* performs compile-time checks and calls recursive implementation defined above */
template< typename T, int ... dimlist >
template< typename ... pack >
vector_1d< T, get_last_parameter< dimlist ... >::val >
tensor< T, dimlist ... >::
get_vector( pack ... p ) const
{
  /* Captain! Remove the -1 from this check? */
  /* validate tuple's cardinality */
  assert( sizeof...( p ) == dims_array.size() - 1 );

  /* This fold expression expands to a boolean statement - ensures correct data types */
  assert( ( std::is_same_v< int, pack > && ... ) );

  int offset = get_vector_recurse( 0, p ... );

  /* Ahoy, Captain! You can simply return the pointer at the given offset here? */
  /* Ahoy, Captain! I think the next thing you have to do is just return element at
     data + offset + get_lat_parameter<>... */

  return vector_1d< T, get_last_parameter< dimlist ... >::val >( data + offset );
}

/* constructor */
template< typename T, int ... dimlist >
tensor< T, dimlist ... >::tensor( T *data_arg )
  :
  dims_array{ dimlist ... },
  data( data_arg )
{
  int factor = 1;

  for( int i = dims_array.size() - 1; i > 0; --i )
  {
    factor *= dims_array[ i ];

    offset_factors[ i - 1 ] = factor;
  }

  return;
}

} // end namespace data

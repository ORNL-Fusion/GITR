#include <cassert>
#include <array>
#include <iostream>
#include <cmath>

/*
  Simple fsm. The "spin()" routine increments this object to the next state
  S = number of states
*/
template< int S >
class round_robin_wheel
{
public:

  round_robin_wheel( int const start_index = 0 ) : size( S ), current_index( start_index ) {}

  void set_index( int const new_index )
  { 
    assert( new_index < size && new_index >= 0 );
    current_index = new_index;
    return; 
  }

  int spin()
  {
    int const n = current_index++;

    if (current_index == size)
      current_index = 0;

    return n;
  }

private:

  int const size;
  /* FSM current state */
  int current_index;
};

/*

This FSM is an ordered tuple of identical round_robin_wheel FSMs.

X = Number of round_robin_wheel FSMs that generate the state-tuple for this FSM
Y = value of the final state for each round_robin_wheel

Finite states are represented as tuples:  
{ ... , t_(i-1), t_i, t_(i+1), ... }. Where the value of each t_i is determined by the state of
an object of type "round_robin_wheel" described above. t_i progresses to its next state when
t_(i+1) completes a full cycle; when t_i completes a full cycle, t_( i - 1 ) progress one state
further.

*/
template< int X, int Y >
class round_robin_nd
{
public:

  /* default initial state { t_i = 0 for all i } */
  round_robin_nd( std::array< int, X > const indices );

  round_robin_nd() = default;

  /* arbitrarily change state */
  void set_indices( std::array< int, X > const &new_indices );

  /* obtain current state */
  std::array< int, X > get_indices() { return indices; }

  /* progress the FSM to the next state */
  void spin();

  /* return the FSM to the previous state */
  void back_spin();

private:

  /* current FSM state */
  std::array< int, X > indices{};

  /* generates next or previous FSM state */
  std::array< round_robin_wheel< Y >, X > wheels;
};

template< int X, int Y >
round_robin_nd< X, Y >::round_robin_nd( std::array< int, X > const indices )
{
  assert( X > 0 );

  set_indices( indices );

  return;
}

template< int X, int Y >
void
round_robin_nd< X, Y >::spin()
{
  for( int i = 0; i < static_cast<int>( indices.size() ); ++i )
  {
    int index = wheels[ i ].spin();
    indices[ i ] = index;
    if( index != 0 ) return;
  }

  return;
}

/* Captain! Replace these returns with break? */
template< int X, int Y >
void
round_robin_nd< X, Y >::back_spin()
{
  for( int i = indices.size() - 1; i >= 0; --i )
  {
    int index = wheels[ i ].spin();
    indices[ i ] = index;
    if( index != 0 ) return;
  }

  return;
}

template< int X, int Y >
void round_robin_nd< X, Y >::set_indices( std::array< int, X > const &new_indices )
{
  indices = new_indices;

  for( int i = 0; i < static_cast<int>(indices.size()); ++i )
  {
    wheels[ i ].set_index( indices[ i ] );
    wheels[ i ].spin();
  }

  return;
}

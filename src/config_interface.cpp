#include "config_interface.h"

std::string config_interface::ahoy() const
{
  return "Ahoy!";
}
/*

What are you doing right now? Setting up ctest to run some ctests. Adding the CLI module you
found, creating tests for all of it

1. Reader object, derived from base class
2. The config object must be constructed with a reader object.
3. Validate and get everything you need in the constructor.
4. Constructors for config modules use the reader object too
5. Getters and setters get and set the config stuff

Questions:

1. How can I read in all variables from a libconfig file? Looks like you can use getLength()
   and iterate. But, is there some way to get the string name of each one?


  Looks like we have an isRoot() method, and a getName. And we have iterators? We can get to
  the bottom of the tree I think. Also a getLength(). We need to get all the names for sure.

  Top level deals with the root object, bottom levels deal with anything below


2. How can I get command line variables set up? Add some basic command line options
   Input directory, output directory, various inputs like config filenames? Or file for
   the root config

3. How can I get a CI setup for GITR using docker images etc? Do this one tomorrow.

*/

/*

command line parser
libconfig backend - receive a mapping and return values based on the mapping

config
config module

How to test this:

1. Setup a command line parser test. Add option for config filename and output directory.
2. ste up a libconfig backend test - read modules and dump keys in a flat manner. Might need
   a smaller example.cfg file. Should probably just adapt one of the test modules from
   libconfig
3. set up a configuration test - define a config module and a config submodule. Get the names
   via scoped enums

4. apply the new functionality to GITR

5. test it via the system test GITR_processing repo. Change that to GITR_system_test repo?

*/

#include <iostream>
#include <string>
#include "libconfig.h++"

/*
Could I make a const class, const constructed? Can a const class have virtual data members?
You should have const accessors to mutable data instead. This class is the flattener. It will
provide 
*/
class libconfig_data_reader
{
};

/* Inherits from a class that defines a mapping etc? Subclass to work with the data reader */
/* The base class should work with the API that the config_interface class will use */
/* subclass once and simply define a constructor that sets the variables a certain way */
/* Include a map of < scoped enum > : config_module that can return a sub-config module */
/* overload the get() method based on what type of scoped enum it is called with */
/* this effectively embeds the type information into the scoped enum string calls */
/* does it effectively add type to an int?? Yep, that's the whole point in fact */
class config_module
{
};

/* Creates a list of config_module that can be queried for objects */
class config_interface
{
  public:

  std::string ahoy() const;
};

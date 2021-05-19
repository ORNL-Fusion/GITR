#include <iostream>
#include <string>
#include <memory>
#include "libconfig.h++"

class libconfig_string_query
{
  public:

  libconfig_string_query( std::string libconfig_file );

  template< typename T >
  void operator()( std::string const query_key, T &query_value ) const
  {
    bool success = cfg.lookupValue( query_key, query_value );

    if( success == false ) 
    {
      std::cout << "invalid query key: " << query_key << std::endl;
      exit(0);
    }
  }

  private:

  libconfig::Config cfg;
};

class config_module_base
{
  public:

  config_module_base( std::string module_path = "" );

  std::string const &module_path() { return module_path }

  /* get config sub-module */
  virtual std::shared_ptr< config_module_base > get( int key );
  /* get a config value */
  virtual void get( int key, double &val );
  virtual void get( int key, float &val );
  virtual void get( int key, int &val );
  virtual void get( int key, std::string &val );

  private:

  std::string module_path;
};

class impurity_particle_source final : public config_module_base
{
  public:

  /* get config sub-module */
  std::shared_ptr< config_module_base > get( int key );

  /* get a config value */
  void get( int key, double &val ) override;
  void get( int key, float &val ) override;
  void get( int key, int &val ) override;
  void get( int key, std::string &val ) override;
  /* define strongly typed enums as well - set them to int? */

  private:

  class libconfig_string_query const &query;
};

class gitr_config_interface : public config_module_base
{
  public:

  /* define the mapping */
  gitr_config_interface();

  /* define scoped enums to refer to the various config modules */

  std::string ahoy() const;
};

/*

Ok... how about instead of that... we just put the root setting in there?

Goals: separate libconfig from this. 


























*/

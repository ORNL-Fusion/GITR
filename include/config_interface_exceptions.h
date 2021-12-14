#include <exception>

class unregistered_config_mapping: public std::exception
{
  public:

    unregistered_config_mapping( int const query_key )
      :
      query_key( query_key )
    { }

    int get_key() const { return query_key; }

    std::string static get_message() { return message; }

    const char* what() const throw()
    {
      return message.c_str();
    }

  private:

    int const query_key;

    std::string static const inline
    message{ "configuration interface error - no mapping exists for key: " };
};

class lookup_failed: public std::exception
{
  public:

    lookup_failed( std::string const query_key )
      :
      query_key( query_key )
    { }

    std::string get_key() const { return query_key; }

    std::string static get_message() { return message; }

    const char* what() const throw()
    {
      return message.c_str();
    }

  private:

    std::string const query_key;

    std::string static const inline
    message{ "configuration interface error - lookup failed for registered key: " };
};

class invalid_key: public std::exception
{
  public:

    invalid_key( std::string const &query_key )
      :
      query_key( query_key )
    { }

    /* Captain! All of these are identical. Template on int vs string for the unregistered
       lookups exception and then replace all this repeated code */
    std::string get_key() const { return query_key; }

    std::string static get_message() { return message; }

    const char* what() const throw() 
    {
      return message.c_str();
    }

  private:

    std::string const query_key;

    std::string static const inline 
    message{ "configuration interface error - key not found: " };
};

/*
class mismatched_dimensions : public std::exception
{
  public:

    mismatched_dimensions( std::string const &query_key )
      :
      query_key( query_key )
    { }

    std::string get_key() const { return query_key; }
};
*/

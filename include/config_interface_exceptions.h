#include <exception>

class unregistered_config_mapping: public std::exception
{
  public:

    unregistered_config_mapping( int const &query_key )
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

    /* Captain! Make this a static base class member */
    std::string static const inline
    message{ "configuration interface error - no mapping exists for key: " };
};

class lookup_failed: public std::exception
{
  public:

    lookup_failed( std::string const &query_key )
      :
      query_key( query_key )
    { }

    std::string get_key() { return query_key; }

    std::string static get_message() { return message; }

    const char* what() const throw()
    {
      return message.c_str();
    }

  private:

    std::string const &query_key;

    std::string static const inline
    message = "configuration interface error - lookup failed for registered key: ";
};

class invalid_key: public std::exception
{
  public:

    invalid_key( std::string const &query_key )
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

    std::string const &query_key;

    std::string static const inline 
    message{ "configuration interface error - key not found: " };
};

#include <exception>

class unregistered_config_mapping: public std::exception
{
  public:

    unregistered_config_mapping( int query_key )
      :
      query_key( query_key )
    { }

    const char* what() const throw()
    {
      std::string message = "configuration interface error - no mapping exists for key: ";
      message += " ";
      message += query_key;

      return message.c_str();
    }

  private:

    int query_key;
};

class lookup_failed: public std::exception
{
  public:

    lookup_failed( std::string query_key )
      :
      query_key( query_key )
    { }

    const char* what() const throw()
    {
      std::string message = "configuration interface error - lookup failed for existing key: ";
      message += " ";
      message += query_key;

      return message.c_str();
    }

  private:

    std::string query_key;
};

class invalid_key: std::exception
{
  public:

    invalid_key( const std::string &query_key )
      :
      query_key( query_key )
    { }

    const char* what() const noexcept override
    {
      std::cout << "configuration interface error - key not found: "; 

      return query_key.c_str();
    }

  private:

    std::string const &query_key;
};

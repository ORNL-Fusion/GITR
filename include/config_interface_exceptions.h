#include <exception>

class unregistered_config_mapping: public std::exception
{
  public:

    unregistered_config_mapping( int const query_key = -1 )
      :
      query_key( query_key )
    { }

    std::string get_key() const { return std::to_string( query_key ); }

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

class config_exception_base : public std::exception
{
public:
    explicit config_exception_base(const std::string& key)
        : query_key(key)
    { }

    std::string get_key() const { return query_key; }

    std::string get_message() const { return message; }

    const char* what() const noexcept override
    {
        return message.c_str();
    }

private:
    const std::string query_key;

    static const inline std::string message{ "configuration interface error - key not found: " };
};

/* Captain! Use config_base above to make a different message for each type of exception? */
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
class not_an_array : public std::exception
{
  public:

  not_a_scalar( std::string const &query_key )
    :
    query_key( query_key )
  { }

  private:
};

class not_a_scalar : public std::exception
{
  public:

  private:
};
*/

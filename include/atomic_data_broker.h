#include <vector>
class atomic_data_broker
{
  public:

  atomic_data_broker();

  std::vector< double > run_2();

  void run_1();

  void run();

  std::vector< double > values;

  std::vector< double > values2;
};

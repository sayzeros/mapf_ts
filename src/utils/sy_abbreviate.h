#include <vector>
#include <tuple>
#include <set>
#include <string>
#include <memory>

template<typename T>
using vector = std::vector<T>;

template<typename ... Args>
using tuple = std::tuple<Args...>;

using std::make_tuple;
using std::tie;

template<typename T>
using set = std::set<T>;

using std::string;
using std::to_string;

template<typename T>
using shared_ptr = std::shared_ptr<T>;

using std::make_shared;
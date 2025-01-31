//
// Created by oleksii on 31.01.25.
//

#ifndef QA2_HELPER_HPP
#define QA2_HELPER_HPP

#include <string>
#include <sstream>

namespace Helper {

template<typename T>
std::string to_string_with_precision(const T a_value, const int n) {
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

}
#endif //QA2_HELPER_HPP

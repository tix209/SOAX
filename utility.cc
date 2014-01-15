#include <iostream>
#include <fstream>
#include <sstream>
#include "utility.h"

namespace soax {

double String2Double(const std::string &s) {
  std::stringstream converter(s);
  double value;
  if (converter >> value)
    return value;
  std::cerr << "Conversion to double failed!" << std::endl;
  return 0.0;
}

unsigned String2Unsigned(const std::string &s) {
  std::stringstream converter(s);
  unsigned value;
  if (converter >> value)
    return value;
  std::cerr << "Conversion to unsigned failed!" << std::endl;
  return 0;
}

unsigned short String2UShort(const std::string &s) {
  std::stringstream converter(s);
  unsigned short value;
  if (converter >> value)
    return value;
  std::cerr << "Conversion to unsigned failed!" << std::endl;
  return 0;
}

std::string GetImagePath(const std::string &snake_path) {
  std::ifstream infile(snake_path.c_str());
  if (!infile) {
    std::cerr << "Couldn't open file: " << infile << std::endl;
    return "";
  }

  std::string line;
  getline(infile, line);
  std::stringstream buffer(line);
  std::string path, dump;
  buffer >> dump >> path;
  return path;
}

double Mean(const DataContainer &data) {
  if (data.empty()) {
    std::cerr << "No data. Mean value is set to 0.0." << std::endl;
    return 0.0;
  }

  double total = 0.0;
  for (unsigned i = 0; i < data.size(); ++i) {
    total += data[i];
  }
  return total / data.size();
}

double StandardDeviation(const DataContainer &data, double mean) {
  if (data.size() <= 1) {
    std::cerr << "StandardDeviation: data size is less than 2!"
              << std::endl;
    return 0.0;
  }
  double squared_sum = 0;
  for (unsigned i = 0; i < data.size(); ++i) {
    squared_sum += std::pow(data[i] - mean, 2);
  }
  squared_sum /= data.size()-1;
  return std::sqrt(squared_sum);
}

double Median(DataContainer &data) {
  if (data.empty()) {
    std::cerr << "No data. Median is set to 0.0." << std::endl;
    return 0.0;
  }
  size_t mid_index = data.size()/2;
  std::nth_element(data.begin(), data.begin() + mid_index, data.end());
  return data[mid_index];
}

double Minimum(const DataContainer &data) {
  if (data.empty()) {
    std::cerr << "No data. Minimum is set to 0.0." << std::endl;
    return 0.0;
  }

  double min_value = data.front();
  DataContainer::const_iterator it = ++data.begin();
  while (it != data.end()) {
    if (*it < min_value)
      min_value = *it;
    ++it;
  }
  return min_value;
}

double Maximum(const DataContainer &data) {
  if (data.empty()) {
    std::cerr << "No data. Maximum is set to 0.0." << std::endl;
    return 0.0;
  }

  double max_value = data.front();
  DataContainer::const_iterator it = ++data.begin();
  while (it != data.end()) {
    if (*it > max_value)
      max_value = *it;
    ++it;
  }
  return max_value;
}

std::string GetImageName(const std::string &snake_path) {
  std::ifstream infile(snake_path.c_str());
  if (!infile) {
    std::cerr << "Couldn't open file: " << infile << std::endl;
    return "";
  }

  std::string line;
  getline(infile, line);

  unsigned found_slash = line.find_last_of("/\\");
  return line.substr(found_slash + 1);
}

void PrintDataContainer(const DataContainer &data) {
  std::cout << "=====================" << std::endl;
  for (DataContainer::const_iterator it = data.begin();
       it != data.end(); ++it) {
    std::cout << *it << "\t";
  }
  std::cout << "\n====================" << std::endl;
}

} // namespace soax

#include <vector>
#include <fstream>
#include <sstream>
#include "droplet_info_dict.h"

namespace soax {

DropletInfoDict::DropletInfoDict() {}

bool DropletInfoDict::LoadInfo(const std::string &filename) {
  std::ifstream infile(filename.c_str());
  if (!infile) return false;

  std::string line;
  std::getline(infile, line); // discard the headline

  while (std::getline(infile, line)) {
    std::stringstream line_buffer(line);
    std::string cell;
    std::vector<std::string> record;
    while (std::getline(line_buffer, cell, ',')) {
      record.push_back(cell);
    }

    assert(record.size() == 6);
    DropletInfo info(std::stod(record[1]), std::stod(record[2]), std::stod(record[3]),
                     std::stod(record[4]), std::stoi(record[5]));
    assert(dict_.find(record.front()) == dict_.end());
    dict_[record.front()] = info;
  }

  return true;
}

}

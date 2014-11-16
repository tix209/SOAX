#ifndef SOAX_DROPLET_INFO_DICT_H_
#define SOAX_DROPLET_INFO_DICT_H_

#include <unordered_map>
#include "global.h"
#include "droplet_info.h"


namespace soax {

class DropletInfoDict {
 public:
  DropletInfoDict();
  bool LoadInfo(const std::string &info_filename);

  bool Has(const std::string &name) {
    return dict_.find(name) != dict_.end();
  }

  DropletInfo Get(const std::string &name) {return dict_[name];}

 private:
  std::unordered_map<std::string, DropletInfo> dict_;
};

} // namespace soax
#endif // SOAX_DROPLET_INFO_DICT_H_

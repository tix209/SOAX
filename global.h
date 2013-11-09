#ifndef SOAX_GLOBAL_H_
#define SOAX_GLOBAL_H_

// #include <limits>
// #include <deque>
#include "itkImage.h"

namespace soax {

// Constants
const unsigned kDimension = 3;
// const double kPi = 3.14159265358979323846264338;
// const double kEpsilon = std::numeric_limits<double>::epsilon();
// const double kPlusInfinity = std::numeric_limits<double>::max();
// const unsigned kBigNumber = std::numeric_limits<unsigned>::max();

// const unsigned kMinimumEvolvingSize = 5;


// typedef std::vector<double> DataContainer;
// typedef itk::Vector<double, kDimension> VectorType;
// typedef std::vector<VectorType>  VectorContainer;
// typedef itk::Point<double, kDimension> PointType;
// typedef std::deque<PointType>  PointContainer;
// typedef PointContainer::iterator PointIterator;
// typedef PointContainer::const_iterator PointConstIterator;
// typedef std::vector<Snake *> SnakeContainer;
typedef itk::Image<double, kDimension> ImageType;

// A macro to disallow the copy constructor and operator= functions
// This should be used in the private: declarations for a class
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)
} // namespace soax

#endif // SOAX_GLOBAL_H_

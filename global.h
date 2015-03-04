#ifndef SOAX_GLOBAL_H_
#define SOAX_GLOBAL_H_

#include <limits>
#include <deque>
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVersorTransform.h"

namespace soax {

// Constants
const unsigned kDimension = 3;
const double kPi = 3.14159265358979323846264338;
const double kEpsilon = std::numeric_limits<double>::epsilon();
const double kPlusInfinity = std::numeric_limits<double>::max();
const double kMinusInfinity = std::numeric_limits<double>::min();
const unsigned kBigNumber = std::numeric_limits<unsigned>::max();

const unsigned kMinimumEvolvingSize = 5;


typedef std::vector<double> DataContainer;
typedef itk::Vector<double, kDimension> VectorType;
typedef std::vector<VectorType>  VectorContainer;
typedef itk::Point<double, kDimension> PointType;
typedef std::deque<PointType>  PointContainer;
typedef std::set<PointType> PointSet;
typedef PointContainer::iterator PointIterator;
typedef PointContainer::const_iterator PointConstIterator;
typedef itk::Image<unsigned short, kDimension> ImageType;
typedef itk::Image<VectorType, kDimension> VectorImageType;

class Snake;
typedef std::vector<Snake *> SnakeContainer;
typedef std::set<Snake *> SnakeSet;
typedef SnakeContainer::iterator SnakeIterator;
typedef SnakeContainer::const_iterator SnakeConstIterator;

typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
// typedef itk::NearestNeighborInterpolateImageFunction<ImageType> InterpolatorType;
typedef itk::VectorLinearInterpolateImageFunction<VectorImageType>
VectorInterpolatorType;
typedef itk::VersorTransform<double> TransformType;


#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)
} // namespace soax

#endif // SOAX_GLOBAL_H_

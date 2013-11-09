#include "multisnake.h"
#include "itkImageFileReader.h"


namespace soax {

Multisnake::Multisnake() {}

Multisnake::~Multisnake() {}

void Multisnake::LoadImage(const std::string &filename) {
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  image_ = reader->GetOutput();
  try {
    reader->Update();
  } catch(itk::ExceptionObject &e) {
    std::cerr << "Exception caught when reading an image!" << std::endl;
    std::cerr << e << std::endl;
  }
  image_filename_ = filename;
  // const ImageType::SizeType &size = image_->GetLargestPossibleRegion().GetSize();
  // std::cout << "Image size: " << size << std::endl;
  // std::cout << image_filename_ << std::endl;
}

} // namespace soax

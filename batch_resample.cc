/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file implements resampling TIFF images to be isotropic (uniform
 * voxel size) in batch mode.
 *
 * Usage: ./batch_resample <input_dir> <output_dir> <xy_z_ratio>
 */

#include <iostream>
#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"


namespace fs = boost::filesystem;
typedef fs::path Path;

void ResampleImages(const Path &input_dir, const Path &output_dir,
                    double zspacing);
void ResampleImage(const std::string &input_filename,
                   const std::string &output_filename,
                   double zspacing);


int main (int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Batch Resample 3.7.0.\n"
		"Usage (for TIFF image only): ./batch_resample "
		"<input_dir> <output_dir> <z_spacing (relative to x/y)>" << std::endl;
    return -1;
  }

  Path input_dir(argv[1]);
  Path output_dir(argv[2]);
  double zspacing = atof(argv[3]);

  if (std::fabs(zspacing - 1) < 1e-6) {
    std::cout << "Images are isotropic. Quit." << std::endl;
    return 0;
  }

  if (fs::is_directory(input_dir)) {
    ResampleImages(input_dir, output_dir, zspacing);
  } else {
    std::cout << "Input dir should be a directory instead of a file. Abort."
              << std::endl;
  }
  return 0;
}

void ResampleImages(const Path &input_dir, const Path &output_dir,
                    double zspacing) {
  fs::directory_iterator end_it;
  for (fs::directory_iterator it(input_dir); it != end_it; ++it) {
    if (it->path().extension().string() == std::string(".tif")) {
      Path output_filename = output_dir;
      output_filename += it->path().stem();
      output_filename += "_iso.tif";
      ResampleImage(it->path().string(), output_filename.string(), zspacing);
      std::cout << "\"" << output_filename << "\" written." << std::endl;
    }
  }
}

void ResampleImage(const std::string &input_filename,
                   const std::string &output_filename,
                   double zspacing) {
  typedef itk::Image<unsigned short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(input_filename);
  try {
    reader->Update();
  } catch(itk::ExceptionObject &e) {
    std::cerr << "Exception caught when reading an image!" << std::endl;
    std::cerr << e << std::endl;
  }

  typedef itk::BSplineInterpolateImageFunction<ImageType, double, double>
      InterpolatorType;
  InterpolatorType::Pointer interp = InterpolatorType::New();
  interp->SetSplineOrder(3);
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInterpolator(interp);
  resampler->SetOutputOrigin(reader->GetOutput()->GetOrigin());
  resampler->SetOutputDirection(reader->GetOutput()->GetDirection());

  ImageType::SpacingType input_spacing;
  input_spacing.Fill(1.0);
  input_spacing[2] = zspacing;
  reader->GetOutput()->SetSpacing(input_spacing);

  const ImageType::SizeType &input_size =
      reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  ImageType::SpacingType output_spacing;
  output_spacing.Fill(1.0);
  ImageType::SizeType output_size;
  const unsigned int kDimension = 3;
  for (unsigned i = 0; i < kDimension; ++i) {
    output_size[i] = static_cast<ImageType::SizeValueType>(
        input_size[i] * input_spacing[i]);
  }

  resampler->SetOutputSpacing(output_spacing);
  resampler->SetSize(output_size);
  resampler->SetDefaultPixelValue(0.0);
  resampler->SetInput(reader->GetOutput());

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(output_filename);
  writer->SetInput(resampler->GetOutput());

  try {
    writer->Update();
  } catch( itk::ExceptionObject & exp ) {
    std::cerr << "Exception caught when write an image!" << std::endl;
    std::cerr << exp << std::endl;
  }
}

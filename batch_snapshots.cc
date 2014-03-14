/*
 * This file implements the batch generation of snapshots of a vtk
 * rendering window without showing it.
 */

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include "multisnake.h"
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageMapper3D.h>
#include <vtkImageActor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkActor.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageToVTKImageFilter.h>


typedef std::vector<vtkSmartPointer<vtkActor> > ActorsType;


std::string GetImageName(const std::string &name, bool suffix = false);
// vtkSmartPointer<vtkImageActor> ImportImage(const std::string &filename);
vtkSmartPointer<vtkActor> ImportSnakesAndJunctions(
    const std::string &filename, ActorsType &junction_actors);
void SaveScene(vtkSmartPointer<vtkImageActor> image_actor,
               vtkSmartPointer<vtkActor> snake_actor,
               const ActorsType &junction_actors,
               const std::string &filename);


int main(int argc, char **argv) {
  try {
    namespace po = boost::program_options;
    po::options_description generic("Generic options");
    generic.add_options()
        ("version,v", "Print version and exit")
        ("help,h", "Print help and exit");

    std::string snake_dir, output_dir;
    po::options_description required("Required options");
    required.add_options()
        ("image,i", po::value<std::string>()->required(),
         "Directory of image files")
        ("snake,s", po::value<std::string>(&snake_dir)->required(),
         "Directory of snake files")
        ("output,o", po::value<std::string>(&output_dir)->required(),
         "Output directory");
    po::options_description all("Allowed options");
    all.add(generic).add(required);
    po::variables_map vm;
    po::store(parse_command_line(argc, argv, all), vm);
    if (vm.count("version")) {
      const std::string version_msg(
          "Batch snapshots 1.0\n"
          "Copyright (C) 2014 Ting Xu, IDEA Lab, Lehigh University.");
      std::cout << version_msg << std::endl;
      return EXIT_SUCCESS;
    }

    if (vm.count("help")) {
      std::cout << "Batch snapshots usage: \n" << all;
      return EXIT_SUCCESS;
    }
    po::notify(vm);

    namespace fs = boost::filesystem;
    fs::path image_path(vm["image"].as<std::string>());
    if (!fs::exists(image_path)) {
      std::cerr << image_path << " does not exist. Abort." << std::endl;
      return EXIT_FAILURE;
    }

    if (snake_dir.at(snake_dir.length()-1) != '/')  snake_dir += '/';
    if (output_dir.at(output_dir.length()-1) != '/')  output_dir += '/';

    typedef std::vector<fs::path> Paths;
    Paths image_paths;
    std::copy(fs::directory_iterator(image_path), fs::directory_iterator(),
              back_inserter(image_paths));
    std::sort(image_paths.begin(), image_paths.end());

    typedef itk::Image<unsigned short, 2> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();

    typedef itk::RescaleIntensityImageFilter<ImageType,
                                             ImageType> RescalerType;
    RescalerType::Pointer rescaler = RescalerType::New();
    rescaler->SetInput(reader->GetOutput());
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);

    typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
    ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(rescaler->GetOutput());

    for (Paths::const_iterator image_it(image_paths.begin());
         image_it != image_paths.end(); ++image_it) {
      // fs::directory_iterator image_end_it;
      // for (fs::directory_iterator image_it(image_path);
      //      image_it != image_end_it; ++image_it) {

      // vtkSmartPointer<vtkImageActor> image_actor = ImportImage(
      //     image_it->string());

      reader->SetFileName(image_it->string());
      reader->UpdateLargestPossibleRegion();
      rescaler->Update();
      connector->Update();

      vtkSmartPointer<vtkImageActor> image_actor =
          vtkSmartPointer<vtkImageActor>::New();
      image_actor->GetMapper()->SetInputData(connector->GetOutput());

      std::string image_name = GetImageName(image_it->string());
      std::string snake_filename = snake_dir + image_name + ".txt";
      std::cout << snake_filename << std::endl;

      ActorsType junction_actors;
      vtkSmartPointer<vtkActor> snake_actor =
          ImportSnakesAndJunctions(snake_filename, junction_actors);
      // std::cout << image_it->filename() << std::endl;
      std::string output_name = output_dir + image_name + ".png";
      SaveScene(image_actor, snake_actor, junction_actors, output_name);
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}

/*
 * Function: GetImageName
 * Usage: image_name = GetImageName(image_path);
 * ---------------------------------------------
 * Extract the image name from image path without the image format
 * suffix. For example, given "/path/to/my/image.png", this function
 * return "image".
 */
std::string GetImageName(const std::string &name, bool suffix) {
  unsigned last_slash_pos = name.find_last_of("/\\");
  if (!suffix) {
    unsigned last_dot_pos = name.find_last_of(".");
    return name.substr(last_slash_pos+1, last_dot_pos-last_slash_pos-1);
  } else {
    return name.substr(last_slash_pos+1);
  }
}

/*
 * Function: ImportImage
 * Usage: image_actor_smart_ptr = ImportImage(image_filename);
 * -----------------------------------------------------------
 * Import an image as a vtkImageActor.
 */
// vtkSmartPointer<vtkImageActor> ImportImage(const std::string &filename) {
//   typedef itk::Image<unsigned short, 2> ImageType;
//   typedef itk::ImageFileReader<ImageType> ReaderType;
//   ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName(filename);
//   reader->Update();

//   typedef itk::RescaleIntensityImageFilter<ImageType,
//                                            ImageType> RescalerType;
//   RescalerType::Pointer rescaler = RescalerType::New();
//   rescaler->SetInput(reader->GetOutput());
//   rescaler->SetOutputMinimum(0);
//   rescaler->SetOutputMaximum(255);
//   rescaler->Update();

//   typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
//   ConnectorType::Pointer connector = ConnectorType::New();
//   connector->SetInput(rescaler->GetOutput());
//   connector->Update();

//   vtkSmartPointer<vtkImageActor> image_actor =
//       vtkSmartPointer<vtkImageActor>::New();
//   vtkImageData *image_data = NULL;
//   image_data->ShallowCopy(connector->GetOutput());
//   // image_actor->GetMapper()->SetInputData(connector->GetOutput());
//   image_actor->GetMapper()->SetInputData(image_data);

//   return image_actor;
// }


/*
 * Function: ImportSnakesAndJunctions
 *
 * Usage: image_actor_smart_ptr =
 * ImportSnakesAndJunctions(image_filename, junction_actors);
 *
 * Import snakes from file as a single vtkActor and construct actors
 * for each junction locations and put them into the junction_actors
 * container.
 */
vtkSmartPointer<vtkActor> ImportSnakesAndJunctions(
    const std::string &filename, ActorsType &junction_actors) {
  return NULL;
}


/*
 * Function: SaveScene
 * Usage: SaveScene(image_actor, snake_actor, output_image_filename);
 * ------------------------------------------------------------------
 * Save the off-screen rendered scene as a png image.
 */
void SaveScene(vtkSmartPointer<vtkImageActor> image_actor,
               vtkSmartPointer<vtkActor> snake_actor,
               const ActorsType &junction_actors,
               const std::string &filename) {
  vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
  if (image_actor) {
    // image_actor->Print(std::cout);
    renderer->AddActor(image_actor);
  }
  if (snake_actor)
    renderer->AddActor(snake_actor);

  for (ActorsType::const_iterator it = junction_actors.begin();
       it != junction_actors.end(); ++it) {
    renderer->AddActor(*it);
  }
  renderer->ResetCamera();

  vtkSmartPointer<vtkRenderWindow> window =
      vtkSmartPointer<vtkRenderWindow>::New();
  window->OffScreenRenderingOn(); // important for removing
  // misplaced image patches!!
  window->AddRenderer(renderer);
  // window->SetAlphaBitPlanes(1); //enable usage of alpha channel
  window->Render();

  vtkSmartPointer<vtkWindowToImageFilter> w2i =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
  w2i->SetInput(window);
  w2i->SetMagnification(4);
  // w2i->SetInputBufferTypeToRGBA();
  w2i->Update();

  vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetInputConnection(w2i->GetOutputPort());
  writer->SetFileName(filename.c_str());
  writer->Write();
}

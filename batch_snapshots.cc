/*
 * File: batch_snapshots.cc
 * ~~~~~~~~~~~~~~~~~~~~~~~~
 * This file implements the batch generation of snapshots of a vtk
 * rendering window without showing it.
 *
 */

#include "vtkSmartPointer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkImageMapper3D.h"
#include "vtkImageActor.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkActor.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkProperty.h"
#include "vtkSphereSource.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageToVTKImageFilter.h"

#include "boost/program_options.hpp"
#include "boost/filesystem.hpp"
#include <iostream>
#include <fstream>
#include <sstream>



typedef itk::Image<unsigned short, 2> ImageType;
typedef std::vector<vtkSmartPointer<vtkPolyData> > PolyDataContainer;


std::string GetImageName(const std::string &name, bool suffix = false);

ImageType::Pointer ReadImage(const std::string &filename);

vtkSmartPointer<vtkPolyData> ImportSnakesAndJunctions(
    const std::string &filename, PolyDataContainer &junctions);

void SaveScene(vtkImageData *image, vtkSmartPointer<vtkPolyData> snakes,
               const PolyDataContainer &junctions,
               const std::string &filename);

vtkSmartPointer<vtkActor> SetupSnakeActor(
    vtkSmartPointer<vtkPolyData> snakes);
vtkSmartPointer<vtkActor> SetupJunction(
    vtkSmartPointer<vtkPolyData> junction);




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
      const std::string version_msg("Batch snapshots 1.0\n"
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

    for (Paths::const_iterator image_it(image_paths.begin());
         image_it != image_paths.end(); ++image_it) {
      ImageType::Pointer itk_image = ReadImage(image_it->string());

      // connecting to VTK
      typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
      ConnectorType::Pointer connector = ConnectorType::New();
      connector->SetInput(itk_image);
      connector->Update();
      // vtkSmartPointer<vtkImageData> vtk_image;
      // vtk_image.TakeReference(connector->GetOutput());
      vtkImageData *vtk_image = connector->GetOutput();

      // Get snake file path and output image path
      std::string image_filename = GetImageName(image_it->string());
      std::string snake_filename = snake_dir + image_filename + ".txt";
      std::string output_filename = output_dir + image_filename + ".png";
      std::cout << output_filename << std::endl;

      PolyDataContainer junctions;
      vtkSmartPointer<vtkPolyData> snakes = ImportSnakesAndJunctions(
          snake_filename, junctions);

      SaveScene(vtk_image, snakes, junctions, output_filename);
      // vtk_image->Delete(); // connector owns it, no need to delete
    }
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/*
 * Function: GetImageName
 * Usage: image_filename = GetImageName(image_path);
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
 * Function: ReadImage
 * Usage: itk_image_smart_ptr = ReadImage(image_filename);
 * -----------------------------------------------------------
 * Import an image as a itkImage and rescale the intensity to 0~255.
 */
ImageType::Pointer ReadImage(const std::string &filename) {
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(filename);
  // reader->Update();

  typedef itk::RescaleIntensityImageFilter<ImageType,
                                           ImageType> RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput(reader->GetOutput());
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->Update();
  return rescaler->GetOutput();
}


/*
 * Function: ImportSnakesAndJunctions
 * Usage: polydata_ptr = ImportSnakesAndJunctions(filename, junctions);
 * --------------------------------------------------------------------
 * Import snakes from file as a single vtkPolyData and junction coordinates
 * as many vtkPolyData and put them into junctions.
 */
vtkSmartPointer<vtkPolyData> ImportSnakesAndJunctions(
    const std::string &filename, PolyDataContainer &junctions) {
  std::ifstream infile(filename.c_str());
  if (!infile.is_open()) {
    std::cerr << "Couldn't open snake file: " << filename << std::endl;
    return NULL;
  }

  vtkSmartPointer<vtkPolyData> curve = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

  bool open = true;
  std::string line;
  unsigned index = 0;
  vtkIdType cell_index[2];
  vtkFloatingPointType coordinates[2];
  unsigned start = 0;
  while(std::getline(infile, line)) {
    if (isalpha(line[0])) {
      continue;
    } else if (line[0] == '#') { // previous is a last point
      if (!open && index) {
        cell_index[0] = index-1;
        cell_index[1] = start;
        cells->InsertNextCell(2, cell_index);
      }
      open = (line[1] != '0');
    } else if (line[0] == '[') {
      std::istringstream buffer(line);
      char padding;
      double x, y, z;
      buffer >> padding >> x >> padding >> y >> padding >> z >> padding;
      vtkSmartPointer<vtkSphereSource> source =
          vtkSmartPointer<vtkSphereSource>::New();
      source->SetCenter(x, y, z);
      source->SetRadius(3.0);
      source->Update();
      junctions.push_back(source->GetOutput());
    } else {
      std::istringstream stream(line);
      int snake_index, point_index;
      double x, y, z;
      stream >> snake_index >> point_index >> x >> y >> z;

      coordinates[0] = x;
      coordinates[1] = y;
      coordinates[2] = z;
      points->InsertPoint(index, coordinates);
      if (point_index) {
        cell_index[0] = index-1;
        cell_index[1] = index;
        cells->InsertNextCell(2, cell_index);
      } else {
        start = index;
      }
      index++;
    }
  }

  curve->SetPoints(points);
  curve->SetLines(cells);
  return curve;
}


/*
 * Function: SaveScene
 * Usage: SaveScene(image, snakes, junctions, output_filename);
 * ------------------------------------------------------------------
 * Save the off-screen rendered scene as a png image.
 */
void SaveScene(vtkImageData *image, vtkSmartPointer<vtkPolyData> snakes,
               const PolyDataContainer &junctions,
               const std::string &filename) {
  vtkSmartPointer<vtkImageActor> image_actor =
      vtkSmartPointer<vtkImageActor>::New();
  image_actor->GetMapper()->SetInputData(image);

  vtkSmartPointer<vtkActor> snake_actor = SetupSnakeActor(snakes);

  vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(image_actor);
  if (snake_actor)
    renderer->AddActor(snake_actor);

  for (PolyDataContainer::const_iterator it = junctions.begin();
       it != junctions.end(); ++it) {
    vtkSmartPointer<vtkActor> junction = SetupJunction(*it);
    renderer->AddActor(junction);
  }

  vtkSmartPointer<vtkRenderWindow> window =
      vtkSmartPointer<vtkRenderWindow>::New();
  window->OffScreenRenderingOn(); // important for removing misplaced
                                  // image patches!!
  window->AddRenderer(renderer);
  // window->SetAlphaBitPlanes(1); //enable usage of alpha channel
  window->Render();

  vtkSmartPointer<vtkWindowToImageFilter> w2i =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
  w2i->Modified();
  w2i->SetInput(window);
  w2i->SetMagnification(8);
  // w2i->SetInputBufferTypeToRGBA();
  w2i->Update();

  vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetInputConnection(w2i->GetOutputPort());
  writer->SetFileName(filename.c_str());
  writer->Write();

  // vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
  //     vtkSmartPointer<vtkRenderWindowInteractor>::New();
  // renderWindowInteractor->SetRenderWindow(window);
  // renderer->ResetCamera();
  // window->Render();
  // renderWindowInteractor->Start();
}


/*
 * Function: SetupSnakeActor
 * Usage: snake_actor_ptr = SetupSnakeActor(snake_polydata)
 *---------------------------------------------------------
 * Turn the snake polydata into a snake actor.
 */
vtkSmartPointer<vtkActor> SetupSnakeActor(
    vtkSmartPointer<vtkPolyData> snakes) {
  if (!snakes) return NULL;
  vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(snakes);
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  actor->GetProperty()->SetInterpolationToPhong();
  actor->GetProperty()->SetOpacity(1.0);
  actor->GetProperty()->SetAmbient(0.2);
  actor->GetProperty()->SetDiffuse(0.7);
  actor->GetProperty()->SetSpecular(0.6);
  actor->GetProperty()->SetSpecularPower(50);
  actor->GetProperty()->SetColor(1.0, 0.0, 1.0);
  actor->GetProperty()->SetLineWidth(2.0);

  return actor;
}

/*
 * Function: SetupJunction
 * Usage: junction_actor = SetupJunction(junction_polydata)
 * -------------------------------------------------------
 * Turn a junction polydata into a sphere actor.
 */
vtkSmartPointer<vtkActor> SetupJunction(
    vtkSmartPointer<vtkPolyData> junction) {
  vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(junction);
  vtkSmartPointer<vtkActor> sphere = vtkSmartPointer<vtkActor>::New();
  sphere->SetMapper(mapper);
  sphere->GetProperty()->SetColor(0, 1, 0);
  return sphere;
}

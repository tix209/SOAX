#include <QApplication>
#include "boost/program_options.hpp"
#include "main_window.h"
#include "multisnake.h"


int main(int argc, char **argv) {
  // GUI mode
  if (argc == 1) {
    QApplication app(argc, argv);
    soax::MainWindow window;
    window.resize(800, 600);
    window.show();
    return app.exec();
  }

  // Commandline mode
  namespace po = boost::program_options;
  po::options_description generic("Generic options");
  generic.add_options()
      ("version,v", "Print version and exit")
      ("help,h", "Print help and exit");

  std::string image_path, parameter_path, snake_path;
  po::options_description required("Required options");
  required.add_options()
      ("image,i", po::value<std::string>(&image_path)->required(),
       "Image file path")
      ("parameters,p", po::value<std::string>(&parameter_path)->
       required(), "Parameter file path")
      ("snake,s", po::value<std::string>(&snake_path)->required(),
       "Snake file path");

  double grad_diff(0.0);
  double stretch(0.0);
  po::options_description optional("Optional options");
  optional.add_options()
      ("grad-diff", po::value<double>(&grad_diff),
       "Ridge threshold for SOAC initialization")
      ("stretch", po::value<double>(&stretch),
       "Magnitude of stretching force");

  po::options_description all("Allowed options");
  all.add(generic).add(required).add(optional);
  po::variables_map vm;
  po::store(parse_command_line(argc, argv, all), vm);

  if (vm.count("version")) {
    const std::string version_msg(
        "SOAX 3.5.0\n"
        "Copyright (C) 2013 Ting Xu, IDEA Lab, Lehigh University.");
    std::cout << version_msg << std::endl;
    return EXIT_SUCCESS;
  }

  if (vm.count("help")) {
    std::cout << "Usage: \n" << all;
    return EXIT_SUCCESS;
  }

  po::notify(vm);
  soax::Multisnake multisnake;
  multisnake.LoadImage(image_path);
  multisnake.LoadParameters(parameter_path);

  if (vm.count("grad-diff")) {
    std::cout << "grad-diff is set to " << grad_diff << std::endl;
    multisnake.set_ridge_threshold(grad_diff);
  }

  if (vm.count("stretch")) {
    std::cout << "stretch is set to " << stretch << std::endl;
    multisnake.set_stretch_factor(stretch);
  }

  multisnake.PrintParameters();
  multisnake.ScaleImageIntensity();
  multisnake.ComputeImageGradient();
  multisnake.InitializeSnakes();
  multisnake.SortSnakesOnLength(multisnake.initial_snakes());
  multisnake.DeformSnakes();
  multisnake.CutSnakesAtTJunctions();
  multisnake.GroupSnakes();
  multisnake.SaveSnakes(multisnake.converged_snakes(), snake_path);
  std::cout << "Segmentation completed.\n" << std::endl;
  return EXIT_SUCCESS;
}

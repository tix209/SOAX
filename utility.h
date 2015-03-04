/**
 * Copyright (c) 2015, Lehigh University
 * All rights reserved.
 * See COPYING for license.
 *
 * This file defines utility functions for SOAX.
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include <string>
#include "./global.h"

namespace soax {

double String2Double(const std::string &s);
unsigned String2Unsigned(const std::string &s);
unsigned short String2UShort(const std::string &s);
std::string GetImagePath(const std::string &snake_path);
double Mean(const DataContainer &data);
double StandardDeviation(const DataContainer &data, double mean);
double Median(DataContainer &data);
double Minimum(const DataContainer &data);
double Maximum(const DataContainer &data);
/*
 * Get the image filename from snake file. Assuming the path is at
 * the first line in a snake file.
 */
std::string GetImageName(const std::string &snake_path);

void PrintDataContainer(const DataContainer &data);
}  // namespace soax

#endif  // UTILITY_H_

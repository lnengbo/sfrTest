﻿#pragma once
#include "opencv2/opencv.hpp"
#include <vector>
#define M_PI 3.1415

int SFRCalculation(cv::Mat &ROI, double gamma,std::vector<double> &resVec,std::vector<double>&,std::vector<double>&);

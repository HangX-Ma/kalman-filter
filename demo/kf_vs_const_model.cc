#include <iostream>
#include <string>
#include <vector>
#include <random>

#include "visual_servo.h"
#include "kalman_filter.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

constexpr const uint32_t kScreenWidth = 2400;
constexpr const uint32_t kScreenHeight = 1344;
/**
 * Camera parameters
 *  - 16:9 S2
 * 	- CMOS size 36.0x24.0 mm
 * 
 * @ref https://cam.start.canon/zh/C013/manual/html/UG-10_Reference_0100.html
 */
constexpr const double kPlaneWidth = 0.036;	 // 36 mm
constexpr const double kPlaneHeight = 0.024; // 24 mm

int main (int argc,char *argv[]) {
    (void)argc;
    (void)argv;


    return 0;
}

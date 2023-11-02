#pragma once

#include <math.h>

namespace ctrl
{

static double ptzFocalLength[] = {
    // TODO: Add focal length value
};
typedef struct PtzProperty {
    int centerPixel;
    double imgPlaneSize;
    double maxViewAngle;
    double maxSpeed;

    double motorSpeed;
    double pixelSpeed;

    unsigned char getNormalizedMotorSpeed(double pixel, double deltaPixel, double period /*sec*/) {
        return (unsigned char) (getMotorSpeed(pixel, deltaPixel, period) / maxSpeed * 255.0);
    }

private:
    /**
     * @brief Get motor speed according to the pixel change.
     * @warning This function has side effect which will update the private 'motorSpeed' variable.
     * @note Rotation speed is always positive
     */
    double getMotorSpeed(double pixel, double deltaPixel, double period /*sec*/) {
        motorSpeed = (getAngle(pixel + deltaPixel) - getAngle(pixel)) / period;
        return abs(motorSpeed) > maxSpeed ? maxSpeed : abs(motorSpeed);
    }

public:
    double getAngle(int pixel) {
        return atan2(pixelToPlane(pixel) * maxViewAngle, (double) centerPixel);
    }

    // convert pixel unit to img plane unit
    double pixelToPlane(int pixel) {
        return imgPlaneSize / (double) centerPixel * 0.5 * (double) pixel;
    }

    // convert img plane unit to pixel unit
    int planeToPixel(double plane) {
        double res = (double) centerPixel * 2.0 / imgPlaneSize * plane;
        return res > 0.0 ? floor(res) : ceil(res);
    }

} PtzProperty_t;

typedef struct PtzParam {
    PtzProperty_t pan; // horizontal
    PtzProperty_t tilt; // vertical

    unsigned int focalSel;
    double ctrlPeriod;
    unsigned int lostCount;

    /** Update focal length property. This will also update the frequently
     * used variable 'unitPan' and 'unitTilt' for speed conversion */
    void updateFocal(unsigned int sel) {
        focalSel = sel;
        unitPan = ptzFocalLength[focalSel] / pan.pixelToPlane(1);
        unitTilt = ptzFocalLength[focalSel] / tilt.pixelToPlane(1);
    }

    double getPtzPixelSpeedPan(double pixelPan, double pixelTilt) {
        double deltaPixelPan = (pixelPan - pan.centerPixel);
        double deltaPixelTilt = (pixelTilt - tilt.centerPixel);
        return pan.motorSpeed * (unitPan + (double) deltaPixelPan * deltaPixelPan / unitPan)
                    - tilt.motorSpeed / unitTilt * (double) deltaPixelPan * deltaPixelTilt;
    }

    double getPtzPixelSpeedTilt(double pixelPan, double pixelTilt) {
        double deltaPixelPan = (pixelPan - pan.centerPixel);
        double deltaPixelTilt = (pixelTilt - tilt.centerPixel);
        return -tilt.motorSpeed * (unitTilt + (double) deltaPixelTilt * deltaPixelTilt / unitTilt)
                + pan.motorSpeed * (double) deltaPixelPan * deltaPixelTilt / unitPan ;
    }

private:
    double unitPan;
    double unitTilt;

} PtzParam_t;

} //ctrl
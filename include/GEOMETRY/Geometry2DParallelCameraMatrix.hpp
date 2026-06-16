#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "MATRIX/Matrix.hpp"
#include "MATRIX/utils.hpp"
#include "PROG/KCTException.hpp"

// Interfaces
#include "GEOMETRY/Geometry2DParallelCameraMatrixI.hpp"

namespace KCT::geometry {

using namespace KCT::matrix;

/**
 *Class to represent projection matrices of parallel rays that are 1x4 matrices.
 */
class Geometry2DParallelCameraMatrix : public Geometry2DParallelCameraMatrixI
{
public:
    Geometry2DParallelCameraMatrix(const double* vector4);
    Geometry2DParallelCameraMatrix(std::initializer_list<double> vector4_copy);
    /**Constructor from the 1x4 Matrix.*/
    Geometry2DParallelCameraMatrix(const Matrix& pm);

    // Initialization from ray direction, detector origin and detector direction
    Geometry2DParallelCameraMatrix(const std::array<double, 3> rayDirection,
                                   const std::array<double, 3> detectorOrigin,
                                   const std::array<double, 3> VX);

    bool operator==(const Geometry2DParallelCameraMatrix& rhs) const;

    // Implements Geometry2DParallelCameraMatrixI
    void directionVectorVR(double* vector3) const;
    void directionVectorVX(double* vector3) const;
    std::array<double, 3> directionVectorVR() const;
    std::array<double, 3> directionVectorVX() const;

    double pixelSpacing() const;

    void backprojectToPosition(const double PX, double* vector3) const;
    std::array<double, 3> backprojectToPosition(const double PX) const;

    void project(const double x0, const double x1, const double x2, double* PX) const;
    void project(const float x0, const float x1, const float x2, float* PX) const;

    void projectionMatrixAsVector4(double* vector4) const;
    std::array<double, 4> projectionMatrixAsVector4() const;

    virtual void projectionMatrixPXAsVector3(double* vector3, double z = 0.0) const;
    virtual std::array<double, 3> projectionMatrixPXAsVector3(double z = 0.0) const;

    std::string toString(std::string name = "P") const;

    /**
     * Auxiliary helper function
     *
     * It is analogous to ASTRA-like parallel initialization in one detector dimension.
     *
     * @param detector_spacing_x X-distance between two adjacent detector pixels
     * @param projection_size_x Count of detector pixels
     * @param angle Angle in radians between positive X axis and the direction vector VR
     *
     * @return Geometry2DParallelCameraMatrix object with given parameters
     */
    static Geometry2DParallelCameraMatrix initializeFromParameters(const double detector_spacing_x,
                                                                   const uint32_t projection_size_x,
                                                                   const double angle);

    /**
     * Auxiliary helper function
     *
     * Initialization from explicit geometry description.
     *
     * @param projection_size_x Count of detector pixels
     * @param rayDirection Direction of incoming rays
     * @param detectorCenter Center of the detector
     * @param VX Direction and spacing vector of detector coordinate PX
     *
     * @return Geometry2DParallelCameraMatrix object with given parameters
     */
    static Geometry2DParallelCameraMatrix
    initializeFromParameters(const uint32_t projection_size_x,
                             const std::array<double, 3> rayDirection,
                             const std::array<double, 3> detectorCenter,
                             const std::array<double, 3> VX);

    /**
     * Auxiliary helper function
     *
     * Intended to distribute directions of incoming rays around [0,pi) for given number of angles.
     *
     * @param detector_spacing_x X-distance between two adjacent detector pixels
     * @param projection_size_x Count of detector pixels
     * @param anglesCount Number of angles in the interval [0,pi)
     * @param angleNum Angle ID
     *
     * @return Geometry2DParallelCameraMatrix object with given parameters
     */
    static Geometry2DParallelCameraMatrix initializeFromParameters(const double detector_spacing_x,
                                                                   const uint32_t projection_size_x,
                                                                   const uint32_t anglesCount,
                                                                   const uint32_t angleNum);

private:
    std::array<double, 4> projectionMatrixVector; // Projection matrix
    std::array<double, 3> rayDirectionVector; // Unit direction of incoming rays
    static constexpr double zeroPrecisionTolerance = 1e-10;
    static constexpr double pi = 3.141592653589793238462643383279502884;
};

} // namespace KCT::geometry

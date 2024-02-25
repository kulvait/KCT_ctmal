#pragma once
// Logging
#include <plog/Log.h>

// Standard libraries
#include <array>
#include <iomanip>
#include <iostream>

// Internal libraries
#include "MATRIX/LUDoolittleForm.hpp"
#include "MATRIX/Matrix.hpp"
#include "MATRIX/SquareMatrix.hpp"
#include "MATRIX/utils.hpp"
#include "PROG/KCTException.hpp"

// Interfaces
#include "GEOMETRY/Geometry3DParallelCameraMatrixI.hpp"

namespace KCT::geometry {

using namespace KCT::matrix;

/**
 *Class to represent projection matrices of parallel rays that are 2x4 matrices.
 */
class Geometry3DParallelCameraMatrix : public Geometry3DParallelCameraMatrixI
{
public:
    Geometry3DParallelCameraMatrix(const double* vector8);
    Geometry3DParallelCameraMatrix(std::initializer_list<double> vector8_copy);
    /**Constructor from the 2x4 Matrix.*/
    Geometry3DParallelCameraMatrix(const Matrix& pm);
    // ASTRA like initialization
    Geometry3DParallelCameraMatrix(const std::array<double, 3> rayDirection,
                                   const std::array<double, 3> detectorOrigin,
                                   const std::array<double, 3> VX,
                                   const std::array<double, 3> VY);
    bool operator==(const Geometry3DParallelCameraMatrix& rhs);

    // Implements Geometry3DParallelCameraMatrixI
    double pixelSkew() const;
    double pixelArea() const;

    void directionVectorVR(double* vector3) const;
    void directionVectorVX(double* vector3) const;
    void directionVectorVY(double* vector3) const;
    /**
     *
     * @return Unit vector in the direction of incomming rays or its oposite.
     */
    std::array<double, 3> directionVectorVR() const;
    std::array<double, 3> directionVectorVX() const;
    std::array<double, 3> directionVectorVY() const;
    void backprojectToPosition(const double pi, const double pj, double* vector3) const;
    std::array<double, 3> backprojectToPosition(const double pi, const double pj) const;
    void project(const double x0, const double x1, const double x2, double* pi, double* pj) const;
    void project(const float x0, const float x1, const float x2, float* pi, float* pj) const;

    void projectionMatrixAsVector8(double* vector8) const;
    std::array<double, 8> projectionMatrixAsVector8() const;
    void projectionMatrixPXAsVector4(double* vector4) const;
    std::array<double, 4> projectionMatrixPXAsVector4() const;
    void projectionMatrixPXAsVector3(double* vector3, double z = 0.0) const;
    std::array<double, 3> projectionMatrixPXAsVector3(double z = 0.0) const;

    std::string toString(std::string name = "P") const;

    /**
     * Auxiliary helper function
     *
     * It is analogous to ASTRA parallel3d initialization, see
     * https://www.astra-toolbox.com/docs/geom3d.html#projection-geometries
     * I am not sure if ASTRA obeys the same convention for angle meaning, needs to be tested.
     *
     * @param detector_spacing_x X-distance between two adjacent pixels on detector
     * @param detector_spacing_y Y-distance between two adjacent pixels on detector
     * @param projection_size_x Count of X detector pixels, col count
     * @param projection_size_y Count of Y detector pixels, row count
     * @param angle Angle in radians between positive X axis and the direction vector VR
     * computed ccw to the X axis, consistent with Astra toolbox implementation, see
     * https://github.com/astra-toolbox/astra-toolbox/blob/fa2ec619edb2994e828897e80c06e7fb35c55c44/src/ParallelProjectionGeometry3D.cpp#L178
     *
     * @return Geometry3DParallelCameraMatrix object with given parameters
     */
    static Geometry3DParallelCameraMatrix initializeFromParameters(const double detector_spacing_x,
                                                                   const double detector_spacing_y,
                                                                   const uint32_t projection_size_x,
                                                                   const uint32_t projection_size_y,
                                                                   const double angle);
    /**
     * Auxiliary helper function
     *
     * It is analogous to ASTRA parallel3d_vec initialization, see
     * https://www.astra-toolbox.com/docs/geom3d.html#projection-geometries
     *
     * @param projection_size_x Count of X detector pixels
     * @param projection_size_y Count of Y detector pixels
     * @param rayDirection See
     * https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html
     * @param detectorCenter See
     * https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html
     * @param VX See
     * https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html
     * @param VY See
     * https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html
     *
     * @return Geometry3DParallelCameraMatrix object with given parameters
     */
    static Geometry3DParallelCameraMatrix
    initializeFromParameters(const uint32_t projection_size_x,
                             const uint32_t projection_size_y,
                             const std::array<double, 3> rayDirection,
                             const std::array<double, 3> detectorCenter,
                             const std::array<double, 3> VX,
                             const std::array<double, 3> VY);

    /**
     *
     * Auxiliary helper function
     *
     * Intended to distribute directions of incomming rays around [0,pi) for given number of angles
     *
     * @param detector_spacing_x X-distance between two adjacent pixels on detector
     * @param detector_spacing_y Y-distance between two adjacent pixels on detector
     * @param projection_size_x Count of X detector pixels
     * @param projection_size_y Count of Y detector pixels
     * @param anglesCount Number of anges in the interval [0,pi)
     * @param angleNum Angle ID.
     *
     * @return Geometry3DParallelCameraMatrix object with given parameters
     */
    static Geometry3DParallelCameraMatrix initializeFromParameters(const double detector_spacing_x,
                                                                   const double detector_spacing_y,
                                                                   const uint32_t projection_size_x,
                                                                   const uint32_t projection_size_y,
                                                                   const uint32_t anglesCount,
                                                                   const uint32_t angleNum);

private:
    std::array<double, 8> projectionMatrixVector; // Projection matrix
    static constexpr double zeroPrecisionTolerance = 1e-10;
    static constexpr double pi = 3.141592653589793238462643383279502884;
};
} // namespace KCT::geometry

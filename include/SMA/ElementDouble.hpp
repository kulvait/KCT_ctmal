#pragma once

#include "littleEndianAlignment.h"
#include "plog/Log.h"

#include <string>

namespace CTL::matrix {
/**
 * Buffered sparse matrix is a structure
 *
 */
class ElementDouble
{
public:
    /**
     */
    ElementDouble(uint8_t* element)
    {
        i = util::nextUint32(element);
        j = util::nextUint32(&element[4]);
        v = util::nextDouble(&element[8]);
    }

    ElementDouble(uint32_t i, uint32_t j, double v)
    {
        this->i = i;
        this->j = j;
        this->v = v;
    }

    /** Position is in the units of the file
     *
     */
    uint32_t getI() const { return i; }

    uint32_t getJ() const { return j; }

    double getVal() const { return v; }

    void addToVal(double x) { v = v + x; }

    // Voxel first sort
    bool operator<(const ElementDouble& e) const
    {
        if(i < e.i)
        {
            return true;
        }
        if(i == e.i)
        {
            return (j < e.j);
        }
        return false;
    }

    // Pixel first sort
    bool operator>(const ElementDouble& e) const
    {
        if(j < e.j)
        {
            return true;
        }
        if(j == e.j)
        {
            return (i < e.i);
        }
        return false;
    }

    uint32_t i, j;
    double v;
};

} // namespace CTL::matrix

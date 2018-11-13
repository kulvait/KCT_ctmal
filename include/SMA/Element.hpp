#pragma once

#include "plog/Log.h"
#include "littleEndianAlignment.h"

#include <string>

namespace CTL::matrix {
/**
 * Buffered sparse matrix is a structure
 *
 */
class Element
{
public:
    /** 
     */
    Element(uint8_t * elementBytes)
    {
	i = 	util::nextUint32(element);
	j = 
		util::nextUint32(&element[4]);
		v=util::nextDouble(&element[8]);
    }

    Element(uint32_t i, uint32_t j, double v)
    {
	this->i = i;
	this->j = j;
	this->v = v;
    }

/** Position is in the units of the file
*
*/
	uint32_t i()
	{
		return i;
	}

uint32_t j()
{
return j;
}

double val()
{
return v;
} 

double addToVal(double x)
{
	v = v + x;
}

protected:
	uint32_t i,j;
	double v;
};

} // namespace CTL::matrix

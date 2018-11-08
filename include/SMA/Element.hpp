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
	std::memcpy((void*)element, (void*)elementBytes, 16);
    }

    Element(uint32_t i, uint32_t j, double v)
    {
	
    	util::putUint32(i, element);
    	util::putUint32(j, &element[4]);
    	util::putDouble(v, &element[8]);
    }

/** Position is in the units of the file
*
*/
	uint32_t i()
	{
		util::nextUint32(element);
	}

uint32_t j()
{

		util::nextUint32(&element[4]);
}

double val()
{

		util::nextDouble(&element[8]);
} 


private:
    uint8_t element[16];
};

} // namespace CTL::matrix

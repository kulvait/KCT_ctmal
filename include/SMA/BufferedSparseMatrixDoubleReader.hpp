#pragma once

#include "plog/Log.h"
#include "rawop.h"
#include "SMA/ElementDouble.hpp"
#include <algorithm>

namespace CTL::matrix {
/**
 * Buffered sparse matrix is a structure
 *
 */
class BufferedSparseMatrixDoubleReader
{
public:
    /**
     *Default buffer size is in number of elements. Default is 1024 elements
     *that means 1024*16 = 16384 bytes.
     */
    BufferedSparseMatrixDoubleReader(std::string triplesFile, uint32_t bufferSize = 1024)
    {
        // Check if the file exists and treat it as a structure in the
        // format (uint32,uint32,double)
        this->triplesFile = triplesFile;
        this->bufferSize = bufferSize;
        if(io::fileExists(triplesFile))
        {
            totalFileSize = io::getFileSize(triplesFile);
            if(totalFileSize % 16 != 0)
            {
                io::throwerr("The file %s is not aligned with 16 byte "
                             "element size. It seems not to be in a "
                             "correct format since the size is %lu and "
                             "size modulo 16 is %d.",
                             triplesFile.c_str(), totalFileSize, totalFileSize % 16);
            }
            numberOfElements = totalFileSize / 16;
            LOGD << io::xprintf("Openned matrix %s, that has %lu elements.", triplesFile.c_str(),
                                totalFileSize);
        } else
        {
            io::throwerr("The file %s does not exist.", triplesFile.c_str());
        }
        buffer = new uint8_t[bufferSize * 16];
        elementsInBuffer = 0;
        currentReadingPos = 0;
    }

    ~BufferedSparseMatrixDoubleReader()
    {
        if(buffer != nullptr)
            delete[] buffer;
    }

    /// Copy constructor
    BufferedSparseMatrixDoubleReader(const BufferedSparseMatrixDoubleReader& b)
        : BufferedSparseMatrixDoubleReader(b.triplesFile, b.bufferSize)
    {
        LOGD << "Caling Copy constructor of BufferedSparseMatrixDoubleReader.";
    }

    // Copy assignment
    BufferedSparseMatrixDoubleReader& operator=(const BufferedSparseMatrixDoubleReader& b)
    {
        LOGD << "Caling Copy assignment constructor of "
                "BufferedSparseMatrixDoubleReader.";
        if(&b != this) // To elegantly solve situation when assigning
                       // to itself
        {
            this->triplesFile = b.triplesFile;
            this->totalFileSize = b.totalFileSize;
            this->numberOfElements = b.numberOfElements;
            this->currentReadingPos = b.currentReadingPos;
            this->elementsInBuffer = 0;
            if(this->buffer != nullptr)
            {
                if(this->bufferSize != b.bufferSize)
                {
                    delete[] this->buffer;
                    this->buffer = nullptr;
                    this->bufferSize = b.bufferSize;
                    this->buffer = new uint8_t[b.bufferSize * 16];
                }
            } else
            {
                this->bufferSize = b.bufferSize;
                this->buffer = new uint8_t[b.bufferSize * 16];
            }
        }
	return *this;
    }

    // Move constructor
    BufferedSparseMatrixDoubleReader(BufferedSparseMatrixDoubleReader&& b)
    {
        LOGD << "Caling Move constructor of BufferedSparseMatrixDoubleReader.";
        this->triplesFile = b.triplesFile;
        this->totalFileSize = b.totalFileSize;
        this->numberOfElements = b.numberOfElements;
        this->currentReadingPos = b.currentReadingPos;
        this->elementsInBuffer = b.elementsInBuffer;
        this->bufferSize = b.bufferSize;
        this->buffer = b.buffer;
        b.buffer = nullptr;
        this->startOfBufferPos = b.startOfBufferPos;
        this->currentBufferOffset = b.currentBufferOffset;
    }

    // Move assignment
    // Non const since I want to set buffer to nullptr
    BufferedSparseMatrixDoubleReader& operator=(BufferedSparseMatrixDoubleReader&& b)
    {
        LOGD << "Caling Move assignment constructor of "
                "BufferedSparseMatrixDoubleReader.";
        if(&b != this) // To elegantly solve situation when assigning
                       // to itself
        {
            this->triplesFile = b.triplesFile;
            this->totalFileSize = b.totalFileSize;
            this->numberOfElements = b.numberOfElements;
            this->currentReadingPos = b.currentReadingPos;
            this->elementsInBuffer = b.elementsInBuffer;
            this->startOfBufferPos = b.startOfBufferPos;
            this->currentBufferOffset = b.currentBufferOffset;
            this->bufferSize = b.bufferSize;
            if(this->buffer != nullptr)
            {
                delete[] this->buffer;
                this->buffer = nullptr;
            }
            this->buffer = b.buffer;
            b.buffer = nullptr;
        }
	return *this;
    }
    // Increments position
    void readNextValue(uint32_t* i, uint32_t* j, double* v)
    {
        if(currentBufferOffset < elementsInBuffer * 16)
        {
            *i = util::nextUint32(&buffer[currentBufferOffset]);
            *j = util::nextUint32(&buffer[currentBufferOffset + 4]);
            *v = util::nextDouble(&buffer[currentBufferOffset + 8]);
            currentBufferOffset += 16;
            currentReadingPos++;
        } else
        {
            // PAGE FAULT
            if(numberOfElements - currentReadingPos < bufferSize)
            {
                elementsInBuffer = numberOfElements - currentReadingPos;
            } else
            {
                elementsInBuffer = bufferSize;
            }
            //	LOGD << io::xprintf("Page fault on reading pos %lu increasing buffer by %d
            //elements.", currentReadingPos, elementsInBuffer);
            io::readBytesFrom(triplesFile, currentReadingPos * 16, buffer, elementsInBuffer * 16);
            currentBufferOffset = 0;
            startOfBufferPos = currentReadingPos;
            readNextValue(i, j, v);
        }
    }

    ElementDouble readNextElement()
    {
        if(currentBufferOffset < elementsInBuffer * 16)
        {
            ElementDouble e(&buffer[currentBufferOffset]);
            currentBufferOffset += 16;
            currentReadingPos++;
            return e;
        } else
        {
            // PAGE FAULT
            if(numberOfElements - currentReadingPos < bufferSize)
            {
                elementsInBuffer = numberOfElements - currentReadingPos;
            } else
            {
                elementsInBuffer = bufferSize;
            }
            //	LOGD << io::xprintf("Page fault on reading pos %lu increasing buffer by %d
            //elements.", currentReadingPos, elementsInBuffer);
            io::readBytesFrom(triplesFile, currentReadingPos * 16, buffer, elementsInBuffer * 16);
            currentBufferOffset = 0;
            startOfBufferPos = currentReadingPos;
            return readNextElement();
        }
    }

    /** Position is in the units of the number of elements
     *
     */
    void seek(uint64_t pos)
    {
        if(elementsInBuffer != 0)
        {
            if(startOfBufferPos <= pos && pos < startOfBufferPos + elementsInBuffer)
            {
                currentBufferOffset = (pos - startOfBufferPos) * 16;
            } else
            {
                elementsInBuffer = 0;
                currentBufferOffset = 0;
            }
        }
        currentReadingPos = pos;
    }

    uint64_t getNumberOfElements() const { return numberOfElements; }

private:
    uint64_t currentReadingPos;
    uint8_t* buffer;
    uint32_t bufferSize;
    uint32_t elementsInBuffer;
    uint64_t startOfBufferPos;
    uint64_t currentBufferOffset; // Offset of the reading buffer.

    std::string triplesFile;
    uint64_t totalFileSize, numberOfElements;
};

} // namespace CTL::matrix

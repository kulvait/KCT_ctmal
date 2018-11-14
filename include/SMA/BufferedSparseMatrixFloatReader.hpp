#pragma once

#include "SMA/ElementFloat.hpp"
#include "plog/Log.h"
#include "rawop.h"
#include <algorithm>

namespace CTL::matrix {
/**
 * Buffered sparse matrix is a structure
 *
 */
class BufferedSparseMatrixFloatReader
{
public:
    /**
     *Default buffer size is in number of elements. Default is 1024 elements
     *that means 1024*12 = 16384 bytes.
     */
    BufferedSparseMatrixFloatReader(std::string triplesFile, uint32_t bufferSize = 1024)
    {
        // Check if the file exists and treat it as a structure in the
        // format (uint32,uint32,double)
        this->triplesFile = triplesFile;
        this->bufferSize = bufferSize;
        if(io::fileExists(triplesFile))
        {
            totalFileSize = io::getFileSize(triplesFile);
            if(totalFileSize % 12 != 0)
            {
                io::throwerr("The file %s is not aligned with 12 byte "
                             "element size. It seems not to be in a "
                             "correct format since the size is %lu and "
                             "size modulo 12 is %d.",
                             triplesFile.c_str(), totalFileSize, totalFileSize % 12);
            }
            numberOfElements = totalFileSize / 12;
            LOGD << io::xprintf("Openned matrix %s, that has %lu elements.", triplesFile.c_str(),
                                totalFileSize);
        } else
        {
            io::throwerr("The file %s does not exist.", triplesFile.c_str());
        }
        buffer = new uint8_t[bufferSize * 12];
        elementsInBuffer = 0;
        currentReadingPos = 0;
    }

    ~BufferedSparseMatrixFloatReader()
    {
        if(buffer != nullptr)
            delete[] buffer;
    }

    /// Copy constructor
    BufferedSparseMatrixFloatReader(const BufferedSparseMatrixFloatReader& b)
        : BufferedSparseMatrixFloatReader(b.triplesFile, b.bufferSize)
    {
        LOGD << "Caling Copy constructor of BufferedSparseMatrixFloatReader.";
    }

    // Copy assignment
    BufferedSparseMatrixFloatReader& operator=(const BufferedSparseMatrixFloatReader& b)
    {
        LOGD << "Caling Copy assignment constructor of "
                "BufferedSparseMatrixFloatReader.";
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
                    this->buffer = new uint8_t[b.bufferSize * 12];
                }
            } else
            {
                this->bufferSize = b.bufferSize;
                this->buffer = new uint8_t[b.bufferSize * 12];
            }
        }
    }

    // Move constructor
    BufferedSparseMatrixFloatReader(BufferedSparseMatrixFloatReader&& b)
    {
        LOGD << "Caling Move constructor of BufferedSparseMatrixFloatReader.";
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
    BufferedSparseMatrixFloatReader& operator=(BufferedSparseMatrixFloatReader&& b)
    {
        LOGD << "Caling Move assignment constructor of "
                "BufferedSparseMatrixFloatReader.";
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
    }
    void getNextValue(uint32_t* i, uint32_t* j, float* v)
    {
        if(currentBufferOffset < elementsInBuffer * 12)
        {
            /*            *i = util::nextUint32(&buffer[currentBufferOffset]);
             *j = util::nextUint32(&buffer[currentBufferOffset + 4]);
             *v = util::nextFloat(&buffer[currentBufferOffset + 8]);
             */
            std::memcpy(i, buffer + currentBufferOffset, 4);
            std::memcpy(j, buffer + currentBufferOffset + 4, 4);
            std::memcpy(v, buffer + currentBufferOffset + 8, 4);

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
            if(elementsInBuffer == 0)
            {
                *i = 256 * 256 * 199;
                *j = 0;
                *v = 0;
                return;
            }
            //	LOGD << io::xprintf("Page fault on reading pos %lu increasing buffer by %d
            // elements.", currentReadingPos, elementsInBuffer);
            io::readBytesFrom(triplesFile, currentReadingPos * 12, buffer, elementsInBuffer * 12);
            currentBufferOffset = 0;
            startOfBufferPos = currentReadingPos;
            getNextValue(i, j, v);
        }
    }

    void increaseCounter()
    {
        currentBufferOffset += 12;
        currentReadingPos++;
    }

    // Increments position
    void readNextValue(uint32_t* i, uint32_t* j, float* v)
    {
        if(currentBufferOffset < elementsInBuffer * 12)
        {
            /*            *i = util::nextUint32(&buffer[currentBufferOffset]);
             *j = util::nextUint32(&buffer[currentBufferOffset + 4]);
             *v = util::nextFloat(&buffer[currentBufferOffset + 8]);
             */
            std::memcpy(i, buffer + currentBufferOffset, 4);
            std::memcpy(j, buffer + currentBufferOffset + 4, 4);
            std::memcpy(v, buffer + currentBufferOffset + 8, 4);

            currentBufferOffset += 12;
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
            if(elementsInBuffer == 0)
            {
                *i = 256 * 256 * 199;
                *j = 0;
                *v = 0;
                return;
            }
            //	LOGD << io::xprintf("Page fault on reading pos %lu increasing buffer by %d
            // elements.", currentReadingPos, elementsInBuffer);
            io::readBytesFrom(triplesFile, currentReadingPos * 12, buffer, elementsInBuffer * 12);
            currentBufferOffset = 0;
            startOfBufferPos = currentReadingPos;
            readNextValue(i, j, v);
        }
    }

    ElementFloat readNextElement()
    {
        if(currentBufferOffset < elementsInBuffer * 12)
        {
            ElementFloat e(&buffer[currentBufferOffset]);
            currentBufferOffset += 12;
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
            // elements.", currentReadingPos, elementsInBuffer);
            io::readBytesFrom(triplesFile, currentReadingPos * 12, buffer, elementsInBuffer * 12);
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
                currentBufferOffset = (pos - startOfBufferPos) * 12;
            } else
            {
                elementsInBuffer = 0;
                currentBufferOffset = 0;
            }
        }
        currentReadingPos = pos;
    }

    bool atEnd() { return (currentReadingPos == numberOfElements); }

    uint64_t getNumberOfElements() const { return numberOfElements; }

private:
    uint64_t currentReadingPos;
    uint8_t* buffer;
    int bufferSize;
    int elementsInBuffer;
    uint64_t startOfBufferPos;
    uint64_t currentBufferOffset; // Offset of the reading buffer.

    std::string triplesFile;
    uint64_t totalFileSize, numberOfElements;
};

} // namespace CTL::matrix

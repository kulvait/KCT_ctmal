#pragma once

#include <mutex>
#include "plog/Log.h"

namespace CTL::matrix {
/**
 * Buffered sparse matrix is a structure
 *
 */
class BufferedSparseMatrixWritter
{
public:
    /**
     */
    BufferedSparseMatrixWritter(std::string triplesFile,
                                uint32_t bufferSize=1024,
                                bool overwrite = false)
    {
        // Check if the file exists and treat it as a structure in the
        // format (uint32,uint32,double)
        this->triplesFile = triplesFile;
        this->bufferSize = bufferSize;
        if(io::fileExists(triplesFile))
        {
            if(overwrite)
            {
                LOGW << io::xprintf("Overwritting the file %s with a empty file.",
                                    triplesFile.c_str());
                io::createEmptyFile(triplesFile, 0, overwrite);
                numberOfElements = 0;
            } else
            {

                uint64_t totalFileSize = io::getFileSize(triplesFile);
                if(fileSize % 16 != 0)
                {
                    io::throwerr("The file %s is not aligned with 16 byte "
                                 "element size. It seems not to be in a "
                                 "correct format since the size is %lu and "
                                 "size modulo 16 is %d.",
                                 triplesFile.c_str(), fileSize, fileSize % 16);
                }
                numberOfElements = totalFileSize / 16;
                LOGD << io::xprintf("Openned matrix %s, that has %lu elements.",
                                    projectionsFile.c_str(), totalFileSize);
            }
        } else
        {
            io::createEmptyFile(triplesFile, 0, overwrite);
            numberOfElements = 0;
        }
        buffer = new uint8_t[bufferSize * 16];
        elementsInBuffer = 0;
        currentBufferPos = 0;
    }

    ~BufferedSparseMatrixWritter()
    {
        if(buffer != nullptr)
        {
            flush();
            delete[] buffer;
        }
    }

    /// Copy constructor
    BufferedSparseMatrixWritter(const BufferedSparseMatrixWritter& b)
    {
        LOGW << "Caling Copy constructor of BufferedSparseMatrixWritter.";
        b.flush();
        this->triplesFile = b.triplesFile;
        this->bufferSize = b.bufferSize;
        this->numberOfElements = b.numberOfElements;
        buffer = new uint8_t[bufferSize * 16];
        elementsInBuffer = 0;
        currentBufferPos = 0;
    }

    // Copy assignment
    BufferedSparseMatrixWritter& operator=(const BufferedSparseMatrixWritter& b)
    {
        LOGW << "Caling Copy assignment constructor of "
                "BufferedSparseMatrixWritter.";
        b.flush();
        if(&b != this) // To elegantly solve situation when assigning
                       // to itself
        {
            this->flush();
            this->triplesFile = b.triplesFile;
            this->numberOfElements = b.numberOfElements;
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
        this->elementsInBuffer = 0;
        this->currentBufferPos = 0;
    }

    // Move constructor
    BufferedSparseMatrixWritter(BufferedSparseMatrixWritter&& b)
    {
        LOGW << "Caling Move constructor of BufferedSparseMatrixWritter.";
        b.flush();
        this->triplesFile = b.triplesFile;
        this->numberOfElements = b.numberOfElements;
        this->buffer = b.buffer;
        b.buffer = nullptr;
        this->bufferSize = b.bufferSize;
        this->elementsInBuffer = 0;
        this->currentBufferPos = 0;
    }

    // Move assignment
    BufferedSparseMatrixWritter& operator=(const BufferedSparseMatrixWritter&& b)
    {
        LOGW << "Caling Move assignment constructor of "
                "BufferedSparseMatrixWritter.";
        this->flush();
        b.flush();
        if(&b != this) // To elegantly solve situation when assigning
                       // to itself
        {
            this->triplesFile = b.triplesFile;
            this->numberOfElements = b.numberOfElements;
            if(this->buffer != nullptr)
            {
                delete[] this->buffer;
                this->buffer = nullptr;
            }
            this->buffer = b.buffer;
            b.buffer = nullptr;
            this->bufferSize = b.bufferSize;
            this->elementsInBuffer = 0;
            this->currentBufferPos = 0;
        }
    }
    // Buffer into disk
    void flush()
    {
        std::lock_guard<std::mutex> guard(
            writingMutex); // Mutex will be released as this goes out of scope.
        if(elementsInBuffer != 0)
        {
            io::appendBytes(triplesFile, buffer, elementsInBuffer);
            currentBufferPos = 0;
            numberOfElements += elementsInBuffer;
            elementsInBuffer = 0;
        }
    }

    void insertValue(uint32_t i, uint32_t j, double v)
    {

        std::lock_guard<std::mutex> guard(
            writingMutex); // Mutex will be released as this goes out of scope.
        if(elementsInBuffer < bufferSize)
        {
            util::putUint32(i, &buffer[currentBufferPos]);
            util::putUint32(j, &buffer[currentBufferPos + 4]);
            util::putDouble(v, &buffer[currentBufferPos + 8]);
            currentBufferPos += 16;
            elementsInBuffer++;
        } else
        {
            // PAGE FAULT
            flush();
            insertValue(i, j, v);// ... this would not work with normal mutex
        }
    }

private:
    uint8_t* buffer;
    int bufferSize;
    int elementsInBuffer;
    int currentBufferPos;

    std::string triplesFile;
    uint64_t numberOfElements;
        mutable std::recursive_mutex writingMutex;
};

} // namespace CTL::matrix

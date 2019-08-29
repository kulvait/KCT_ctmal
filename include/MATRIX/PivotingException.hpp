#pragma once

#include <exception>

namespace CTL {
namespace matrix {
    class PivotingException : public std::exception
    {
        std::string msg;

    public:
        PivotingException(std::string msg)
            : msg(msg)
        {
        }

        const char* what() { return msg.c_str(); }
    };

} // namespace matrix
} // namespace CTL

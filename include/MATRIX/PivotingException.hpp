#pragma once

#include <exception>

namespace KCT {
namespace matrix {
    class PivotingException : public std::exception
    {
        const std::string msg;

    public:
        PivotingException(std::string msg)
            : msg(msg)
        {
        }

        const char* what() const throw() { return msg.c_str(); }
    };

} // namespace matrix
} // namespace KCT

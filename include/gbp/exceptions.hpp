/*  This file is part of libGBP - http://www.libgbp.org/
 *
 *  Copyright (c) 2006-2011, The libGBP authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */


/// \file
/// \brief Defines the Exception class and macros for throwing exceptions and doing assertions


#ifndef EXCEPTIONS_HPP_
#define EXCEPTIONS_HPP_


#include <exception>
#include <stdexcept>
#include <string>
#include <iostream>


/// Used by GBP_THROW
#define GBP_QUOTE(x) #x

/// Used by GBP_THROW
#define GBP_TOSTRING(x) GBP_QUOTE(x)

/// Macro that simplifies throwing an exception with a useful default error message. 
/** The error message consists of a description of the exception, the source 
 *  code file and line number where the exception has been thrown.
 *  \param cod Corresponds to one of the enum values of gbp::Exception::Code
 *
 *  \par Example:
 *  \code
 *  GBP_THROW(NOT_IMPLEMENTED);
 *  \endcode
 */
#if defined __GNUG__ // GNU C++
  #define FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined _MSC_VER // Visual Studio
  #define FUNCTION_NAME __FUNCTION__
#else // other compilers
  #define FUNCTION_NAME __func__
#endif
#define GBP_THROW(cod) throw gbp::Exception(gbp::Exception::cod, __FILE__, FUNCTION_NAME, GBP_TOSTRING(__LINE__), "")

/// Macro that simplifies throwing an exception with a user-defined error message.
/** \param cod Corresponds to one of the enum values of gbp::Exception::Code
 *  \param msg Detailed error message that will be written to std::cerr.
 *
 *  \par Example:
 *  \code
 *  GBP_THROWE(NOT_IMPLEMENTED,"Detailed error message");
 *  \endcode
 */
#define GBP_THROWE(cod,msg) throw gbp::Exception(gbp::Exception::cod, __FILE__, FUNCTION_NAME, GBP_TOSTRING(__LINE__), msg)

/// Assertion mechanism, similar to the standard assert() macro. It is always active, even if NDEBUG is defined
#define GBP_ASSERT(condition) ((condition) ? ((void)0) : GBP_THROWE(ASSERTION_FAILED, std::string("Assertion \"" #condition "\" failed")))

// Assertion only if GBP_DEBUG is defined
#ifdef GBP_DEBUG
/// Assertion mechanism similar to GBP_ASSERT which is only active if GBP_DEBUG is defined
#define GBP_DEBASSERT(x) do {GBP_ASSERT(x);} while(0)
#else
#define GBP_DEBASSERT(x) do {} while(0)
#endif


namespace gbp {


/// Error handling in libGBP is done by throwing an instance of the Exception class.
/** The Exception class inherits from std::runtime_error. It defines several types of exceptions
 *  and corresponding error messages. The recommended way to throw an instance of the Exception
 *  class is by using the #GBP_THROW or #GBP_THROWE macros.
 */
class Exception : public std::runtime_error {
    public:
        /// Enumeration of exceptions used in libGBP
        enum Code {NOT_IMPLEMENTED,
                   ASSERTION_FAILED,
                   IMPOSSIBLE_TYPECAST,
                   OBJECT_NOT_FOUND,
                   BELIEF_NOT_AVAILABLE,
                   UNKNOWN_ENUM_VALUE,
                   UNKNOWN_GBP_ALGORITHM,
                   UNKNOWN_PARAMETER_ESTIMATION_METHOD,
                   UNKNOWN_PROPERTY_TYPE,
                   UNKNOWN_PROPERTY,
                   MALFORMED_PROPERTY,
                   NOT_ALL_PROPERTIES_SPECIFIED,
                   INVALID_ALIAS,
                   CANNOT_READ_FILE,
                   CANNOT_WRITE_FILE,
                   INVALID_FACTORGRAPH_FILE,
                   INVALID_EVIDENCE_FILE,
                   INVALID_EMALG_FILE,
                   NOT_NORMALIZABLE,
                   MULTIPLE_UNDO,
                   FACTORGRAPH_NOT_CONNECTED,
                   INTERNAL_ERROR,
                   RUNTIME_ERROR,
                   OUT_OF_MEMORY,
                   NUM_ERRORS};  // NUM_ERRORS should be the last entry

        /// Constructor
        Exception( Code code, const char *filename, const char *function, const char *line, const std::string& detailedMsg ) :
            std::runtime_error(ErrorStrings[code] + (detailedMsg.empty() ? "" : (": " + detailedMsg)) + " [File " + filename + ", line " + line + ", function: " + function + "]"), 
            _errorcode(code), _detailedMsg(detailedMsg), _filename(filename), _function(function), _line(line) {}

        /// Destructor
        ~Exception() throw () {}

        /// Returns error code of this exception
        Code getCode() const { return _errorcode; }

        /// Returns error code of this exception
        /** \deprecated Please use gbp::Exceptions::getCode() instead
         */
        Code code() const { return getCode(); }

        /// Returns short error message of this exception
        const std::string& getMsg() const { return ErrorStrings[_errorcode]; }

        /// Returns detailed error message of this exception
        const std::string& getDetailedMsg() const { return _detailedMsg; }

        /// Returns filename where this exception was thrown
        const std::string& getFilename() const { return _filename; }

        /// Returns function name in which this exception was thrown
        const std::string& getFunction() const { return _function; }

        /// Returns line number where this exception was thrown
        const std::string& getLine() const { return _line; }

        /// Returns error message corresponding to an error code
        const std::string& message( const Code c ) const { return ErrorStrings[c]; }

    private:
        /// Contains the error code of this exception
        Code _errorcode;

        /// Contains the detailed message of this exception, if any
        std::string _detailedMsg;
        
        /// Contains the filename where this exception was thrown
        std::string _filename;

        /// Contains the function name in which this exception was thrown
        std::string _function;

        /// Contains the line number where this exception was thrown
        std::string _line;

        /// Error messages corresponding to the exceptions enumerated above
        static std::string ErrorStrings[NUM_ERRORS];
};


}


#endif

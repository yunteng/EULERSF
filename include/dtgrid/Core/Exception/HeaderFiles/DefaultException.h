/*************************************************************************************************
 *
 * Copyright (c) 2009, Michael Bang Nielsen (nielsenmb@gmail.com), Aarhus University, Denmark.
 * All rights reserved.
 *
 *************************************************************************************************/

#ifndef _core_exception_defaultexception_h_
#define _core_exception_defaultexception_h_

#include <string>
#include <sstream>
#include <iostream>

namespace Core
{

    class DefaultException : public std::exception
    {
    public:
	~DefaultException() throw ()
	{
	}

	DefaultException(const std::string& message) throw ()
	{
	    this->message = message;
	}

	virtual const char *what( ) const throw ()
	{
	    return message.c_str();
	}

    protected:
	std::string message;
    };


    inline void throwDefaultException(const std::string& message, const std::string& file, int line)
    {
	std::stringstream ss;

	std::cerr << "Error at line " << line << " in file " << file << ": " << message << std::endl; 
	ss << "Error at line " << line << " in file " << file << ": " << message << std::endl; 

	throw DefaultException(ss.str());
    }
}

#endif

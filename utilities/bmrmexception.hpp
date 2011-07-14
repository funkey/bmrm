#ifndef _BMRMEXCEPTION_HPP_
#define _BMRMEXCEPTION_HPP_

#include <exception>
#include <string>

class CBMRMException : public std::exception 
{
public:
	CBMRMException(const std::string& theMessage, const std::string& theThrower);
	virtual ~CBMRMException() throw();
	virtual const std::string& Report();

protected:
	std::string message;
	std::string thrower;      
};


#endif

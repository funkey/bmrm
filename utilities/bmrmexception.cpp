#ifndef _BMRMEXCEPTION_CPP_
#define _BMRMEXCEPTION_CPP_

#include "bmrmexception.hpp"
#include <sstream>

using namespace std;

/**   BMRM Exception class constructor
 *
 *    \param theMessage [read] The error message
 *    \param theThrower [read] Description of the object throwing this error
 */
CBMRMException::CBMRMException(const string& theMessage, const string& theThrower="NONAME")
{
   message = theMessage;
   thrower = theThrower;
}

/**   BMRM Exception class destructor
 */
CBMRMException::~CBMRMException() throw()
{}


/**   Report error message and related information
 */
const string& CBMRMException::Report()
{
   ostringstream ostr;
   ostr << "[BMRM error message] (Thrower: " 
	<< thrower 
	<< ")\n" 
	<< message 
	<< "\n";

   message = ostr.str();
   return message;
}

#endif

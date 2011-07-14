#ifndef _COMMON_CPP_
#define _COMMON_CPP_ 

#include <string>
#include <vector>
#include <cctype>
#include <math.h>

/** Check if a string is blank.
 *
 *  \param line [read] The line to be checked if it is blank
 */
bool IsBlankLine(std::string &line)
{
    size_t n = line.size();
    for(size_t i=0; i<n; i++)
        if(!std::isspace(line[i]))
            return false;
    return true;
}


void trim(std::string& str)
{  
    typedef std::string::size_type str_pos;
    str_pos pos = str.find_last_not_of(" \t\r");
    if(pos != std::string::npos) 
    {
        str.erase(pos + 1);
        pos = str.find_first_not_of(" \t\r");
        if(pos != std::string::npos) 
            str.erase(0, pos);
    }
    else 
        str.erase(str.begin(), str.end());
}


void tokenize(const std::string& str, 
              std::vector<std::string>& tokens, 
              const std::string& delimiter = " ") 
{
    typedef std::string::size_type str_pos;
    
    // Skip delimiter at beginning.
    str_pos lastPos = str.find_first_not_of(delimiter, 0);
    // Find first "non-delimiter".
    str_pos pos = str.find_first_of(delimiter, lastPos);
    
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiter.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiter, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiter, lastPos);
    }
}


// adopted from Numerical Recipe

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <=0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


float gasdev(long *idum)
{
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;
    
    if(*idum < 0) iset=0;
    if (iset == 0) {
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}

#endif

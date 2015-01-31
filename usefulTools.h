/* 
 * File:   usefulTools.h
 * Author: Choi
 *
 * Created on February 18, 2012, 5:41 PM
 */

#ifndef USEFULTOOLS_H
#define	USEFULTOOLS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>


class usefulTools{
    public:
    static std::string clean(const std::string seq);
    static std::string trim(const std::string seq);
    static void tokenizer(const std::string seq, const std::string separators, std::vector<std::string>* result);
    static bool fileExists(const std::string fileName) ;
    static bool checkIfNumeric(const std::string seq);
	static double dnorm(const double x);
	static double qnorm(const double p);
};

#endif	/* USEFULTOOLS_H */

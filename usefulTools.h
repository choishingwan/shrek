// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.

#ifndef USEFULTOOLS_H
#define	USEFULTOOLS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>

/** \class usefulTools
 *  \brief a collection of small utility functions gathered in different places
 */
class usefulTools{
    public:
    /** Remove special characters from the end of the string (by Thomas, my FYP supervisor) */
    static std::string clean(const std::string seq);
    /** Remove special characters from the front and end of the string */
    static std::string trim(const std::string seq);
    /** Tokenize the string using based on the separator (by Thomas, my FYP supervisor) */
    static void tokenizer(const std::string seq, const std::string separators, std::vector<std::string>* result);
    /** Check if file exists */
    static bool fileExists(const std::string fileName) ;
    /** Check if the string is a number */
    static bool isNumeric(const std::string seq);
    /** dnorm from R */
	static double dnorm(const double x);
	/** qnorm from R */
	static double qnorm(const double p);
    /* Black magic got from http://stackoverflow.com/a/109025 */
    /* Replacing pop_count, should not be platform dependent */
    static int NumberOfSetBits(uint32_t i);

    /** Functions for getting the sign of a value input, obtained from  <a href="http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c">stackoverflow</a>  */
	template <typename T> inline constexpr
	static int signum(T x, std::false_type is_signed) { return T(0) < x; }
    /** Functions for getting the sign of a value input, obtained from  <a href="http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c">stackoverflow</a>  */
	template <typename T> inline constexpr
	static int signum(T x, std::true_type is_signed) { return (T(0) < x) - (x < T(0)); }
    /** Functions for getting the sign of a value input, obtained from  <a href="http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c">stackoverflow</a>  */
	template <typename T> inline constexpr
	static int signum(T x) { return signum(x, std::is_signed<T>()); }
};

#endif	/* USEFULTOOLS_H */


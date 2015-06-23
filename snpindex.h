// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.
#ifndef SNPINDEX_H
#define SNPINDEX_H

#include <map>
#include <string>
#include <iostream>

/** \class SnpIndex
 *	\brief a special storage class
 *
 *	Essentially, this is only a map structure. But we use a class
 *	to store this information so that in the future, if someone
 *	would like to optimize the code, they can quickly change this
 *	slow map usage into something faster e.g. Hash
 */

class SnpIndex
{
	public:
		/** Default constructor */
		SnpIndex();
		/**Default destructor */
		virtual ~SnpIndex();
		/** Get the value of the current position of m_iter */
		size_t value() const;
		/** Get the value of a specific key */
		size_t value(std::string key);
        /** Get the current size of the map */
		size_t size() const;
		/** Set the value of the key */
		void set(std::string key, size_t value);
		/** Initialize the m_iter to the front of the map */
		void init();
		/** Add one to the value of the current key */
		void increment(std::string key);
		/** Find if the key is within the map */
		bool contains(std::string key);
		/** Check if we have reach the end of the map */
		bool valid();
		/** Got to the next key */
		bool next();
        /** Print the map */
		void print();
	protected:
	private:
		bool m_isInitialized;
		std::map<std::string, size_t> m_index;
		std::map<std::string, size_t>::iterator m_iter;
};

#endif // SNPINDEX_H

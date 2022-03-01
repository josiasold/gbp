/*  This file is part of libGBP - http://www.libGBP.org/
 *
 *  Copyright (c) 2006-2011, The libGBP authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */

/// \file
/// \brief Defines class Var, which represents a discrete random variable.

#ifndef VAR_HPP_
#define VAR_HPP_

#include <iostream>
#include <sstream>
#include "exceptions.hpp"

namespace gbp
{

        /// Represents a discrete random variable.
        /** A Var stores the \a label of the variable (an unsigned integer-valued
         *  unique ID) and the number of possible values (\a states) of that variable. 
         *  Two Var objects with the same label are assumed to be identical (i.e., it
         *  is assumed that they have the same number of possible states).
         *
         *  In the documentation, we use the following notational conventions. The discrete
         *  random variable with label \f$l\f$ is denoted as \f$x_l\f$, and the number
         *  of possible values of this variable as \f$S_l\f$; this is represented in
         *  code by the object Var(\f$l\f$,\f$S_l\f$). The set of possible values of
         *  variable \f$x_l\f$ is denoted \f$X_l := \{0,1,\dots,S_l-1\}\f$.
         */
        class Var
        {
        private:
                /// Label of the variable (its unique ID)
                size_t p_label;

                /// Number of possible values
                size_t p_states;

        public:
                /// Default constructor (creates a variable with label 0 and 0 states)
                Var() : p_label(0), p_states(0) {}
                /// Constructs a variable with a given label and number of states
                Var(size_t label, size_t states) : p_label(label), p_states(states) {}

                /// Returns the label
                size_t label() const { return p_label; }
                /// Returns reference to label
                size_t &label() { return p_label; }

                /// Returns the number of states
                size_t states() const { return p_states; }
                /// Returns reference to number of states
                size_t &states() { return p_states; }

                /// Smaller-than operator (only compares labels)
                bool operator<(const Var &n) const
                {
#ifdef GBP_DEBUG
                        if (p_label == n.p_label)
                                GBP_ASSERT(p_states == n.p_states);
#endif
                        return (p_label < n.p_label);
                }

                /// Larger-than operator (only compares labels)
                bool operator>(const Var &n) const
                {
#ifdef GBP_DEBUG
                        if (p_label == n.p_label)
                                GBP_ASSERT(p_states == n.p_states);
#endif
                        return (p_label > n.p_label);
                }

                /// Smaller-than-or-equal-to operator (only compares labels)
                bool operator<=(const Var &n) const
                {
#ifdef GBP_DEBUG
                        if (p_label == n.p_label)
                                GBP_ASSERT(p_states == n.p_states);
#endif
                        return (p_label <= n.p_label);
                }

                /// Larger-than-or-equal-to operator (only compares labels)
                bool operator>=(const Var &n) const
                {
#ifdef GBP_DEBUG
                        if (p_label == n.p_label)
                                GBP_ASSERT(p_states == n.p_states);
#endif
                        return (p_label >= n.p_label);
                }

                /// Not-equal-to operator (only compares labels)
                bool operator!=(const Var &n) const
                {
#ifdef GBP_DEBUG
                        if (p_label == n.p_label)
                                GBP_ASSERT(p_states == n.p_states);
#endif
                        return (p_label != n.p_label);
                }

                /// Equal-to operator (only compares labels)
                bool operator==(const Var &n) const
                {
#ifdef GBP_DEBUG
                        if (p_label == n.p_label)
                                GBP_ASSERT(p_states == n.p_states);
#endif
                        return (p_label == n.p_label);
                }

                /// Writes a Var to an output stream
                friend std::ostream &operator<<(std::ostream &os, const Var &n)
                {
                        return (os << "x" << n.label());
                }

                /// Formats a Var as a string
                std::string toString() const
                {
                        std::stringstream ss;
                        ss << *this;
                        return ss.str();
                }
        };

} // end of namespace gbp

#endif

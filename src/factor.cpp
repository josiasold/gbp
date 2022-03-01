/*  This file is part of libGBP - http://www.libGBP.org/
 *
 *  Copyright (c) 2006-2011, The libGBP authors. All rights reserved.
 *
 *  Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.
 */

#include <gbp/factor.hpp>

namespace gbp
{

    using namespace std;

    Factor createFactorIsing(const Var &n, Real h)
    {
        GBP_ASSERT(n.states() == 2);
        Real buf[2];
        buf[0] = std::exp(-h);
        buf[1] = std::exp(h);
        return Factor(n, &buf[0]);
    }

    Factor createFactorIsing(const Var &n1, const Var &n2, Real J)
    {
        GBP_ASSERT(n1.states() == 2);
        GBP_ASSERT(n2.states() == 2);
        GBP_ASSERT(n1 != n2);
        Real buf[4];
        buf[0] = (buf[3] = std::exp(J));
        buf[1] = (buf[2] = std::exp(-J));
        return Factor(VarSet(n1, n2), &buf[0]);
    }

    Factor createFactorExpGauss(const VarSet &ns, Real beta)
    {
        Factor fac(ns);
        for (size_t t = 0; t < fac.nrStates(); t++)
            fac.set(t, std::exp(rnd_stdnormal() * beta));
        return fac;
    }

    Factor createFactorPotts(const Var &n1, const Var &n2, Real J)
    {
        Factor fac(VarSet(n1, n2), 1.0);
        GBP_ASSERT(n1.states() == n2.states());
        for (size_t s = 0; s < n1.states(); s++)
            fac.set(s * (n1.states() + 1), std::exp(J));
        return fac;
    }

    Factor createFactorDelta(const Var &v, size_t state)
    {
        Factor fac(v, 0.0);
        GBP_ASSERT(state < v.states());
        fac.set(state, 1.0);
        return fac;
    }

    Factor createFactorDelta(const VarSet &vs, size_t state)
    {
        Factor fac(vs, 0.0);
        GBP_ASSERT((BigInt)state < vs.nrStates());
        fac.set(state, 1.0);
        return fac;
    }

} // end of namespace GBP

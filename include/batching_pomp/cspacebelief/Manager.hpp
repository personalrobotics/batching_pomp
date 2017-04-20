/***********************************************************************
Copyright (c) 2017, Carnegie Mellon University
All rights reserved.
Authors: Shushman Choudhury <shushmanchoudhury@gmail.com>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

#ifndef BATCHING_POMP_CSPACEBELIEF_MANAGER_HPP_
#define BATCHING_POMP_CSPACEBELIEF_MANAGER_HPP_

namespace batching_pomp {
namespace cspacebelief {

/// Abstract class that represents the configuration space belief manager.
/// It is templated on the type of data structure used to represent
/// a configuration space point and its result
template <class _T>
class Manager
{
public:

    virtual ~Manager() = default;

    /// Add a data point to the belief manager
    virtual void addPoint(const _T &data) = 0;

    /// Remove a data point from the belief manager (rarely used)
    virtual void removePoint(const _T &data) = 0;

    /// Estimate the configuration space belief of the query
    virtual double estimate(const _T &query) const = 0;

};

} //namespace cspacebelief
} //namespace ompl_pomp

#endif //BATCHING_POMP_CSPACEBELIEF_MANAGER_HPP_
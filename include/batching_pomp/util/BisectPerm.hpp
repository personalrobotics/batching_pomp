/***********************************************************************
Copyright (c) 2015, Carnegie Mellon University
All rights reserved.
Authors: Chris Dellin <cdellin@gmail.com>

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

#ifndef BATCHING_POMP_UTIL_BISECTPERM_HPP_
#define BATCHING_POMP_UTIL_BISECTPERM_HPP_

namespace batching_pomp {
namespace util {
 
//! Generates a Van Der Corput sequence ordering for states to check along an edge.

class BisectPerm
{
public:
  BisectPerm() {}

  /// For an edge that has n states, we generate a sequence of n fractional
  /// positions along the edge to check for collision, based on Van Der Corput Sequences
  /// We start at 1/2, then 1/4, and 3/4 and so on upto n states
  /// \param[in] n The number of states along the edge
  /// \return A map from integer index to fractional position along edge.
  const std::vector< std::pair<int,int> > & get(int n)
  {
    std::map<int, const std::vector< std::pair<int,int> > >::iterator it;
    it = cache.find(n);
    if (it != cache.end())
      return it->second;

    int i;
    int last_true;
    int max_i;
    int max_val;
    std::vector< std::pair<int,int> > perm;
    std::vector<bool> done(n, false);
    std::vector<int> dist(n);

    for (;;)
    {
      last_true = -1;
      for (i=0; i<n; i++)
      {
        if (done[i]) last_true = i;
        dist[i] = (i-last_true);
      }
      last_true = n;
      for (i=n-1; i>=0; i--)
      {
        if (done[i]) last_true = i;
        dist[i] = (last_true-i) < dist[i] ? (last_true-i) : dist[i];
      }
      max_val = 0;
      max_i = 0;
      for (i=0; i<n; i++) if (max_val < dist[i])
      {
        max_val = dist[i];
        max_i = i;
      }
      if (!max_val)
        break;
      perm.push_back(std::make_pair(max_i,max_val));
      done[max_i] = true;
    }

    cache.insert(std::make_pair(n,perm));
    return cache[n];
  }

private:

  std::map<int, const std::vector< std::pair<int,int> > > cache;
};

} // namespace util
} // namespace batching_pomp

#endif // BATCHING_POMP_UTIL_BISECTPERM_HPP_
import networkx as nx
import numpy
import sys
import yaml
import math
import copy
import argparse

# Typically the first d primes are used as generators
# for a d-dimensional Halton sequence. The first 20
# primes are provided here for convenience.
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]

# The joint limits for the WAM arm (http://www.barrett.com/products-arm.htm)
# Obtained from https://github.com/personalrobotics/herb_description/blob/master/robots/wam.urdf.xacro
WAM_LOWER = numpy.array([ 0.55,-2.,-2.8,-0.9,-4.76,-1.6,-3.])
WAM_HIGHER = numpy.array([ 5.74,2.,2.8,3.1,1.24,1.6,3.])


def get_halton_value(_index,base):
  """
  Compute the halton value at an index for a base b
  Source - https://en.wikipedia.org/wiki/Halton_sequence
  @param index The index in the Halton sequence for which to generate the point
  @param base The prime generator used for the sequence
  @return res The float value of the index'th member of the sequence with base b
  """
  assert type(_index) is int
  assert type(base) is int

  res = 0
  f = 1
  index = _index
  while index > 0:
    f = f*1.0/base
    res = res + f*(index%base)
    index = index/base
  return res

def wraparound_unit(_jointvals):
  """
  Wrap the individual values of an array around 1
  @param _jointvals The array to wrap around
  @return jointvals The wrapped around numpy array
  """
  assert type(_jointvals) is numpy.ndarray

  jointvals = numpy.copy(_jointvals)
  ndims = len(jointvals)
  for i in range(ndims):
    if jointvals[i] >=1.0:
      jointvals[i] = jointvals[i]-1.0

  return jointvals

def get_halton_sequence(n,bases,lower,upper,offset,discard):
  """
  Generate a Halton sequence of n points using the d-dimensional vector
  of prime bases as generators. Then scale the sequence from the unit
  hypercube to be between the lower and upper limits. Apply some random
  offset in each dimension
  @param n The number of points of the Halton sequence to generate
  @param bases The array of prime bases for the Halton sequence
  @param lower The lower bounds of the state space
  @param upper The upper bounds of the state space
  @param offset An array of offsets to apply to each index of each point in the sequence
  @param discard The first X points in the sequence to discard
  @return list_of_scaled_pts A Python list of n points of the Halton sequence
  """
  assert (type(n) is int and n > 0)
  assert len(bases)==len(lower) and len(bases)==len(upper) and len(bases)==len(offset)
  assert discard >= 0

  diff = upper - lower
  list_of_points = [wraparound_unit(numpy.array([get_halton_value(i,base) for base in bases])+offset) for i in range(1,n+1+discard)]
  list_of_scaled_pts = [lower + numpy.multiply(pt,diff) for pt in list_of_points[discard:]]

  return list_of_scaled_pts

def get_random_sequence(n,lower,upper):
  """
  Generate a random i.i.d set of d-dimensional points between lower and upper bounds
  @param n The number of random i.i.d points to generate
  @param lower The lower bounds of the space
  @param upper The upper bounds of the space
  @return list_of_rand_pts A Python list of n random i.i.d points
  """
  assert (type(n) is int and n > 0)
  assert len(lower)==len(upper)

  diff = upper - lower
  ndims = len(lower)
  list_of_rand_pts = [lower + numpy.multiply(numpy.random.rand(ndims,),diff) for i in range(n)]

  return list_of_rand_pts

def array_to_stringstream(arr):
  """
  Convert a d-dimensional float array to a stream of characters.
  Each number has only 2 digits of precision to save space for large dense high-dimensional roadmaps
  @param arr The array to convert to character stream
  @return result The string that represents the array, with index values separated by whitespaces
  """
  assert len(arr) > 0

  result = '%.2f' % arr[0]
  dims = len(arr)
  for i in range(1,dims):
    result = result+' '+('%.2f' % arr[i])

  return result


# The script generates multi-dimensional points of either the Halton or 
# random I.I.D sequence. Then a networkx graph is created where the vertices
# have a "state" attribute, and each vertex represents a point in the sequence.
# The graph is then saved as a .graphml file which is readable by Boost Graph
# and usable by the batching_pomp planner. Edges are not saved as we assume the
# complete roadmap will be used by the planner
if __name__=="__main__":

  parser = argparse.ArgumentParser(description='Generate Halton Sequence and save .graphml file with vertices')
  parser.add_argument('--num_nodes',type=int,required=True,help='The number of points of the sequence to generate')
  parser.add_argument('--dimensions',type=int,required=True,help='The dimensionality of the points')
  parser.add_argument('--bounds',type=str,required=True,help='The kind of bounds. Options are unit or wam')
  parser.add_argument('--offset',action='store_true',help='Whether to have a random offset for each point of a Halton sequence')
  parser.add_argument('--sequence',type=str,required=True,help='The kind of sequence. Options are halton or random')
  parser.add_argument('--discard',type=int,default=500,help='The burn-in for Halton sequence. How many points to discard at beginning.')
  parser.add_argument('--outfile',type=str,required=True,help='The .graphml file to save sequence to')
  args = parser.parse_args()

  if args.bounds == 'unit':
    # For the unit hypercube
    lower = numpy.array([0.0 for i in range(args.dimensions)])
    upper = numpy.array([1.0 for i in range(args.dimensions)])
  elif args.bounds == 'wam':
    # For the Barrett WAM arm
    lower = WAM_LOWER
    upper = WAM_HIGHER

  # Use the first d prime numbers as bases
  bases = PRIMES[0:args.dimensions]

  # Only generate a random offset if the argument is given, else all zeroes
  offset = numpy.zeros(args.dimensions,)
  if args.offset:
    offset = numpy.random.random_sample(args.dimensions,)

  if args.sequence == 'halton':
    point_seq = get_halton_sequence(args.num_nodes,bases,lower,upper,offset,args.discard)
  elif args.sequence == 'random':
    point_seq = get_random_sequence(args.num_nodes,lower,upper)

  # Create networkx graph
  G = nx.Graph()

  # Add node attributes
  G.add_nodes_from(range(args.num_nodes))
  states = {i : array_to_stringstream(pt) for (i,pt) in enumerate(point_seq)}
  nx.set_node_attributes(G, 'state', states)

  # Save graph to .graphml file
  nx.write_graphml(G, args.outfile)
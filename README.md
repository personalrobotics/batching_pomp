# batching_pomp

A motion planning algorithm implemented in [OMPL](http://ompl.kavrakilab.org/), also using [Boost Graph Library](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/index.html). It is released under the [BSD License](https://opensource.org/licenses/BSD-2-Clause). The code is based on our upcoming IJRR submission. Doxygen-generated documentation is available [here](http://web.stanford.edu/~shushman/documentation/batching_pomp/index.html).
In the meantime, please cite the following papers if you are interested in using this code:

[1] *Choudhury, S., Salzman, O., Choudhury, S., & Srinivasa, S. S. (2017, May). Densification strategies for anytime motion planning over large dense roadmaps. In Robotics and Automation (ICRA), 2017 IEEE International Conference on (pp. 3770-3777). IEEE.*

[2] *Choudhury, S., Dellin, C. M., & Srinivasa, S. S. (2016, October). Pareto-optimal search over configuration space beliefs for anytime motion planning. In Intelligent Robots and Systems (IROS), 2016 IEEE/RSJ International Conference on (pp. 3742-3749). IEEE.*

This is an implementation of an algorithmic framework for anytime motion planning on large dense roadmaps. There are two key components:

- A sequence of subgraphs of the entire roadmap is generated, using some densification strategy.
- After each subgraph is added, the current roadmap is searched using an anytime roadmap planning algorithm called POMP (Pareto-Optimal Motion Planner).

### Disclaimer
Currently the densification strategies use radius of connectivity (where appropriate) to induce subgraphs (check the densification paper for details). The results are obtained from analysis done for unit hypercubes. For other state spaces, please ensure that the distance function is appropriately defined so that the r-disk roadmaps have the appropriate expected average degree. Instead, a k-NN version of our framework could be used (and was started to be implemented in the `kNN_bpomp` branch but was not completed). Contact the author for more details if the r-disk version cannot be used reliably.

# Dependencies

## For the library
- Ubuntu 14.04 or later (Has also been built on OSX)
- [C++11](https://isocpp.org/wiki/faq/cpp11) or higher
- [cmake](https://cmake.org/download/)
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [FLANN](http://www.cs.ubc.ca/research/flann/)
- [OMPL](http://ompl.kavrakilab.org/) built with FLANN
- [Boost Graph Library](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/index.html)
## For the example
- [Networkx](https://networkx.github.io/)
- [OpenCV](http://opencv.org/)

# Build
The `CMakeLists.txt` file supports [catkin tools](https://catkin-tools.readthedocs.io/en/latest/). Once you have created and [initialized your workspace](https://catkin-tools.readthedocs.io/en/latest/quick_start.html#initializing-a-new-workspace), you should be able to clone `batching_pomp` into your workspace and do `catkin build batching_pomp`.

*Please note that the build process has not been extensively tested for various different systems and the authors cannot commit to provide support for the same. Feel free to raise a Github issue though.*

# Usage
`batching_pomp` is implemented as an [OMPL Geometric Multi-Query Planner](http://ompl.kavrakilab.org/planners.html#geometric_planners) and can be invoked on a [`RealVectorStateSpace`](http://ompl.kavrakilab.org/classompl_1_1base_1_1RealVectorStateSpace.html) problem.
It returns a sequence of successively shorter feasible paths on the provided roadmap, for the given motion planning problem. 

The following are the OMPL parameters that the planner supports:

- `increment`(required): The step-wise increment for the value of the `alpha` tradeoff parameter in POMP. `alpha` varies from 0 (only considers edge collision measure) to 1 (only considers edge length)
- `batching_type`(required): The kind of batching method to be used by the `BatchingManager`. Options are vertex, edge, hybrid and single. The first three are defined in paper [1] and the last refers to an explicit roadmap with vertices and edges, that is added in a single batch.
- `roadmap_filename`(required): The path to the `.graphml` file that will be read in by the `BatchingManager`. If the corresponding `batching_type` is single batching, then the vertices should have underlying states and the edges along with their end-point vertices should be present. For any other kind of batching, only the vertices with their states should be present.
- `graph_type`: The kind of sample sequence used to generate the vertex states. Either `halton` or `rgg`. If edge or hybrid batching is chosen, then graph type *must* be provided.
- `prune_threshold`(default=1.05): The threshold of the ratio of previous best solution cost to new best solution cost, above which pruning of existing vertices should be done. This prunes more conservatively.
- `start_goal_radius`(required only for Single Batching): The radius of connectivity to use for connecting the start and goal to the `singlebatching` roadmap
- `selector_type`(optional): Arranges the order of edges to evaluate on a candidate path. If not included, the full candidate path will be evaluated each time, rather than just the first unevaluated edge.

# Example
A `Python` script for generating a complete roadmap has been provided with some options, along with a `C++` example of using the planner. To use them, checkout the `example` branch of the repository. Now we will run `batching_pomp` on a 2D point planning example on a B/W map. All code is assumed to be run from the `batching_pomp` level of the repository. Check the imports in the `scripts/save_sequence_to_graph.py` file to see what `Python` libraries you need. You will also need `OpenCV` to run the example.

The first step is to generate a complete roadmap for the unit 2D grid (the locations in the unit grid will be scaled to the square image, where (0,0) is the upper left and (1,1) is the lower right of the image). To do this, run
```shell
$ python scripts/save_sequence_to_graph.py --num_nodes 1000 --dimensions 2 --bounds unit --sequence halton --outfile data/halton_2d_1000.graphml
```
Open up the `.graphml` file if you have never seen one before and take a look at how the states are associated with the vertex IDs. Note that there are no edges defined between vertex IDs as these are assumed to be complete roadmaps, with the edges based on radii of connectivity, handled by the batching manager.

Now you can run the example in `scripts/example_2d_image.cc` on the map in `data/example_map.png`, once you have built the package and example. The example assumes you have a catkin workspace in which `batching_pomp` resides but you can also just use the typical `cmake` procedure 
```shell
$ catkin build batching_pomp # Assuming you have a catkin workspace
$ ../../devel/lib/batching_pomp/example_2d -f data/halton_2d_5000.graphml -m data/example_map.png -v halton -s 0.1 0.1 -g 0.9 0.9 -b hybrid # There's probably an easier way to call the executable
```
You should see an OpenCV window with a solution path pop up. Keep pressing the 'Enter' key over the window to get shorter and shorter paths.

The above example assumes you have a `catkin` workspace but if you do not, you should still be able to do it through the `cmake` and `make` route after editing the `CMakeLists.txt` file appropriately. Contact the author for help with that if needed.

@def title = "A new array type for domain decomposition in Julia"
@def hascode = true
@def date = "3/15/2022"
@def tags = ["julia", "domain_decomposition"]

# A new array type for domain decomposition in Julia

Part of my research as a grad student was creating a high-performance Fortran code to simulate compressible flow. This code was parallelized using the domain decomposition method, where a large problem is split into many smaller problems that can be solved semi-independently. This is typically done using MPI too facilitate the communication of data between each domain. Neighbor information is exchanged using "halo" cells, which are fictitious grid cells that make communication simple. The tricky part of this is knowing what the neighbor domains are and how to efficiently pass data back and forth.

I have recently converted this Fortran code ([Cato](https://github.com/smillerc/cato)) over to Julia. Initially the parallelization was accomplished using multi-threaded loops. This method of parallelization doesn't scale as well as domain-decomposition however. MPI is available in Julia, but I needed to create a high level library that facilitated halo exchange. [`MPIHaloArrays.jl`](https://github.com/smillerc/MPIHaloArrays.jl) is a new package that provides the `MPIHaloArray` array type, which is a subtype of `AbstractArray`, that does just this. While there are other existing libraries that have similar decomposition functionality, such as [`MPIArrays.jl`](https://github.com/barche/MPIArrays.jl), [`PencilArrays.jl`](https://github.com/jipolanco/PencilArrays.jl), and [`ImplicitGlobalGrid.jl`](https://github.com/eth-cscs/ImplicitGlobalGrid.jl), I wanted to create an array type specifically for this task. Each library only provides only a portion of the functionality I needed, so I decided to try creating a new one. 

Domain-decomposition splits a large problem into small subdomains that live on different MPI ranks. Stencil operations commonly used in solving PDEs require neighbor information, and this is done with halo cells along the edge of each subdomain. In the image below, the domain is decomposed into 4 subdomains, with a single layer of halo cells on all sides. Here process 4 is shown exchanging neighbor information with processes 3 and 2. Domain decomposition can happen in 1D, 2D, and 3D.

![](/assets/images/halo_exchange.png)

## Features

- 1D, 2D, and 3D domain decomposition.
- Neighbor exchange with arbitrarily sized halo cell regions. Note, this is currently fixed across all dimensions. 
- Communication is currently handled by `MPI.ISend` and `MPI.IRecv` underneath, but future versions will give the option for one-sided communication with `MPI.Put` and `MPI.Get`.
- Halo exchange can be orthogonal-only (e.g. `[i+1,j,k]`) , or it can include corners as well (e.g `[i+1,j+1,k]`)

## TODO

- GPU testing – I don’t currently have a local GPU to test this with, but future work will include GPU-aware MPI
- Fix the limitation of 1D, 2D, or 3D arrays. For example, an array `U` could have dimensions for `[q,i,j]`, but halo exchanges are only done on `i` and `j`.
- Optimization of neighbor communication (caching, etc…)
- Finish onesided implementation

## Examples

**Quickstart**

```julia
using MPI, MPIHaloArrays

MPI.Init()
const comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const root = 1

# Create a topology type the facilitates neighbor awareness and global size
topology = CartesianTopology(comm, 
                             [4,2], # domain decompose in a 4x2 grid of ranks
                             [false, false]) # no periodic boundaries

nhalo = 2

# Create some data on the local rank
local_data = rand(10,20) * rank
A = MPIHaloArray(local_data, topology, nhalo)
fillhalo!(A, -1)

# Perform halo exchange
updatehalo!(A)

# Or start with a large global array
B_global = rand(512,512)

# divide by rank and return B::MPIHaloArray
B_local = scatterglobal(B_global, root, nhalo, topology; 
                        do_corners = false) # if your algorithm doesn't need 
                        # corner info, this will save communication
# do some work
# ....
updatehalo!(B_local)
# do some work
# ....

# and pull it back to the root rank
B_result = gatherglobal(B_local; root=root)

GC.gc()
MPI.Finalize()
```

**2D Diffusion**

See in the example [here](https://github.com/smillerc/MPIHaloArrays.jl/blob/main/docs/examples/04-diffusion2d.jl) in the github repository.

### Documentation

See the docs [here](https://smillerc.github.io/MPIHaloArrays.jl/stable/)
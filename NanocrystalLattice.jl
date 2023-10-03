# This file contains the functions for describing the lattice of nanocrystals.
# The details of each function can be found in the comment of each function below.

# ------------------------------------------

"""
Function that returns the indices of the nanocrystal in a orthorhombic lattice
Indexing goes along x-axis -> y-axis -> z-axis
```julia
julia> NanocrystalIndexing(10, 10, 1)
```
## INPUT
`nx`: number of nanocrystal along x direction 
`ny`: ____________//_____________ y direction
`nz`: ____________//_____________ z direction
"""
function NanocrystalIndexing(nx, ny, nz)
  result = zeros(Int32, nx, ny, nz)
  for inx=1:nx
   for iny=1:ny
    for inz=1:nz
     result[inx,iny,inz] = inx+nx*(iny-1)+nx*ny*(inz-1)
    end
   end
  end
  return result
end

# ------------------------------------------

"""
Function that returns the 3D position vectors of nanocrystals in a cubic superlattice
corresponding to their indices, which are computed from function "NanocrystalIndexing"
```julia
julia> NanocrystalIndexing(10, 10, 1)
```
## INPUT
`nx`: number of nanocrystal along x direction 
`ny`: ____________//_____________ y direction
`nz`: ____________//_____________ z direction
`length`:  distance between the centers of two adjacent nanocrystals
           each nanocrystal encapsulated by ligand is assumed to be a cube
`NCIndex`: 3D array of dimension nx*ny*nz containing indices of the nanocrystals
"""
function NanocrystalPositions(nx, ny, nz, length, NCIndex)
  tot_num_NC = nx*ny*nz
  result = zeros(3, tot_num_NC)
  for inx=1:nx
   for iny=1:ny
    for inz=1:nz
     current_index = NCIndex[inx,iny,inz]
     result[1:3, current_index] = length * [inx-0.5; iny-0.5; inz-0.5]
    end
   end
  end
  return result
end




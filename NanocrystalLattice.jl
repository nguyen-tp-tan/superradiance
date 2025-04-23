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
julia> NanocrystalPositions(2, 2, 1, 2.5, [1 2 3; 4 5 6])
```
## INPUT
`nx`: number of nanocrystal along x direction 

`ny`: ____________//_____________ y direction

`nz`: ____________//_____________ z direction

`length`: distance between the centers of two adjacent nanocrystals
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

# ------------------------------------------

"""
Function that returns the 3D position vectors of nanocrystals in a cubic superlattice
corresponding to their indices, which are computed from function "NanocrystalIndexing"
```julia
julia> NanocrystalPositionsOrthorhombic(3, 2, 1, [2.5, 4.0, 2.5], [1 2 3; 4 5 6])
```
## INPUT
`nx`: number of nanocrystal along x direction 

`ny`: ____________//_____________ y direction

`nz`: ____________//_____________ z direction

`basis_lengths`: distances between the centers of two adjacent nanocrystals along 3 axes
                 each nanocrystal encapsulated by ligand is assumed to be a cube

`NCIndex`: 3D array of dimension nx*ny*nz containing indices of the nanocrystals
"""
function NanocrystalPositionsOrthorhombic(nx, ny, nz, basis_lengths, NCIndex)
  tot_num_NC = nx*ny*nz
  result = zeros(3, tot_num_NC)
  if (size(basis_lengths)[1] != 3)
    return result
  end

  for inx=1:nx
    for iny=1:ny
      for inz=1:nz
        current_index = NCIndex[inx,iny,inz]
        result[1:3, current_index] = 
          [basis_lengths[1] * inx; basis_lengths[2]*iny; basis_lengths[3]*inz]
      end
    end
  end
  return result
end

# ------------------------------------------
"""
Function that returns the 3D positions of a 2D array of nanocrystals
in xy-plane with a rectangular empty region in the middle
## INPUT
`nx` : max number of nanocrystals along x direction
`ny` : _______________//_______________ y direction
`percent_empty` : percentage of the empty area
`dx_c2c` : center-to-center distance along x direction
`dy_c2c` : _______________//______________
`z_coor` : z coordinates of the plane, zCoor = 0 if the array lies in xy-plane
"""
function NanocrystalPositions2DRectangularWithHole(
  nx::Int, ny::Int, dx_c2c::Float64, dy_c2c::Float64, z_coor::Float64, 
  hole_ratio_x::Float64, hole_ratio_y::Float64)
  result = zeros(3, 0)
  if (hole_ratio_x > 1 || hole_ratio_y > 1 || (hole_ratio_x == 1 && hole_ratio_y == 1))
    return result
  end

  is_empty = (i::Int, n::Int, percent::Float64) -> begin
    upper_limit = 1 + 0.5*(n-1)*(1+percent)
    lower_limit = 1 + 0.5*(n-1)*(1-percent)
    return (i >= lower_limit) && (i <= upper_limit)
  end

  for iy=1:ny
    for ix=1:nx
      if (is_empty(ix, nx, hole_ratio_x) && is_empty(iy, ny, hole_ratio_y))
        continue;
      end
      current_position = [ix*dx_c2c; iy*dy_c2c; z_coor]
      result = [result current_position]
    end
  end
  return result
end

# ------------------------------------------
"""
Function that returns a matrix of positions (in 3D) of 2D array of nanocrystals
of which its plane is parallel to xy-plane and the a_1 primitive vector points 
along the x direction
## INPUT
`a` : length of the primitive translational vectors
`n1` : number of points along a line parallel to the a_1 primitive vector
`n2` : number of points along a line parallel to the a_2 primitive vector
`z_coor` : z coordinates of the plane, zCoor = 0 if the array lies in xy-plane
"""
function NanocrystalPositions2DHexagonal(a, n1, n2, z_coor)
  HALF_FRACION = 1//2
  SQUARE_ROOT_3 = sqrt(3.0)
  basis_a1 = [1, 0]
  basis_a2 = [-HALF_FRACION, HALF_FRACION]

  hex_coordinates = zeros(3, n1*n2)
  hex_coordinates[3, :] .= z_coor
  for i2=1:n2
    for i1=1:n1
      idx_2d = (i2-1)*n1 + i1
      hex_coordinates[1, idx_2d] = mod((i1-1)*basis_a1[1] + (i2-1)*basis_a2[1], n1) * a
      hex_coordinates[2, idx_2d] = ((i1-1)*basis_a1[2] + SQUARE_ROOT_3*((i2-1)*basis_a2[2])) * a
    end
  end
  return hex_coordinates
end

# ------------------------------------------
"""
Function that returns a matrix of positions (in 3D) of nanocrystals arranged 
in a 1-dimensional ring contained in the xy-plane
## INPUT
`N` : total number of nanocrystals in the ring
`LNC` : average size of nanocrystals
`Lshell` : shell length
"""
function NanocrystalPositions1DRing(n, length, z_coor)
  ring_pos = zeros(Float64, 3, n)
  ring_pos[3,:] .= z_coor

  ring_rad = (n//2) *BigFloat(length/pi) 
  angles = BigFloat(pi) * 2*range(0,n-1)//n
  
  ring_pos[1,:], ring_pos[2,:] = ring_rad*cos.(angles), ring_rad*sin.(angles)

  if (n%2 == 0)
    ring_pos[2, n ÷ 2 + 1] = 0.0
  end

  if (n%4 == 0)
    ring_pos[1, n ÷ 4 + 1] = 0.0
    ring_pos[1, 3 * n ÷ 4 + 1] = 0.0
  end

  return ring_pos
end


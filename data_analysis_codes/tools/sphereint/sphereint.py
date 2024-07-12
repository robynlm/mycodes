"""This module provides the SphereIntegrate class.

Copyright (C) 2022  Robyn L. Munoz

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

You may contact the author at : robyn.munoz@yahoo.fr
"""

import numpy as np
import scipy as sc
import itertools 

class SphereIntegrate:
    def __init__(self, N, L, centre):
        """Initialises grid parameters.

        This class provides a weight for each grid position based on 
        whether or not it is in, out, or partially in a sphere of a 
        given radius (with ) or in a given domain (with ) which is 
        delimited by a sign change (only one sign change 
    
        A cubic cell is placed around each grid position and the volume of 
        the cell in the sphere (assuming a flat suface in the cell) 
        is calculated and normalised by the cell volume to obtain the weight.
        
        Parameters:
        N : int, number of data points in each direction of the grid
        L : float, size of grid
            Here I assume the data box is a cube that 
            ranges from -L/2 to L/2.
        centre : (3) list of integers
                 Indexes of the grid position where the centre 
                 of the sphere or domain is placed.

        Attributes:
        dx : float, grid spacing
        volume_cell : float, volume of one grid cell
        area_cell : float, area of one side of a grid cell
        radius : float, radius of sphere to calculate the volume or area 
                 of with the get_box_weights_radius function.
                 The user should feel free to overwrite this.
        data_function : function, used to interpolate boundary points.
                        It is intially defined as the radius_function 
                        to be used to integrate in/on a sphere. 
                        To integrate in/on a given domain, this is overwritten 
                        with scipy.interpolate.RegularGridInterpolator.
                        The user shouldn't bother with this.
        
        Notes:
        - The input grid is assumed to be periodic.
        """

        # input
        self.N = N
        self.L = L
        self.centre = centre
        
        self.dx = L / N
        self.x1d = np.arange(-self.L/2, self.L/2, self.dx)
        self.volume_cell = self.dx**3
        self.area_cell = self.dx**2
        
        self.radius_max = 0.0
        self.data_function = self.radius_function
        self.center_sign = +1
        
        # nbr of points considered in Monte Carlo
        self.Nrand_volume = 10000 # when calculating the volume
        self.Nrand_area = 10 # when calculating the area
        self.geometrical = True 
        # if False, Monte Carlo method is used exclusively
        
        self.verbose = False
        self.nbr_cas_de_merde = 0 
    
    def get_box_weights_radius(self, domain='volume'):
        """Compute weight at every grid position for given radius.
        
        Parameters : 
            domain : either 'volume' (default) or 'area' 
        
        Returns :
            (N, N, N) array_like of floats of weight value            
        """
        
        
        if self.radius_max > self.L / 2:
            # Radius can't be bigger than data box.
            # Because I use periodic boundary conditions 
            # the sphere would fold on itself.
            print('ERROR: that radius is too big for me sorry')
        else:
            x, y, z = np.meshgrid(self.x1d, self.x1d, self.x1d, indexing='ij')
            data = self.shift_grid0_to_igrid(
                self.radius_function_grid(x, y, z))
            
            weight = self.get_box_weights_data(
                data, 
                data_function = self.radius_function, 
                domain = domain)
            
            return weight
   
    def get_box_weights_data(self, data, data_function = None, domain='volume'):
        """Compute weight at every grid position for given boundary.
        
        Parameters : 
            domain : either 'volume' (default) or 'area' 
        
        Returns :
            (N, N, N) array_like of floats of weight value            
        """
        
        # shift, input grid to grid around 0
        # provided center location -> 0
        data = self.shift_igrid_to_grid0(data)
        
        # Define function representing the data
        if data_function is None:
            self.data_function = sc.interpolate.RegularGridInterpolator(
                (self.x1d, self.x1d, self.x1d), 
                data, method='cubic')
        else:
            self.data_function = data_function
        
        # Cells inside the domain
        self.center_sign = np.sign(data[
            int(self.N/2), int(self.N/2), int(self.N/2)]) 
        cell_in_domain = np.sign(data) == self.center_sign
        
        # Am I calculating the volume or area?
        weight = np.zeros((self.N, self.N, self.N))
        if domain=='volume':
            weight[cell_in_domain] = 1.0
            cell_weight = self.volume_cell
            if self.geometrical:
                compute_weight = self.compute_volumes_geometrical
            else:
                compute_weight = self.compute_volumes_MonteCarlo
        else:
            cell_weight = self.area_cell
            if self.geometrical:
                compute_weight = self.compute_areas_geometrical
            else:
                compute_weight = self.compute_areas_MonteCarlo
        
        # initialise counter
        self.nbr_cas_de_merde = 0
           
        # ====================================
        # Calculate boundary cells
         
        # Create grid of vertices of the 
        # cells containing the data points
        x1dnodes = self.x1d[:-1] + (self.dx / 2)
        xn, yn, zn = np.array(np.meshgrid(
            x1dnodes, x1dnodes, x1dnodes, indexing='ij'))
        
        # Interpolate to get data value
        xyznodes = np.array([xn.ravel(), yn.ravel(), zn.ravel()]).T
        data_nodes = self.data_function(xyznodes).reshape(
            (self.N-1), (self.N-1), (self.N-1))

        # Difference in sign in various directions
        xdiff = abs(np.diff(np.sign(data_nodes), axis=0))
        ydiff = abs(np.diff(np.sign(data_nodes), axis=1))
        zdiff = abs(np.diff(np.sign(data_nodes), axis=2))
        # Combine to see if data point cell experiences a sign change
        signchange = (xdiff[:, :-1, :-1] + xdiff[:, 1:, 1:] 
                      + ydiff[:-1, :, :-1] + ydiff[1:, :, 1:]
                      + zdiff[:-1, :-1, :] + zdiff[1:, 1:, :]) 
        # Introduce a buffer 
        signchange += np.append(
            signchange[1:, :, :], 
            np.array([signchange[0, :, :]]), 
            axis=0)
        signchange += np.append(
            signchange[:, 1:, :], 
            np.array([signchange[:, 0, :]]).reshape(
                self.N - 2, 1, self.N - 2), 
            axis=1)
        signchange += np.append(
            signchange[:, :, 1:], 
            np.array([signchange[:, :, 0]]).reshape(
                self.N - 2, self.N - 2, 1), 
            axis=2)
        # Cell boundary domain
        cell_boundary = signchange != 0    
       
        # ====================================
        # Loop over boundary cells
        
        counter = 0
        total_count = np.sum(cell_boundary)
        if self.verbose:
            print(total_count, ' boundary cells') 
        for ix, iy, iz in zip(*np.where(cell_boundary)):
            # vertices indices are offset by 1 compared to data indices
            ix += 1
            iy += 1
            iz += 1
            # data point location
            xyz = np.array([self.x1d[ix], self.x1d[iy], self.x1d[iz]])
            # coordinates of its cell's corners (or vertices)
            points = np.array(self.cell_corner_positions(xyz))
            # calculate weight
            
            # number of points in the domain
            corner_sign = np.sign(self.data_function(points))
            nbr_points_in_sphere = np.sum([psign==self.center_sign 
                                           for psign in corner_sign])
            
            if 0 < nbr_points_in_sphere and nbr_points_in_sphere < 8:
                weight[ix,iy,iz] = compute_weight(points) / cell_weight
            elif nbr_points_in_sphere == 8 and domain=='volume':
                weight[ix, iy, iz] = 1
            
            counter += 1
            if self.verbose and (counter % int(total_count / 50) == 0):
                print('percent done : ', counter * 100 / total_count )
                
        # shift back, grid around 0 to input grid
        # 0 -> provided center location
        weight = self.shift_grid0_to_igrid(weight)
        
        if self.verbose:
            print(self.nbr_cas_de_merde, ' cas de merde')
        return weight

    ########################################################################
    ########################################################################
    ########################################################################
    #            Helpers
    ########################################################################
    ########################################################################
    ########################################################################
        
    def radius_function(self, coords):
        """If coordinates are inside (negative) or outside (positive) a sphere.

        Parameters:
        coord (list of tuples or array-like): A list or array of coordinates,
            where each coordinate is a tuple or list of (x, y, z) cartesian
            values on the grid around 0.

        Returns:
        np.ndarray: An array of condition values, where the value is negative
            if the radius is smaller than `self.maximum_radius` and positive
            otherwise.
            
        Note:
        - This funtion's input and output formats are the same as
            scipy.interpolate.RegularGridInterpolator.
        """
        return np.sqrt(np.sum(coords**2, axis=1)) - self.radius_max
    
    def radius_function_grid(self, x, y, z):
        return np.sqrt(x**2 + y**2 + z**2) - self.radius_max 
    
    def radius(self, points):
        if type(points)==list:
            points = np.array(points)
        return np.sqrt(np.sum(points**2))  
       
    def shift_grid0_to_igrid(self, data):
        """Shift values to be around wanted grid center.
        
        Parameters :
            phi : (N, N, N) array_like
                  Value to be shifted.
        
        Returns :
            (N, N, N) array_like
        """
        xshift = int(self.N/2) - self.centre[0]
        yshift = int(self.N/2) - self.centre[1]
        zshift = int(self.N/2) - self.centre[2]
        data = np.append(data[xshift:, :, :], data[:xshift, :, :], axis=0)
        data = np.append(data[:, yshift:, :], data[:, :yshift, :], axis=1)
        data = np.append(data[:, :, zshift:], data[:, :, :zshift], axis=2)
        return data
        
    def shift_igrid_to_grid0(self, data):
        """Shifts a 3D np.array around the zero in a periodic grid."""
        xshift = int(self.N/2) - self.centre[0]
        yshift = int(self.N/2) - self.centre[1]
        zshift = int(self.N/2) - self.centre[2]
        data = np.append(data[-xshift:,:,:], data[:-xshift,:,:], axis=0)
        data = np.append(data[:,-yshift:,:], data[:,:-yshift,:], axis=1)
        data = np.append(data[:,:,-zshift:], data[:,:,:-zshift], axis=2)
        return data
        
    def cell_corner_positions(self, xyz):
        """ Provide cell corner positions.
        
        A cubic cell of size dx is placed around the grid point,
        the coordinate position of the corners (or vertices)
        are provided here.
        
        Parameters :
            xyz : (3) array_like
                  Coordinate position of grid point.
            
        Returns :
            (8) list
            Each element is a (3) array_like coordinate position
               3--------1
              /|       /|
             2--------0 |
             | |      | |
             | 5------|-7
             |/       |/
             4--------6
        """
        dxmax = self.dx / 2
        x, y, z = xyz[0], xyz[1], xyz[2]
        return [np.array([x + dxmax, y + dxmax, z + dxmax]),
                np.array([x - dxmax, y + dxmax, z + dxmax]),
                np.array([x + dxmax, y - dxmax, z + dxmax]),
                np.array([x - dxmax, y - dxmax, z + dxmax]),
                np.array([x + dxmax, y - dxmax, z - dxmax]),
                np.array([x - dxmax, y - dxmax, z - dxmax]),
                np.array([x + dxmax, y + dxmax, z - dxmax]),
                np.array([x - dxmax, y + dxmax, z - dxmax])]
    
    def interpolate_get_distance(self, point, neighbourpoint):
        """Compute distance between cell corner and sphere boundary on cell edge
        
        1) Compute intersection between sphere and direction passing 
        the cell corner in the sphere and its' neighbour 
        that is outside the sphere. 
        Each neighbouring point share two coordinate values with 
        the cell corner in the sphere and the third coordinate 
        is + or - dx different. 
        The intersecting point shares the first two coordinate values, 
        but the last one needs to be computed.
        
        2) Then compute the distance between the cell corner 
        and the intersecting point.
        
        Parameters :
            points : (3) array_like 
                     Coordinate position of the cell corner.
            neighbourpoint : (?) list
                             Each element is a (3) array_like coordinate 
                             position of the cell corner neighbours 
                             that are outside of the sphere, 
                             dimension can go from 1 to 3.
                             
        Returns :
            list : depths, ixyz
                   depths : (?) array_like
                            Distances between the cell corner in the sphere
                            and the point intersecting the sphere 
                            and the neighbouring edge, 
                            dimension can go from 1 to 3.
                   ixyz : (?) list
                          Direction along which intersecting point lies, 
                          with 0 -> x, 1 -> y, 2 -> z, 
                          dimension can go from 1 to 3.
        """
        distances_to_point = []
        ixyz = []
        fA = self.data_function(np.array([point]))[0]
        for p in neighbourpoint:
            # Two of the coordinates are the same
            interpolated_point = p.copy()
            
            # Index of coordinate that needs to change
            ctochange = np.where(point - p != 0)
            ixyz += [ctochange]
            
            # New coordinate
            new_coord = self.intersect(
                point[ctochange], fA, 
                p[ctochange] , self.data_function(np.array([p]))[0])
            
            interpolated_point[ctochange] = new_coord
            
            # Distance between interpolated point and cell corner in sphere
            d_to_point = np.sqrt(np.sum( (point - interpolated_point)**2 ))
            distances_to_point += [d_to_point]
            
            # Check that distance is smaller than cell size
            if self.verbose and (d_to_point - self.dx)/self.dx > 1e-14:
                print('WARNING: (depth - dx) / dx > 1e-14, dx =', 
                      self.dx, 'and depth =', d_to_point)
        return np.array(distances_to_point), ixyz
    
    def intersect(self, xA, fA, xB, fB):
        slope = (fB - fA) / (xB - xA)
        ordonnee = fA - slope * xA
        return - ordonnee / slope 
    
    def intersect_3D(self, pointA, pointB):
        nbr_overlap = np.sum(pointA == pointB)
        # if these two points are the same
        if nbr_overlap == 3:
            return None
        else:
            fA = self.data_function(np.array([pointA]))[0]
            fB = self.data_function(np.array([pointB]))[0]
            # if there is a sign change betweent these points
            if np.sign(fA)!=np.sign(fB):
                # depending on if they have coordinates in common
                if nbr_overlap == 0:
                    # intersect in a volume
                    # Not doing that, I'm using scipy minimize instead
                    interpolated_point = None
                elif nbr_overlap == 1:
                    # intersect on a plane
                    ctochange = np.where(pointA - pointB != 0)
                    rA = self.radius(pointA[ctochange])
                    rB = self.radius(pointB[ctochange])
                    new_r = self.intersect(rA, fA, rB, fB)
                    interpolated_point = pointA.copy()
                    interpolated_point[ctochange] = pointA[ctochange] + (
                        (pointB[ctochange] - pointA[ctochange]) 
                        * (new_r - rA) / (rB - rA) )
                elif nbr_overlap == 2:
                    # intersect on a line
                    ctochange = np.where(pointA - pointB != 0)
                    interpolated_point = pointA.copy()
                    interpolated_point[ctochange] = self.intersect(
                        pointA[ctochange], fA, 
                        pointB[ctochange] , fB)
                return interpolated_point
            else:
                return None
    
    def data_function_minimize(self, coord):
        return np.abs(self.data_function(np.array([coord]))[0])
        
    def triangle_area_from_points(self, pointA, pointB, pointC):
        triangle_sides = np.array([self.radius(pointA - pointB),
                                   self.radius(pointB - pointC),
                                   self.radius(pointC - pointA)])
        return self.triangle_area(triangle_sides)

    def triangle_area_from_depth(self, depth):
        triangle_sides = np.sqrt(np.array([depth[0]**2 + depth[1]**2,
                                           depth[0]**2 + depth[2]**2,
                                           depth[1]**2 + depth[2]**2]))
        
        return self.triangle_area(triangle_sides)
    
    def triangle_area(self, triangle_sides):
        # Heron's formula
        semi_perimeter = np.sum(triangle_sides)/2
        area = np.sqrt(
            np.max([0, semi_perimeter
                       * np.prod(semi_perimeter - triangle_sides)]))
        return area 

    ########################################################################
    ########################################################################
    ########################################################################
    #            Volume
    ########################################################################
    ########################################################################
    ########################################################################
    
    def compute_volumes_MonteCarlo(self, points):
        xvals = list(points[:,0])
        yvals = list(points[:,1])
        zvals = list(points[:,2])
        random_points = np.array(
            [np.random.uniform(low=np.min(xvals), high=np.max(xvals), 
                               size=self.Nrand_volume), 
             np.random.uniform(low=np.min(yvals), high=np.max(yvals), 
                               size=self.Nrand_volume), 
             np.random.uniform(low=np.min(zvals), high=np.max(zvals), 
                               size=self.Nrand_volume)]).T
        fraction_in_domain = np.sum(np.sign(
            self.data_function(random_points)
            ) == self.center_sign) / self.Nrand_volume
        volume = fraction_in_domain * self.volume_cell
        return volume
        
    def compute_volumes_geometrical(self, points):
        """Compute volume of the cell contained in the sphere.
        
        Parameters :
            points : (8) list
                     Each element is a (3) array_like coordinate 
                     position of the cell corner.
            points_in_sphere : (8) list
                               Each element is a boolean:
                                - True : cell corner in sphere
                                - False : cell corner not in sphere
            
        Returns :
            float, volume of the cell contained in the sphere
        """
        corner_sign = np.sign(self.data_function(points))
        points_in_sphere = [psign==self.center_sign
                            for psign in corner_sign]
                
        nbr_points_in_sphere = np.sum(points_in_sphere)
        if nbr_points_in_sphere == 8:
            volume = self.volume_cell
        elif nbr_points_in_sphere == 0:
            volume = 0.0
        else:
            points = np.array(points)
            # Identify the neighbouring corner point of each corner.
            # For example, the corner : points[0] (x + dx/2, y + dx/2, z + dx/2)
            # Has the neighbours: 
            # points[1] (x - dx/2, y + dx/2, z + dx/2), neighbour along x
            # points[2] (x + dx/2, y - dx/2, z + dx/2), neighbour along y
            # points[6] (x + dx/2, y + dx/2, z - dx/2), neighbour along z
            neighbouring_points = np.array([[points[1], points[2], points[6]], 
                                            [points[0], points[3], points[7]], 
                                            [points[0], points[3], points[4]], 
                                            [points[1], points[2], points[5]], 
                                            [points[2], points[5], points[6]], 
                                            [points[3], points[4], points[7]], 
                                            [points[0], points[4], points[7]], 
                                            [points[1], points[5], points[6]]])
            
            sphere_mask = np.where(points_in_sphere)
            not_sphere_mask = np.where(~np.array(points_in_sphere))
            
            # Compute volume.
            # The radius, corner positions and neighbouring points (masked depending
            # on the sphere) are passed to the volume_*_points function 
            # according to the number of corners that are in the sphere.
            if nbr_points_in_sphere == 7:
                volume = self.volume_1_point(
                    points[not_sphere_mask], neighbouring_points[not_sphere_mask])
                volume = self.volume_cell - volume
            elif nbr_points_in_sphere == 6:
                volume = self.volume_2_points(
                    points[not_sphere_mask], neighbouring_points[not_sphere_mask])
                volume = self.volume_cell - volume
            elif nbr_points_in_sphere == 5:
                volume = self.volume_3_points(
                    points[not_sphere_mask], neighbouring_points[not_sphere_mask])
                volume = self.volume_cell - volume
            elif nbr_points_in_sphere == 4:
                volume = self.volume_4_points(
                    points[sphere_mask], neighbouring_points[sphere_mask])
            elif nbr_points_in_sphere == 3:
                volume = self.volume_3_points(
                    points[sphere_mask], neighbouring_points[sphere_mask])
            elif nbr_points_in_sphere == 2:
                volume = self.volume_2_points(
                    points[sphere_mask], neighbouring_points[sphere_mask])
            else: # nbr_points_in_sphere == 1
                volume = self.volume_1_point(
                    points[sphere_mask], neighbouring_points[sphere_mask])
        return volume
    
    def cas_de_merde_volume(self, points, pointneighbours):
        self.nbr_cas_de_merde += 1
        xvals = list(points[:,0]) + list(itertools.chain(*pointneighbours[:,:,0]))
        yvals = list(points[:,1]) + list(itertools.chain(*pointneighbours[:,:,1]))
        zvals = list(points[:,2]) + list(itertools.chain(*pointneighbours[:,:,2]))
        points = np.array([xvals, yvals, zvals]).T
        return self.compute_volumes_MonteCarlo(points)
    
    def volume_1_point(self, point, pointneighbour):  
        """Compute cell volume in the sphere when 1 corner is in the sphere.
        
        Compute the volume of a trirectangular tetrahedron.
        
        Parameters :
            point : (3) array_like 
                    Coordinate position of the cell corner.
            neighbourpoint : (3) array_like
                             Each element is a (3) array_like coordinate 
                             position of the cell corner neighbours.  
            
        Returns :
            float, volume            
        """
        depth, ixyz = self.interpolate_get_distance(
            point[0], pointneighbour[0])
        volume = np.prod(depth) / 6
        volume = np.min([volume, self.volume_cell / 6])
        return volume
        
        
    def volume_2_points(self, points, pointneighbours):  
        """Compute cell volume in the sphere when 2 corners are in the sphere.
        
        Compute the volume of a trirectangular tetrahedron 
        that extends larger than the cell size in one direction. 
        Then remove that extension that correspondes to 
        a smaller trirectangular tetrahedron
        such that we only consider the part in the cell.
        To find the side that needs to be extended, 
        the area of each triangular base, the smallest one is extended.
        If the two areas are equal then we have a right triangular prism.
        
        Parameters :
            point : (2, 3) array_like 
                    Coordinate positions of the 2 cell corners.
            neighbourpoint : (2, 3) array_like
                             Coordinate positions of the cell corner neighbours.  
            
        Returns :
            float, volume     
        """
        # Cell corners in the sphere
        point1 = points[0]
        point2 = points[1]
        
        # Neighbouring points that are outside the sphere
        pointneighbour1 = [neighbour 
                           for neighbour in pointneighbours[0] 
                           if list(neighbour)!=list(point2)]
        pointneighbour2 = [neighbour 
                           for neighbour in pointneighbours[1] 
                           if list(neighbour)!=list(point1)]
        
        # Distances between corner cell in sphere and sphere boundary.
        depth1, ixyz1 = self.interpolate_get_distance(
            point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_get_distance(
            point2, pointneighbour2)
         
        len_depth = np.sort([len(d) for d in [depth1, depth2]])
        if np.sum(len_depth == [2, 2]) == len(len_depth):
            # Area of each base
            area1 = np.prod(depth1) / 2
            area2 = np.prod(depth2) / 2
            
            if abs(area1 / area2 - 1)>1e-10: 
                # Base triangles different so trirectangular tetrahedron considered.
                # To find this: volume = big_tetrahedron - small_tetrahedron
                #                      = entire_shape - extended_part
                # I'm making sure the depth/width of the triangles overlap
                # Create dict with {coord_direction: depth_val}
                dict1 = {ixyz1[i][0][0]:depth1[i] 
                        for i in range(len(ixyz1))}
                dict2 = {ixyz2[i][0][0]:depth2[i] 
                        for i in range(len(ixyz2))}
                idict = list(dict1.keys())
                if area1 > area2: # Point1 has the bigger base
                    big_width = dict1[idict[0]]
                    big_depth = dict1[idict[1]]
                    small_width = dict2[idict[0]]
                    small_depth = dict2[idict[1]]
                else: # Point2 has the bigger base
                    big_width = dict2[idict[0]]
                    big_depth = dict2[idict[1]]
                    small_width = dict1[idict[0]]
                    small_depth = dict1[idict[1]]
                    
                extended_height = (small_width * self.dx 
                                / ( big_width - small_width ))
                
                vol_big = ( (self.dx + extended_height) 
                        * big_width * big_depth) / 6
                vol_small = ( extended_height * small_width * small_depth) / 6
                volume = vol_big - vol_small
            else:
                # This is a right triangular prism
                volume = self.dx * np.average([area1, area2])
                
            volume = np.min([volume, self.volume_cell / 2])
        else:
            volume = self.cas_de_merde_volume(points, pointneighbours)
        return volume
    
    def volume_3_points(self, points, pointneighbours):
        """Compute cell volume in the sphere when 3 corners are in the sphere.
        
        Compute the volume of a trirectangular tetrahedron 
        that extends larger than the cell size in two directions. 
        Then remove those extensions that correspond to 
        smaller trirectangular tetrahedrons
        such that we only consider the part in the cell.
        
        Parameters :
            point : (3, 3) array_like 
                    Coordinate positions of the 3 cell corners.
            neighbourpoint : (3, 3) array_like
                             Coordinate positions of the cell corner neighbours.  
            
        Returns :
            float, volume     
        """
        # Cell corners in the sphere
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        
        # Neighbouring points that are outside the sphere
        pointneighbour1 = [neighbour 
                           for neighbour in pointneighbours[0] 
                           if (list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point3))]
        pointneighbour2 = [neighbour 
                           for neighbour in pointneighbours[1] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point3))]
        pointneighbour3 = [neighbour 
                           for neighbour in pointneighbours[2] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point2))]
        
        # Distances between corner cell in sphere and sphere boundary.
        depth1, ixyz1 = self.interpolate_get_distance(
            point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_get_distance(
            point2, pointneighbour2)
        depth3, ixyz3 = self.interpolate_get_distance(
            point3, pointneighbour3)
        depth = [depth1, depth2, depth3]
        ixyz = [ixyz1, ixyz2, ixyz3]
       
        len_depth = np.sort([len(d) for d in depth])
        if np.sum(len_depth == [1, 2, 2]) == len(len_depth): 
            # ID points
            # Point with only one depth value is A
            iA = np.where(np.array(list(map(len, depth)))==1)[0][0]
            depthA = depth[iA]
            # The two other points are B and C 
            # they will each have a tetrahedron extended through their side.
            iB, iC = np.delete(np.arange(3), iA)
            ixyzB = ixyz[iB]
            ixyzC = ixyz[iC]
            depthB = depth[iB]
            depthC = depth[iC]
                
            # ID depths
            A_width = depthA[0]
            # {coord_direction: depth_val}
            dictB = {ixyzB[i][0][0]:depthB[i] for i in range(len(ixyzB))}
            dictC = {ixyzC[i][0][0]:depthC[i] for i in range(len(ixyzC))}
            # All corners have a width
            key_width = list(dictB.keys() & dictC.keys())[0]
            B_width = dictB[key_width]
            C_width = dictC[key_width]
            # Remaining coordinates correspond to depth and height
            B_depth = [dictB[key] for key in dictB if key!=key_width][0]
            C_height = [dictC[key] for key in dictC if key!=key_width][0]
            
            # Compute extended part
            B_height = np.max([B_width * self.dx / (A_width - B_width), 0])
            C_depth = np.max([C_width * self.dx / (A_width - C_width), 0])
                
            # Compute tetrahedron volume
            vol_B = B_height * B_width * B_depth / 6
            vol_C = C_height * C_width * C_depth / 6
            vol_big = (self.dx + B_height) * A_width * (self.dx + C_depth) / 6
            volume = vol_big - vol_B - vol_C
            
            volume_max = self.volume_cell * 4 / 6
            rel_diff = (volume - volume_max)/volume_max
            if rel_diff > 1e-14:
                volume = self.cas_de_merde_volume(points, pointneighbours)
        else:
            volume = self.cas_de_merde_volume(points, pointneighbours)
        return volume
    
    def volume_4_points(self, points, pointneighbours):
        """Compute cell volume in the sphere when 4 corners are in the sphere.
        
        Two cases :
            
        1 ) A truncated right square prism
        
        2 ) Compute the volume of a trirectangular tetrahedron 
            that extends larger than the cell size in three directions. 
            Then remove those extensions that correspond to 
            smaller trirectangular tetrahedrons
            such that we only consider the part in the cell.
        
        Parameters :
            point : (4, 3) array_like 
                    Coordinate positions of the 3 cell corners.
            neighbourpoint : (4, 3) array_like
                             Coordinate positions of the cell corner neighbours.  
            
        Returns :
            float, volume     
        """
        # Cell corners in the sphere
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        point4 = points[3]
        
        # Neighbouring points that are outside the sphere
        pointneighbour1 = [neighbour 
                           for neighbour in pointneighbours[0] 
                           if (list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point3) 
                               and list(neighbour)!=list(point4))]
        pointneighbour2 = [neighbour 
                           for neighbour in pointneighbours[1] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point3) 
                               and list(neighbour)!=list(point4))]
        pointneighbour3 = [neighbour 
                           for neighbour in pointneighbours[2] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point4))]
        pointneighbour4 = [neighbour 
                           for neighbour in pointneighbours[3] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point3))]
        
        # Distances between corner cell in sphere and sphere boundary.
        depth1, ixyz1 = self.interpolate_get_distance(
            point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_get_distance(
            point2, pointneighbour2)
        depth3, ixyz3 = self.interpolate_get_distance(
            point3, pointneighbour3)
        depth4, ixyz4 = self.interpolate_get_distance(
            point4, pointneighbour4)
        depth = [depth1, depth2, depth3, depth4]
        
        len_depth = np.sort([len(d) for d in depth])
        if sum(len_depth) == 4:
            # A truncated right square prism
            volume = ( np.min([np.max(depth), self.dx]) 
                      + np.min([np.min(depth), self.dx]) ) * self.area_cell / 2
        elif np.sum(len_depth == [0, 2, 2, 2]) == len(len_depth):  
            # A trirectangular tetrahedron
            # ID points
            # Point whose neighbours are all in the sphere is A
            iA = np.where(np.array(list(map(len, depth)))==0)[0][0]
            # The two other points are B, C and D
            # they will each have a tetrahedron extended through their side.
            ixyz = [ixyz1, ixyz2, ixyz3, ixyz4]
            iB, iC, iD = np.delete(np.arange(4), iA)
            ixyzB = ixyz[iB]
            ixyzC = ixyz[iC]
            ixyzD = ixyz[iD]
            depthB = depth[iB]
            depthC = depth[iC]
            depthD = depth[iD]
            
            # ID depths
            # {coord_direction: depth_val}
            dictB = {ixyzB[i][0][0]:depthB[i] for i in range(len(ixyzB))}
            dictC = {ixyzC[i][0][0]:depthC[i] for i in range(len(ixyzC))}
            dictD = {ixyzD[i][0][0]:depthD[i] for i in range(len(ixyzD))}
            # ID the depths
            key_width = list(dictB.keys() & dictC.keys())[0]
            B_width = dictB[key_width] 
            C_width = dictC[key_width] 
            key_depth = list(dictB.keys() & dictD.keys())[0]
            B_depth = dictB[key_depth] 
            D_depth = dictD[key_depth] 
            key_height = list(dictC.keys() & dictD.keys())[0]
            C_height = dictC[key_height] 
            D_height = dictD[key_height] 
            
            # Compute extended part
            B_height = B_depth * (self.dx - C_height) / (self.dx - B_depth)
            C_depth = C_width * (self.dx - D_depth) / (self.dx - C_width)
            D_width = D_height * (self.dx - B_width) / (self.dx - D_height)
            
            # Compute tetrahedron volume
            vol_big = ((self.dx + B_height)
                       * (self.dx + C_depth)
                       * (self.dx + D_width)) / 6
            vol_B = B_height * B_depth * B_width / 6
            vol_C = C_height * C_depth * C_width / 6
            vol_D = D_height * D_depth * D_width / 6
            volume = vol_big - vol_B - vol_C - vol_D
            volume = np.min([volume, self.volume_cell])      
        else:
            volume = self.cas_de_merde_volume(points, pointneighbours)  
        return volume

    ########################################################################
    ########################################################################
    ########################################################################
    #         Area
    ########################################################################
    ########################################################################
    ########################################################################
    
    def compute_areas_MonteCarlo(self, points):
        """

        Args:
            points (_type_): _description_

        Returns:
            _type_: _description_
            
        Notes:
        Currently does not support having multiple boundaries in the 
        same cell, although it doesn't return an error.
        """
        
        # ========================== Intersecting points
        xvals = list(points[:,0])
        yvals = list(points[:,1])
        zvals = list(points[:,2])
        bounds = ((np.min(xvals), np.max(xvals)),
                  (np.min(yvals), np.max(yvals)),
                  (np.min(zvals), np.max(zvals)))
        
        intersect_points = []
        # intersections, starting from vertices
        for i, p1 in enumerate(points):
            # along axes and faces
            for p2 in points[i+1:]:
                # make sure it's on axis, or plane
                intersect_coord = self.intersect_3D(p1, p2)
                if intersect_coord is not None:  
                    intersect_points += [list(intersect_coord)]
            
        # extra random points with minimization
        random_points = np.array(
            [np.random.uniform(low=bounds[0][0], high=bounds[0][1], 
                               size=self.Nrand_area), 
             np.random.uniform(low=bounds[1][0], high=bounds[1][1], 
                               size=self.Nrand_area), 
             np.random.uniform(low=bounds[2][0], high=bounds[2][1], 
                               size=self.Nrand_area)]).T
        
        # find points on boundary wall with minimize function
        all_points = np.append(points, random_points, axis=0)
        for p in all_points:
            intersect_points += [list(sc.optimize.minimize(
                self.data_function_minimize, p, bounds=bounds).x)]
        
        # delete repeated points
        intersect_points.sort()
        intersect_points = np.array(
            [intersect_points 
             for intersect_points,_ in itertools.groupby(intersect_points)])
        
        # ========================== Project on 2D for triangle mesh
        
        # ID points with biggest distance 
        distAB = 0
        for i, p1 in enumerate(intersect_points):
            p1 = np.array(p1)
            p2 = np.array(intersect_points[i+1:])
            p1p2 = p2 - p1
            distp1p2 = np.sqrt(np.sum(
                p1p2**2, axis=np.min([1, np.shape(p1p2)[0]]) ) )
            j = np.argmax(distp1p2)
            if distAB < distp1p2[j]:
                ABvec = p1p2[j]
                pointA = p1
                pointB = p2[j]
                distAB = distp1p2[j]
                iA = i
                iB = i+1 + j
        
                    
        # ID 3rd point that makes the biggest traingle
        pointC = np.copy(intersect_points)
        mask = np.ones(pointC.shape[0], dtype=bool)
        mask[iA] = False
        mask[iB] = False
        pointC = pointC[mask, :] 
        distAC = np.sqrt(np.sum((pointC - pointA)**2, axis=1))
        distBC = np.sqrt(np.sum((pointC - pointB)**2, axis=1))
        areaABtri = self.triangle_area(
            [np.array([distAB]*len(distBC)), distAC, distBC])
        pointC = pointC[np.argmax(areaABtri)]
        ACvec = pointC - pointA
        
        # Define vector normal to plane defined by point A, B and C
        n = np.cross(ABvec, ACvec)
        n /= np.linalg.norm(n)
        
        # project all points onto the plane
        pointD = np.copy(intersect_points)
        pointD[iA] = np.ones(3) * 0.99 # avoid nan issues
        pointD[iB] = np.ones(3) * 0.99 # avoid nan issues
        ADvec = pointD - pointA
        ndot = np.dot(np.array([np.dot(ADvec, n)]).T, np.array([n]))
        pointDproj = pointD + ndot / np.linalg.norm(n)**2
        
        # 2D coordinate on plane with A at the origin and phi=0 along AB
        ADpvec = pointDproj - pointA
        r = np.sqrt(np.sum(ADpvec**2, axis=1))
        comb = ( np.dot(ABvec,  ADpvec.T) 
                / (np.linalg.norm(ABvec) 
                   * np.linalg.norm(ADpvec, axis=1)) )
        comb[np.where(comb<-1)] += 2e-15
        comb[np.where(1<comb)] -= 2e-15
        phi = np.arccos( comb )
        
        # fix pointB
        r[iB] = distAB
        phi[iB] = 0.0
        
        # coord on 2d plane
        xcoord = r * np.cos(phi)
        ycoord = r * np.sin(phi)
        
        # fix pointA
        xcoord[iA] = 0.0
        ycoord[iA] = 0.0
        
        # create mesh of triangles
        tri = sc.spatial.Delaunay(np.array([xcoord, ycoord]).T).simplices
        
        # ========================== Calculate area
        
        area = np.sum([self.triangle_area_from_points(
            intersect_points[i[0]], 
            intersect_points[i[1]], 
            intersect_points[i[2]]) 
                       for i in tri])
        
        # ========================== Check
        
        area_max = np.pi * self.area_cell # sphere in cell
        rel_diff = (area - area_max)/area_max
        if self.verbose and rel_diff > 1e-14:
            print('WARNING: area too big, '
                  + 'calculated {}, max {}, relative difference {}'.format(
                      area, area_max, rel_diff))
        return area
        
    def compute_areas_geometrical(self, points):
        """Compute area of the cell contained in the sphere.
        
        Parameters :
            points : (8) list
                     Each element is a (3) array_like coordinate 
                     position of the cell corner.
            points_in_sphere : (8) list
                               Each element is a boolean:
                                - True : cell corner in sphere
                                - False : cell corner not in sphere
            
        Returns :
            float, area of the cell contained in the sphere
        """
        
        corner_sign = np.sign(self.data_function(points))
        points_in_sphere = [psign==self.center_sign
                            for psign in corner_sign]
        
        nbr_points_in_sphere = np.sum(points_in_sphere)
        if nbr_points_in_sphere == 8 or nbr_points_in_sphere == 0:
            area = 0.0
        else:
            points = np.array(points)
            # Identify the neighbouring corner point of each corner.
            # For example, the corner : points[0] (x + dx/2, y + dx/2, z + dx/2)
            # Has the neighbours: 
            # points[1] (x - dx/2, y + dx/2, z + dx/2), neighbour along x
            # points[2] (x + dx/2, y - dx/2, z + dx/2), neighbour along y
            # points[6] (x + dx/2, y + dx/2, z - dx/2), neighbour along z
            neighbouring_points = np.array([[points[1], points[2], points[6]], 
                                            [points[0], points[3], points[7]], 
                                            [points[0], points[3], points[4]], 
                                            [points[1], points[2], points[5]], 
                                            [points[2], points[5], points[6]], 
                                            [points[3], points[4], points[7]], 
                                            [points[0], points[4], points[7]], 
                                            [points[1], points[5], points[6]]])
            
            sphere_mask = np.where(points_in_sphere)
            not_sphere_mask = np.where(~np.array(points_in_sphere))
            
            # Compute area.
            # The radius, corner positions and neighbouring points (masked depending
            # on the sphere) are passed to the area_*_points function 
            # according to the number of corners that are in the sphere.
            if nbr_points_in_sphere == 7:
                area = self.area_1_point(points[not_sphere_mask],
                                         neighbouring_points[not_sphere_mask])
            elif nbr_points_in_sphere == 6:
                area = self.area_2_points(points[not_sphere_mask], 
                                          neighbouring_points[not_sphere_mask])
            elif nbr_points_in_sphere == 5:
                area = self.area_3_points(points[not_sphere_mask], 
                                          neighbouring_points[not_sphere_mask])
            elif nbr_points_in_sphere == 4:
                area = self.area_4_points(points[sphere_mask], 
                                          neighbouring_points[sphere_mask])
            elif nbr_points_in_sphere == 3:
                area = self.area_3_points(points[sphere_mask],
                                          neighbouring_points[sphere_mask])
            elif nbr_points_in_sphere == 2:
                area = self.area_2_points(points[sphere_mask],
                                          neighbouring_points[sphere_mask])
            else: # nbr_points_in_sphere == 1
                area = self.area_1_point(points[sphere_mask], 
                                         neighbouring_points[sphere_mask])
        return area
    
    def cas_de_merde_area(self, points, pointneighbours):
        self.nbr_cas_de_merde += 1
        xvals = list(points[:,0]) + list(itertools.chain(*pointneighbours[:,:,0]))
        yvals = list(points[:,1]) + list(itertools.chain(*pointneighbours[:,:,1]))
        zvals = list(points[:,2]) + list(itertools.chain(*pointneighbours[:,:,2]))
        points = np.array([xvals, yvals, zvals]).T
        return self.compute_areas_MonteCarlo(points)
    
    def area_1_point(self, point, pointneighbour):  
        """Compute cell area in the sphere when 1 corner is in the sphere.
        
        Compute the area of a trirectangular tetrahedron.
        
        Parameters :
            point : (3) array_like 
                    Coordinate position of the cell corner.
            neighbourpoint : (3) array_like
                             Each element is a (3) array_like coordinate 
                             position of the cell corner neighbours.  
            
        Returns :
            float, area            
        """
        depth, ixyz = self.interpolate_get_distance(
            point[0], pointneighbour[0])
        area = self.triangle_area_from_depth(depth)
        area_max = self.triangle_area_from_depth([self.dx]*3)
        if self.verbose and (area - area_max)/area_max > 1e-14:
            print('WARNING: area 1 point too big')
        return area
        
        
    def area_2_points(self, points, pointneighbours):  
        """Compute cell area in the sphere when 2 corners are in the sphere.
        
        Compute the area of a trirectangular tetrahedron 
        that extends larger than the cell size in one direction. 
        Then remove that extension that correspondes to 
        a smaller trirectangular tetrahedron
        such that we only consider the part in the cell.
        To find the side that needs to be extended, 
        the area of each triangular base, the smallest one is extended.
        If the two areas are equal then we have a right triangular prism.
        
        Parameters :
            point : (2, 3) array_like 
                    Coordinate positions of the 2 cell corners.
            neighbourpoint : (2, 3) array_like
                             Coordinate positions of the cell corner neighbours.  
            
        Returns :
            float, area     
        """
        # Cell corners in the sphere
        point1 = points[0]
        point2 = points[1]
        
        # Neighbouring points that are outside the sphere
        pointneighbour1 = [neighbour 
                           for neighbour in pointneighbours[0] 
                           if list(neighbour)!=list(point2)]
        pointneighbour2 = [neighbour 
                           for neighbour in pointneighbours[1] 
                           if list(neighbour)!=list(point1)]
        
        # Distances between corner cell in sphere and sphere boundary.
        depth1, ixyz1 = self.interpolate_get_distance(
            point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_get_distance(
            point2, pointneighbour2)
        
        len_depth = np.sort([len(d) for d in [depth1, depth2]])
        if np.sum(len_depth == [2, 2]) == len(len_depth):
            # Area of each base
            area1 = np.prod(depth1) / 2
            area2 = np.prod(depth2) / 2
            
            if abs(area1 / area2 - 1)>1e-10: 
                # Base triangles different so trirectangular tetrahedron considered.
                # To find this: area = big_tetrahedron - small_tetrahedron
                #                      = entire_shape - extended_part
                # I'm making sure the depth/width of the triangles overlap
                # Create dict with {coord_direction: depth_val}
                dict1 = {ixyz1[i][0][0]:depth1[i] 
                        for i in range(len(ixyz1))}
                dict2 = {ixyz2[i][0][0]:depth2[i] 
                        for i in range(len(ixyz2))}
                idict = list(dict1.keys())
                if area1 > area2: # Point1 has the bigger base
                    big_width = dict1[idict[0]]
                    big_depth = dict1[idict[1]]
                    small_width = dict2[idict[0]]
                    small_depth = dict2[idict[1]]
                else: # Point2 has the bigger base
                    big_width = dict2[idict[0]]
                    big_depth = dict2[idict[1]]
                    small_width = dict1[idict[0]]
                    small_depth = dict1[idict[1]]
                    
                extended_height = (small_width * self.dx 
                                / ( big_width - small_width ))
                if extended_height > self.L or extended_height < 0:
                    area = np.average([area1, area2])
                else:
                    area_big = self.triangle_area_from_depth([self.dx + extended_height,  
                                                            big_width, big_depth])
                    area_small = self.triangle_area_from_depth([extended_height, 
                                                                small_width,
                                                                small_depth])
                    area = area_big - area_small
            else:
                # This is a right triangular prism
                area = np.average([area1, area2])
        else:
            area = self.cas_de_merde_area(points, pointneighbours)
            
        area_max = self.area_cell * np.sqrt(2)
        rel_diff = (area - area_max)/area_max
        if self.verbose and rel_diff > 1e-14:
            print('WARNING: area 2 points too big, '
                  + 'calculated {}, max {}, relative difference {}'.format(
                      area, area_max, rel_diff))
            print()
        return area
    
    def area_3_points(self, points, pointneighbours):
        """Compute cell area in the sphere when 3 corners are in the sphere.
        
        Compute the area of a trirectangular tetrahedron 
        that extends larger than the cell size in two directions. 
        Then remove those extensions that correspond to 
        smaller trirectangular tetrahedrons
        such that we only consider the part in the cell.
        
        Parameters :
            point : (3, 3) array_like 
                    Coordinate positions of the 3 cell corners.
            neighbourpoint : (3, 3) array_like
                             Coordinate positions of the cell corner neighbours.  
            
        Returns :
            float, area     
        """
        # Cell corners in the sphere
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        
        # Neighbouring points that are outside the sphere
        pointneighbour1 = [neighbour 
                           for neighbour in pointneighbours[0] 
                           if (list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point3))]
        pointneighbour2 = [neighbour 
                           for neighbour in pointneighbours[1] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point3))]
        pointneighbour3 = [neighbour 
                           for neighbour in pointneighbours[2] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point2))]
        
        # Distances between corner cell in sphere and sphere boundary.
        depth1, ixyz1 = self.interpolate_get_distance(
            point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_get_distance(
            point2, pointneighbour2)
        depth3, ixyz3 = self.interpolate_get_distance(
            point3, pointneighbour3)
        depth = [depth1, depth2, depth3]
        ixyz = [ixyz1, ixyz2, ixyz3]
       
        len_depth = np.sort([len(d) for d in depth])
        if np.sum(len_depth == [1, 2, 2]) == len(len_depth):
            # ID points
            # Point with only one depth value is A
            iA = np.where(np.array(list(map(len, depth)))==1)[0][0]
            ixyzA = ixyz[iA]
            depthA = depth[iA][0]
            # The two other points are B and C 
            # they will each have a tetrahedron extended through their side.
            iB, iC = np.delete(np.arange(3), iA)
            ixyzB = ixyz[iB]
            ixyzC = ixyz[iC]
            depthB = depth[iB]
            depthC = depth[iC]
                
            # ID depths
            A_width = depthA
            # {coord_direction: depth_val}
            dictB = {ixyzB[i][0][0]:depthB[i] for i in range(len(ixyzB))}
            dictC = {ixyzC[i][0][0]:depthC[i] for i in range(len(ixyzC))}
            # All corners have a width
            key_width = list(dictB.keys() & dictC.keys())[0]
            B_width = dictB[key_width]
            C_width = dictC[key_width]
            # Remaining coordinates correspond to depth and height
            B_depth = [dictB[key] for key in dictB if key!=key_width][0]
            C_height = [dictC[key] for key in dictC if key!=key_width][0]
            
            # Compute extended part
            B_height = B_width * self.dx / (A_width - B_width)
            C_depth = C_width * self.dx / (A_width - C_width)
                
            # Compute tetrahedron area
            area_B = self.triangle_area_from_depth([B_height, B_width, B_depth])
            area_C = self.triangle_area_from_depth([C_height, C_width, C_depth])
            area_big = self.triangle_area_from_depth([self.dx + B_height,  A_width, 
                                                    self.dx + C_depth])
            area = area_big - area_B - area_C
        else:
            area = self.cas_de_merde_area(points, pointneighbours)
        
        area_max = self.area_cell * np.sqrt(2)
        if self.verbose and (area - area_max)/area_max > 1e-14:
            print('WARNING: area 3 points too big')
        return area
    
    def area_4_points(self, points, pointneighbours):
        """Compute cell area in the sphere when 4 corners are in the sphere.
        
        Two cases :
            
        1 ) A truncated right square prism
        
        2 ) Compute the area of a trirectangular tetrahedron 
            that extends larger than the cell size in three directions. 
            Then remove those extensions that correspond to 
            smaller trirectangular tetrahedrons
            such that we only consider the part in the cell.
        
        Parameters :
            point : (4, 3) array_like 
                    Coordinate positions of the 3 cell corners.
            neighbourpoint : (4, 3) array_like
                             Coordinate positions of the cell corner neighbours.  
            
        Returns :
            float, area     
        """
        # Cell corners in the sphere
        point1 = points[0]
        point2 = points[1]
        point3 = points[2]
        point4 = points[3]
        
        # Neighbouring points that are outside the sphere
        pointneighbour1 = [neighbour 
                           for neighbour in pointneighbours[0] 
                           if (list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point3) 
                               and list(neighbour)!=list(point4))]
        pointneighbour2 = [neighbour 
                           for neighbour in pointneighbours[1] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point3) 
                               and list(neighbour)!=list(point4))]
        pointneighbour3 = [neighbour 
                           for neighbour in pointneighbours[2] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point4))]
        pointneighbour4 = [neighbour 
                           for neighbour in pointneighbours[3] 
                           if (list(neighbour)!=list(point1) 
                               and list(neighbour)!=list(point2) 
                               and list(neighbour)!=list(point3))]
        
        # Distances between corner cell in sphere and sphere boundary.
        depth1, ixyz1 = self.interpolate_get_distance(
            point1, pointneighbour1)
        depth2, ixyz2 = self.interpolate_get_distance(
            point2, pointneighbour2)
        depth3, ixyz3 = self.interpolate_get_distance(
            point3, pointneighbour3)
        depth4, ixyz4 = self.interpolate_get_distance(
            point4, pointneighbour4)
        depth = [depth1, depth2, depth3, depth4]
        
        len_depth = np.sort([len(d) for d in depth])
        if sum(len_depth) == 4:
            # A truncated right square prism
            idel = np.argmax([np.sqrt(sum((point-point1)**2))
                              for point in points[1:]])
            depth2, depth3 = [depth[i][0]
                              for i in np.delete([0,1,2,3], idel+1)[1:]]
            
            side1 = np.sqrt(self.dx**2 + abs(depth1[0]-depth2)**2)
            side2 = np.sqrt(self.dx**2 + abs(depth1[0]-depth3)**2)
            
            area = side1 * side2
        elif np.sum(len_depth == [0, 2, 2, 2]) == len(len_depth):  
            # A trirectangular tetrahedron
            # ID points
            # Point whose neighbours are all in the sphere is A
            iA = np.where(np.array(list(map(len, depth)))==0)[0][0]
            # The two other points are B, C and D
            # they will each have a tetrahedron extended through their side.
            ixyz = [ixyz1, ixyz2, ixyz3, ixyz4]
            iB, iC, iD = np.delete(np.arange(4), iA)
            ixyzB = ixyz[iB]
            ixyzC = ixyz[iC]
            ixyzD = ixyz[iD]
            depthB = depth[iB]
            depthC = depth[iC]
            depthD = depth[iD]
            
            # ID depths
            # {coord_direction: depth_val}
            dictB = {ixyzB[i][0][0]:depthB[i] for i in range(len(ixyzB))}
            dictC = {ixyzC[i][0][0]:depthC[i] for i in range(len(ixyzC))}
            dictD = {ixyzD[i][0][0]:depthD[i] for i in range(len(ixyzD))}
            # ID the depths
            key_width = list(dictB.keys() & dictC.keys())[0]
            B_width = dictB[key_width] 
            C_width = dictC[key_width] 
            key_depth = list(dictB.keys() & dictD.keys())[0]
            B_depth = dictB[key_depth] 
            D_depth = dictD[key_depth] 
            key_height = list(dictC.keys() & dictD.keys())[0]
            C_height = dictC[key_height] 
            D_height = dictD[key_height] 
            
            # Compute extended part
            B_height = B_depth * (self.dx - C_height) / (self.dx - B_depth)
            C_depth = C_width * (self.dx - D_depth) / (self.dx - C_width)
            D_width = D_height * (self.dx - B_width) / (self.dx - D_height)
            
            # Compute tetrahedron area
            area_big = self.triangle_area_from_depth([self.dx + B_height,
                                                      self.dx + C_depth,
                                                      self.dx + D_width])
            area_B = self.triangle_area_from_depth([B_height, B_depth, B_width])
            area_C = self.triangle_area_from_depth([C_height, C_depth, C_width])
            area_D = self.triangle_area_from_depth([D_height, D_depth, D_width])
            area = area_big - area_B - area_C - area_D   
        else:
            area = self.cas_de_merde_area(points, pointneighbours)  
            
        area_max = self.area_cell * np.sqrt(2)
        if self.verbose and (area - area_max)/area_max > 1e-14:
            print('WARNING: area 4 points too big')
        return area
    

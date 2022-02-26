#######################################################################
#             Computational and differential Geometry                 #
#                         Second Homework                             #
#                 Andres Felipe Florian Quitian                       #
#######################################################################

import numpy as np
import pandas as pd
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial.distance import euclidean, cdist
import matplotlib.path as mplPath
import matplotlib.pyplot as plt


class FlightOperations:

    """ Class FlightOperations that represents a set of airports
        and its corresponding Voronoi Diagram

        Attributes: - Coordinates: Longitude and Latitude for each airport
                    - Altitude: Airport altitude
                    - Cities: Cities where the airports are located
                    - Names: Airport names
                    - x_coord: Points in the x-axis for the borders (map)
                    - y_coord: Points in the y-axis fot the borders (map)
                    - vor: Voronoi Diagram of the coordinates
    """

    def __init__(self, borders, airports):

        """ Constructor: Class atributes initialization

            Arguments: - borders: str. Name of the file storing
                                       the borders data (map)
                       - airports: str. Name of the file storing
                                        the airports data
        """

        airp_file = pd.read_csv(airports, sep=r'\s+',
                                names=['Latitude', 'Longitude', 'Altitude',
                                       'City', 'Departament', 'Name'])

        self.coordinates = np.array([[airp_file['Longitude'][i],
                                      airp_file['Latitude'][i]]
                                    for i in range(len(airp_file))])

        self.altitude = airp_file['Altitude']
        self.cities = airp_file['City']
        self.names = airp_file['Name']
        self.vor = Voronoi(self.coordinates)

        bor_file = pd.read_csv(borders, sep=r'\s+',
                               names=['Latitude', 'Longitude'])

        self.x_coord = np.array([bor_file['Longitude'][i]
                                for i in range(len(bor_file))])
        self.y_coord = np.array([bor_file['Latitude'][i]
                                for i in range(len(bor_file))])

    def _area(self, region):

        """ Private method that calculates the area of a given region
            of the Voronoi Diagram

            Input: - region: list that represents the region
                             for the Voronoi Subdivision
            Output: float. Area of the region; or Infinite Area if -1 in region
        """

        if -1 in region:
            return "Infinite Area"
        vertices = [self.vor.vertices[i].tolist() for i in region]
        area = 0
        for i in range(len(vertices)):
            if i == len(vertices)-1:
                x = vertices[0][1]*vertices[i][0]
                x -= vertices[0][0]*vertices[i][1]
            else:
                x = vertices[i][0]*vertices[i+1][1]
                x -= vertices[i][1]*vertices[i+1][0]
            area += x
        return abs(area)/2

    def _inside_outside(self, poly, vertex):

        """ Private method that verifies if a vertex is in a polygon

            Input: - poly: list with the vertices of the polygon
            Output True or False. vertex in or out the polygon
        """

        poly_path = mplPath.Path(np.array(poly))
        return poly_path.contains_point(vertex)

    def _sort(self, x_coords, y_coords):

        """ Private method that sort in counterclockwise order
            the given coordinates

            Input: - x_coords: array with the coordinates in x
                   - y_coords: array with the coordinates in y
            Output: x_new, y_new coordinates sorted
        """

        x_0 = np.mean(x_coords)
        y_0 = np.mean(y_coords)

        r = np.sqrt((x_coords-x_0)**2 + (y_coords-y_0)**2)

        points_angles = np.where((y_coords-y_0) > 0,
                                 np.arccos((x_coords-x_0)/r),
                                 2*np.pi-np.arccos((x_coords-x_0)/r))

        sort_ = np.argsort(points_angles)

        new_x = x_coords[sort_]
        new_y = y_coords[sort_]

        return new_x, new_y

    def plot_vor_airports(self):

        """ Method that plots the Voronoi subdivision and the borders """

        fig = voronoi_plot_2d(self.vor, show_vertices=False, line_alpha=0.8)

        plt.plot(self.x_coord, self.y_coord, 'k')

        plt.title('Voronoi Map')
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.show()

    def max_min_areas(self, plot=False):

        """ Method that returns the airports where is reported the max
            and min area in the subdivision for region with finite areas

            Input(optional): - plot: boolean for plotting
            Output: (str, str). The airports with min and max coverage
                                area respectively
        """

        regions = self.vor.regions
        min_index = max_index = 0
        point_reg = self.vor.point_region
        max_area = min_area = self._area(regions[point_reg[0]])
        for i in range(1, len(point_reg)):
            region = regions[point_reg[i]]
            if region != []:
                new_area = self._area(region)
            else:
                new_area = "Infinite Area"
            if new_area != "Infinite Area":
                if new_area > max_area:
                    max_area = new_area
                    max_index = i
                if new_area < min_area:
                    min_area = new_area
                    min_index = i
        min_airport = self.names[min_index]
        max_airport = self.names[max_index]

        if plot:
            fig = voronoi_plot_2d(self.vor, show_vertices=False,
                                  line_alpha=0.8)
            plt.plot(self.x_coord, self.y_coord, 'k')

            vertex_1 = self.coordinates[min_index]
            vertex_2 = self.coordinates[max_index]
            min_region = regions[point_reg[min_index]]
            max_region = regions[point_reg[max_index]]
            min_ver = [self.vor.vertices[i] for i in min_region]
            max_ver = [self.vor.vertices[i] for i in max_region]

            plt.fill(*zip(*min_ver), color='#f5cba7')
            plt.fill(*zip(*max_ver), color='#f5cba7')

            plt.plot(vertex_1[0], vertex_1[1], marker='*',
                     c='green', label='Min Area: ' + min_airport)
            plt.plot(vertex_2[0], vertex_2[1], marker='*',
                     c='yellow', label='Max Area: ' + max_airport)

            fig.legend(loc='lower right')
            plt.title('Maximum and minimun areas')
            plt.xlabel('Longitude (degrees)')
            plt.ylabel('Latitude (degrees)')
            plt.show()

        return min_airport, max_airport

    def build_airport(self, plot=False):

        """ Method that returns the center and radius of the
            largest circle centered within the convex hull and
            enclosing none of the already existing airports

            Input(optional): - plot: boolean for plotting
            Output: (list, float). list with the coordinates
                                   of center and its radius
        """

        vertices = self.vor.vertices
        point_reg = self.vor.point_region
        regions = self.vor.regions
        index = -1
        max_dis = -1
        max_vertex = []
        borders = [(self.x_coord[i], self.y_coord[i])
                   for i in range(len(self.x_coord))]

        edges = self.vor.ridge_vertices
        inf_edges = []
        for j in range(len(edges)):
            if -1 in edges[j]:
                index = self.vor.ridge_points[j]
                inf_edges.append([self.vor.points[index[0]],
                                  self.vor.points[index[1]]])

        vertices_hull = []
        for ver in inf_edges:
            coord_1 = ver[0].tolist()
            coord_2 = ver[1].tolist()
            if coord_1 not in vertices_hull:
                vertices_hull.append(coord_1)
            if coord_2 not in vertices_hull:
                vertices_hull.append(coord_2)

        x = np.array([i[0] for i in vertices_hull])
        y = np.array([i[1] for i in vertices_hull])
        x, y = self._sort(x, y)
        hull = [(x[i], y[i]) for i in range(len(x))]

        for i in range(len(vertices)):
            vertex = vertices[i].tolist()
            check_1 = self._inside_outside(hull, vertex)
            check_2 = self._inside_outside(borders, vertex)
            if check_1 and check_2:
                for j, r in enumerate(regions):
                    if i in r:
                        index = self.vor.point_region.tolist().index(j)
                        pos = self.vor.points[index]
                        dis = euclidean(vertex, pos)
                        if dis > max_dis:
                            max_dis = dis
                            max_vertex = vertex
        center = max_vertex[0], max_vertex[1]

        if plot:

            fig, ax = plt.subplots()

            ax.plot(x, y, c='red')
            ax.plot([x[len(x)-1], x[0]], [y[len(x)-1], y[0]], c='red')

            voronoi_plot_2d(self.vor, ax=ax, show_vertices=False,
                            line_alpha=0.8)

            ax.plot(self.x_coord, self.y_coord, 'k')
            ax.plot(max_vertex[0], max_vertex[1], marker='*',
                    c='green', label='Center of largest circle')

            circle = plt.Circle(center, max_dis,
                                fill=True, color='purple')
            ax.add_patch(circle)

            fig.legend(loc='lower right')
            plt.title('Build New Airport')
            plt.xlabel('Longitude (degrees)')
            plt.ylabel('Latitude (degrees)')
            plt.show()

        return center, max_dis

    def most_less_crowded(self, plot=False):

        """ Method that returns the airports with most and least neighbors
            respectively

            Input(optional): - plot: boolean for plotting
            Output: (str, str). The airports with least and
                                most neighbors respectively
        """

        total_airports = len(self.coordinates)
        edges = self.vor.ridge_points
        neighbors = total_airports*[0]
        for i in edges:
            neighbors[i[0]] += 1
            neighbors[i[1]] += 1
        min_index = neighbors.index((min(neighbors)))
        max_index = neighbors.index((max(neighbors)))
        most_airport = self.names[max_index]
        least_airport = self.names[min_index]

        if plot:
            fig = voronoi_plot_2d(self.vor, show_vertices=False,
                                  line_alpha=0.8)

            plt.plot(self.x_coord, self.y_coord, 'k')

            vertex_1 = self.coordinates[min_index]
            vertex_2 = self.coordinates[max_index]
            point_reg = self.vor.point_region
            min_region = self.vor.regions[point_reg[min_index]]
            max_region = self.vor.regions[point_reg[max_index]]
            min_ver = [self.vor.vertices[i] for i in min_region]
            max_ver = [self.vor.vertices[i] for i in max_region]

            plt.fill(*zip(*min_ver), color='#f5cba7')
            plt.fill(*zip(*max_ver), color='#f5cba7')

            plt.plot(vertex_1[0], vertex_1[1], marker='*', c='green',
                     label='Most Crowded Airport: ' + most_airport)
            plt.plot(vertex_2[0], vertex_2[1], marker='*', c='yellow',
                     label='Least Crowded Airport: ' + least_airport)

            fig.legend(loc='lower right')
            plt.title('Most and least Crowded Airports')
            plt.xlabel('Longitude (degrees)')
            plt.ylabel('Latitude (degrees)')
            plt.show()

        return least_airport, most_airport

    def merge_airport(self, plot=False):

        """ Method that returns the airports that may be
            deleted and the position of the new airport

            Input(optional): - plot: boolean for plotting
            Output: (str, str, array). The airports that may be deleted,
                                       the coordinates of the mean between them
        """

        points = self.coordinates.tolist()
        min_dis = euclidean(points[0], points[1])
        coords = points[0], points[1]
        for i in range(len(points)):
            x = points[i]
            for j in range(i+1, len(points)):
                y = points[j]
                new_dis = euclidean(x, y)
                if new_dis < min_dis:
                    min_dis = new_dis
                    coords = x, y
        coor_1 = points.index(coords[0])
        coor_2 = points.index(coords[1])
        airport_1 = self.names[coor_1]
        airport_2 = self.names[coor_2]
        pos = np.mean(np.array([coords[0], coords[1]]), axis=0)

        if plot:
            fig = voronoi_plot_2d(self.vor, show_vertices=False,
                                  line_alpha=0.8)

            plt.plot(self.x_coord, self.y_coord, 'k')

            vertex_1 = self.coordinates[coor_1]
            vertex_2 = self.coordinates[coor_2]

            plt.plot(vertex_1[0], vertex_1[1], marker='*', c='green',
                     label=airport_1 + ' out')
            plt.plot(vertex_2[0], vertex_2[1], marker='*', c='yellow',
                     label=airport_2 + ' out')
            plt.plot(pos[0], pos[1], marker='*', c='red',
                     label='New Airport')

            fig.legend(loc='lower right')
            plt.title('Merge Airport')
            plt.xlabel('Longitude (degrees)')
            plt.ylabel('Latitude (degrees)')
            plt.show()

        return airport_1, airport_2, pos

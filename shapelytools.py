from shapely.geometry import (box, LineString, MultiLineString, MultiPoint, 
    Point, Polygon)
import shapely.ops
import shapefile
import itertools
import pdb
import progress as pbar




def write_shp(filename, geometry, records=[], fields=[]):
        """Write a single shapely MultiLineString or Polygon to a shapefile.
        Argument geometry may also be a list of LineString objects. In that case,
        each entry may have associated data in the optional records list. If
        records is given, fields is the list of fieldnames for the value list in
        each record.
        Usage:
                write_shp(filename, geometry, records=[], fields=[])
        Arguments:
                filename    filename of shapefile
                geometry    a shapely geometry (MultiLineString, Polygon, ...)
                records     optional list of list of values, one per geometry
                fields      optional (implied by records) list of fieldnames
        """
    
        # SINGLE MULTILINESTRING
        if isinstance(geometry, MultiLineString):
            
            sw = shapefile.Writer(shapefile.POLYLINE)
            
            # fields
            sw.field("length")
            sw.field("start-x")
            sw.field("start-y")
            sw.field("end-x")
            sw.field("end-y")
            sw.field("npoints")
            
            # geometry and record
            for line in geometry:
                    sw.line([list(line.coords)])
                    sw.record(line.length,
                                        line.coords[0][0],
                                        line.coords[0][1],
                                        line.coords[-1][0],
                                        line.coords[-1][1],
                                        len(line.coords))
            sw.save(filename)
            
        # SINGLE POLYGON
        elif isinstance(geometry, Polygon):
            # data
            parts = [list(geometry.exterior.coords)]
            parts.extend(list(interior.coords) for interior in geometry.interiors)
            
            sw = shapefile.Writer(shapefile.POLYGON)
            sw.field("area")
            sw.poly(parts)
            sw.record(geometry.area)
            sw.save(filename)
            
        # LISTS
        elif isinstance(geometry, list):
            # check for fields and records
            if (fields and not records) or (not fields and records):
                raise ValueError('Arguments records and fields must both be provided.')
#                
            if records and (len(fields) != len(records[0])):
                raise ValueError('Length of fields and records do not match.')
#                
            # derive field types based on record values
            field_type = {}
            precision = {}
            for i, field in enumerate(fields):
#                print(i,field)
                if all(type(record[i]) in [int, int] for record in records):
                    field_type[field] = 'N' # integers
                    precision[field] = 0
                elif all(type(record[i]) in [int, int, float] for record in records):
                    field_type[field] = 'N' # numeric
                    precision[field] = 0#5
                else:
                    field_type[field] = 'C' # string (characters)
                    precision[field] = 0
#                print(i,field,field_type[field])
            # LIST OF LINESTRINGS
            if isinstance(geometry[0], LineString):
                
                print(shapefile.POLYLINE)
                sw = shapefile.Writer(filename)
                # fields
                for field in fields:
                    sw.field(field, field_type[field], decimal=precision[field])
                if not fields:
                    # add dummy length field and prepare records for it
                    records = [[line.length] for line in geometry]
                    sw.field('length', 'N', decimal=5)
                # geometry and records
                for line, record in zip(geometry, records):
                    sw.line([list(line.coords)])
                    sw.record(*record)
                sw.close()
            # LIST OF POLYGONS
            elif isinstance(geometry[0], Polygon):
                sw = shapefile.Writer(shapefile.POLYGON)
                # fields
                for field in fields:
                    sw.field(field, field_type[field], decimal=precision[field])
                # geometry and records
                for polygon, record in itertools.izip(geometry, records):
                    sw.poly([list(polygon.exterior.coords)])
                    sw.record(*record)
                sw.save(filename)
            # LIST OF POINTS
            elif isinstance(geometry[0], Point):
                sw = shapefile.Writer(shapefile.POINT)
                
                # fields
                for field in fields:
                    sw.field(field, field_type[field], decimal=precision[field])
                    
                if not fields:
                    # add dummy length field and prepare records for it
                    records = [[point.x, point.y] for point in geometry]
                    sw.field('x', 'N', decimal=5)
                    sw.field('y', 'N', decimal=5)
                    
                # shapes and records
                for point, record in itertools.izip(geometry, records):
                    sw.point(point.x, point.y)
                    sw.record(*record)
                    
                # save
                sw.save(filename)
            else:
                raise NotImplementedError
                
        else:
            raise NotImplementedError
            
            
def read_shp(filename):
    """Read contents of a shapefile to a shapely geometry object.
    Usage:
        geometries = read_shp(filename)
        (geometries, records, fields) = read_shp(filename)
    Arguments:
        filename    shapefile name
    Returns:
        geometries  list of shapely geometries (Polygon, LineString, ...)
        records     list of records (list of values), one per geometry
        fields      list of fieldnames of a record
    """
    sr = shapefile.Reader(filename)
    
    if sr.shapeType == shapefile.POLYGON:
        shapes = sr.shapes()
        geometries = [Polygon(shape.points) for shape in shapes]
        
        fields = sr.fields[:]
        if fields[0][0] == 'DeletionFlag':
            fields.pop(0)
        fields = [field[0] for field in fields] # extract field name only
        
        records = []
        for record in sr.records():
            for i, value in enumerate(record):
                try:
                    record[i] = float(value) # convert record values to numeric...
                except ValueError:
                    pass # ... if possible
                    
            records.append(record)
            
        return (geometries, records, fields)
    
    elif sr.shapeType == shapefile.POLYLINE:
        shapes = sr.shapes()
        geometries = [LineString(shape.points) for shape in shapes]
        
        fields = sr.fields[:] # [:] = duplicate field list
        if fields[0][0] == 'DeletionFlag':
            fields.pop(0)
        fields = [field[0] for field in fields] # extract field name only
        
        records = []
        for record in sr.records():
            for i, value in enumerate(record):
                try:
                    record[i] = float(value) # convert record values to numeric...
                except ValueError:
                    pass # ... if possible
                    
            records.append(record)
            
        return (geometries, records, fields)
    
    
    elif sr.shapeType == shapefile.MULTIPOINT:
        raise NotImplementedError
        
    else:
        raise NotImplementedError


def endpoints_from_lines(lines):
    """Return list of terminal points from list of LineStrings."""
    
    all_points = []
    for line in lines:
        for i in [0, -1]: # start and end point
            all_points.append(line.coords[i])
    
    unique_points = set(all_points)
    
    return [Point(p) for p in unique_points]
    
def vertices_from_lines(lines):
    """Return list of unique vertices from list of LineStrings."""
    count = len(lines)
    print("Getting vertices 1/3")
    pb = pbar.ProgressBar(count)
    vertices = []
#    print("getting vertices from line")
    for line in lines:
        pb +=1
        vertices.extend(list(line.coords))
    del pb
    return [Point(p) for p in set(vertices)]


def prune_short_lines(lines, min_length):
    """Remove lines from a LineString DataFrame shorter than min_length.
    
    Deletes all lines from a list of LineStrings or a MultiLineString
    that have a total length of less than min_length. Vertices of touching 
    lines are contracted towards the centroid of the removed line.
    
    Args:
        lines: list of LineStrings or a MultiLineString
        min_length: minimum length of a single LineString to be preserved
        
    Returns:
        the pruned pandas DataFrame
    """   
    pruned_lines = [line for line in lines] # converts MultiLineString to list
    to_prune = []
    
    for i, line in enumerate(pruned_lines):
        if line.length < min_length:
            to_prune.append(i)
            for n in neighbors(pruned_lines, line):
                contact_point = line.intersection(pruned_lines[n])
                pruned_lines[n] = bend_towards(pruned_lines[n], 
                                               where=contact_point,
                                               to=line.centroid)
                
    return [line for i, line in enumerate(pruned_lines) if i not in to_prune] 


def neighbors(lines, of):
    """Find the indices in a list of LineStrings that touch a given LineString.
    
    Args:
        lines: list of LineStrings in which to search for neighbors
        of: the LineString which must be touched
        
    Returns:
        list of indices, so that all lines[indices] touch the LineString of
    """
    return [k for k, line in enumerate(lines) if line.touches(of)]
    

def bend_towards(line, where, to):
    """Move the point where along a line to the point at location to.
    
    Args:
        line: a LineString
        where: a point ON the line (not necessarily a vertex)
        to: a point NOT on the line where the nearest vertex will be moved to
    
    Returns:
        the modified (bent) line 
    """
    
    if not line.contains(where) and not line.touches(where):
        raise ValueError('line does not contain the point where.')
        
    coords = line.coords[:]
    # easy case: where is (within numeric precision) a vertex of line
    for k, vertex in enumerate(coords):
        if where.almost_equals(Point(vertex)):
            # move coordinates of the vertex to destination
            coords[k] = to.coords[0]
            return LineString(coords)
    
    # hard case: where lies between vertices of line, so
    # find nearest vertex and move that one to point to
    _, min_k = min((where.distance(Point(vertex)), k) 
                           for k, vertex in enumerate(coords))
    coords[min_k] = to.coords[0]
    return LineString(coords)


def snappy_endings(lines, max_distance):
    """Snap endpoints of lines together if they are at most max_length apart.
    
    Args:
        lines: a list of LineStrings or a MultiLineString
        max_distance: maximum distance two endpoints may be joined together 
    """
    
    # initialize snapped lines with list of original lines
    # snapping points is a MultiPoint object of all vertices
    snapped_lines = [line for line in lines]
    snapping_points = vertices_from_lines(snapped_lines)
    
    
    # isolated endpoints are going to snap to the closest vertex
    isolated_endpoints = find_isolated_endpoints(snapped_lines)  
    
    
    # only move isolated endpoints, one by one
#    accum=1
    
    
    count = len(isolated_endpoints)
    print("Performing line snapping 3/3")
    pb = pbar.ProgressBar(count)
    
    for endpoint in isolated_endpoints:
        pb += 1
        # find all vertices within a radius of max_distance as possible
        target = nearest_neighbor_within(snapping_points, endpoint, 
                                         max_distance)
        
        # do nothing if no target point to snap to is found
        if not target:
            continue       
        
        # find the LineString to modify within snapped_lines and update it        
        for i, snapped_line in enumerate(snapped_lines):
            if endpoint.touches(snapped_line):
                snapped_lines[i] = bend_towards(snapped_line, where=endpoint, 
                                                to=target)
                break
        
        # also update the corresponding snapping_points
        for i, snapping_point in enumerate(snapping_points):
            if endpoint.equals(snapping_point):
                snapping_points[i] = target
                break

    # post-processing: remove any resulting lines of length 0
    snapped_lines = [s for s in snapped_lines if s.length > 0]
    del pb
    return snapped_lines
    
    
def nearest_neighbor_within(others, point, max_distance):
    """Find nearest point among others up to a maximum distance.
    
    Args:
        others: a list of Points or a MultiPoint
        point: a Point
        max_distance: maximum distance to search for the nearest neighbor
        
    Returns:
        A shapely Point if one is within max_distance, None otherwise
    """
    search_region = point.buffer(max_distance)
    interesting_points = search_region.intersection(MultiPoint(others))
    
    if not interesting_points:
        closest_point = None
    elif isinstance(interesting_points, Point):
        closest_point = interesting_points
    else:            
        distances = [point.distance(ip) for ip in interesting_points
                     if point.distance(ip) > 0]
        closest_point = interesting_points[distances.index(min(distances))]
    
    return closest_point


def find_isolated_endpoints(lines):
    """Find endpoints of lines that don't touch another line.
    
    Args:
        lines: a list of LineStrings or a MultiLineString
        
    Returns:
        A list of line end Points that don't touch any other line of lines
    """
        
    isolated_endpoints = []
    count = len(lines)
    print("Finding isolated end points 2/3")
    pb = pbar.ProgressBar(count)
    for i, line in enumerate(lines):
        pb += 1
        other_lines = lines[:i] + lines[i+1:]
        for q in [0,-1]:
            endpoint = Point(line.coords[q])
            if any(endpoint.touches(another_line) 
                   for another_line in other_lines):
                continue
            else:
                isolated_endpoints.append(endpoint)
    del pb
    return isolated_endpoints
    
def closest_object(geometries, point):
    """Find the nearest geometry among a list, measured from fixed point.
    
    Args:
        geometries: a list of shapely geometry objects
        point: a shapely Point
       
    Returns:
        Tuple (geom, min_dist, min_index) of the geometry with minimum distance 
        to point, its distance min_dist and the list index of geom, so that
        geom = geometries[min_index].
    """    
    min_dist, min_index = min((point.distance(geom), k) 
                              for (k, geom) in enumerate(geometries))
    
    return geometries[min_index], min_dist, min_index
    
    
def project_point_to_line(point, line_start, line_end):
    """Find nearest point on a straight line, measured from given point.
    
    Args:
        point: a shapely Point object
        line_start: the line starting point as a shapely Point
        line_end: the line end point as a shapely Point
    
    Returns:
        a shapely Point that lies on the straight line closest to point
    
    Source: http://gis.stackexchange.com/a/438/19627
    """
    line_magnitude = line_start.distance(line_end)
    
    u = ((point.x - line_start.x) * (line_end.x - line_start.x) +
         (point.y - line_start.y) * (line_end.y - line_start.y)) \
         / (line_magnitude ** 2)

    # closest point does not fall within the line segment, 
    # take the shorter distance to an endpoint
    if u < 0.00001 or u > 1:
        ix = point.distance(line_start)
        iy = point.distance(line_end)
        if ix > iy:
            return line_end
        else:
            return line_start
    else:
        ix = line_start.x + u * (line_end.x - line_start.x)
        iy = line_start.y + u * (line_end.y - line_start.y)
        return Point([ix, iy])
        
def pairs(lst):
    """Iterate over a list in overlapping pairs.
    
    Args:
        lst: an iterable/list
        
    Returns:
        Yields a pair of consecutive elements (lst[k], lst[k+1]) of lst. Last 
        call yields (lst[-2], lst[-1]).
        
    Example:
        lst = [4, 7, 11, 2]
        pairs(lst) yields (4, 7), (7, 11), (11, 2)
       
    Source:
        http://stackoverflow.com/questions/1257413/1257446#1257446
    """
    i = iter(lst)
    prev = next(i)
    for item in i:
        yield prev, item
        prev = item


def project_point_to_object(point, geometry):
    """Find nearest point in geometry, measured from given point.
    
    Args:
        point: a shapely Point
        geometry: a shapely geometry object (LineString, Polygon)
        
    Returns:
        a shapely Point that lies on geometry closest to point
    """
    nearest_point = None
    min_dist = float("inf")
    
    if isinstance(geometry, Polygon):
        for seg_start, seg_end in pairs(list(geometry.exterior.coords)):
            line_start = Point(seg_start)
            line_end = Point(seg_end)
        
            intersection_point = project_point_to_line(point, line_start, line_end)
            cur_dist =  point.distance(intersection_point)
        
            if cur_dist < min_dist:
                min_dist = cur_dist
                nearest_point = intersection_point
    
    elif isinstance(geometry, LineString):
        for seg_start, seg_end in pairs(list(geometry.coords)):
            line_start = Point(seg_start)
            line_end = Point(seg_end)
        
            intersection_point = project_point_to_line(point, line_start, line_end)
            cur_dist =  point.distance(intersection_point)
        
            if cur_dist < min_dist:
                min_dist = cur_dist
                nearest_point = intersection_point
    else:
        raise NotImplementedError("project_point_to_object not implemented for"+
                                  " geometry type '" + geometry.type + "'.")
    return nearest_point
    

def one_linestring_per_intersection(lines):
    """ Move line endpoints to intersections of line segments.
    
    Given a list of touching or possibly intersecting LineStrings, return a
    list LineStrings that have their endpoints at all crossings and
    intersecting points and ONLY there.
    
    Args:
        a list of LineStrings or a MultiLineString
        
    Returns:
        a list of LineStrings
    """
    lines_merged = shapely.ops.linemerge(lines)

    # intersecting multiline with its bounding box somehow triggers a first
    bounding_box = box(*lines_merged.bounds)

    # perform linemerge (one linestring between each crossing only)
    # if this fails, write function to perform this on a bbox-grid and then
    # merge the result
    lines_merged = lines_merged.intersection(bounding_box)
    lines_merged = shapely.ops.linemerge(lines_merged)
    return lines_merged


def linemerge(linestrings_or_multilinestrings):
    """ Merge list of LineStrings and/or MultiLineStrings.
    
    Given a list of LineStrings and possibly MultiLineStrings, merge all of
    them to a single MultiLineString.
    
    Args:
        list of LineStrings and/or MultiLineStrings
    
    Returns:
        a merged LineString or MultiLineString
    """
    lines = []
    for line in linestrings_or_multilinestrings:
        if isinstance(line, MultiLineString):
            # line is a multilinestring, so append its components
            lines.extend(line)
        else:
            # line is a line, so simply append it
            lines.append(line)
        
    return shapely.ops.linemerge(lines)
    
    
import numpy as np



#########################    UTILITY FUNCTIONS FOR PLOTTING   ################################
def position_to_coords(position_name, x_range, y_range):
    """
    Converts a position name (e.g., 'upper right') to real coordinates based on the x and y range of the plot.

    Parameters:
        position_name (str): Name of the position (e.g., 'lower left', 'upper right', 'center').
        x_range (tuple): The x range of the plot as (xmin, xmax).
        y_range (tuple): The y range of the plot as (ymin, ymax).

    Returns:
        tuple: Real coordinates (x, y).
    """
    # Dictionary mapping position names to normalized coordinates
    position_to_normalized = {
        'lower left': (0.0, 0.0),
        'lower right': (1.0, 0.0),
        'upper left': (0.0, 1.0),
        'upper right': (1.0, 1.0),
        'center': (0.5, 0.5),
        'center left': (0.0, 0.5),
        'center right': (1.0, 0.5),
        'lower center': (0.5, 0.0),
        'upper center': (0.5, 1.0),
    }

    # Get the normalized coordinates for the given position name
    if position_name in position_to_normalized:
        nx, ny = position_to_normalized[position_name]
    else:
        raise ValueError(f"Unknown position name: {position_name}. Valid options are: {list(position_to_normalized.keys())}")

    # Convert normalized coordinates to real coordinates
    xmin, xmax = x_range
    ymin, ymax = y_range
    x = xmin + nx * (xmax - xmin)
    y = ymin + ny * (ymax - ymin)

    return x, y


def closest_point_on_segment(p, a, b):

    ap = np.array(p) - np.array(a)  # Vector from a to p
    ab = np.array(b) - np.array(a)  # Vector from a to b

    # Project ap onto ab
    t = np.dot(ap, ab) / np.dot(ab, ab)

    # Clamp t to the range [0, 1] to ensure the point lies on the segment
    t = max(0, min(1, t))

    # Calculate the closest point
    closest = np.array(a) + t * ab
    return tuple(closest)


def closest_point_on_rectangle(p, rectangle):

    # Extract rectangle coordinates
    xmin, ymin = rectangle.get_xy()  # Bottom-left corner
    width = rectangle.get_width()
    height = rectangle.get_height()
    xmax = xmin + width
    ymax = ymin + height

    # Define the four sides of the rectangle as line segments
    sides = [
        ((xmin, ymin), (xmax, ymin)),  # Bottom side
        ((xmax, ymin), (xmax, ymax)),  # Right side
        ((xmax, ymax), (xmin, ymax)),  # Top side
        ((xmin, ymax), (xmin, ymin)),  # Left side
    ]
    

    # Find the closest point on each side
    closest_points = [closest_point_on_segment(p, a, b) for a, b in sides]

    # Find the closest point overall
    distances = [np.hypot(p[0] - cp[0], p[1] - cp[1]) for cp in closest_points]
    closest_index = np.argmin(distances)

    return closest_points[closest_index]

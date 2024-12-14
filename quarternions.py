import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def quaternion_rotate(vector, axis, angle):
    """Rotate a vector using a quaternion."""
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)  # Ensure axis is a unit vector
    half_angle = angle / 2
    w = np.cos(half_angle)
    x, y, z = axis * np.sin(half_angle)
    q = np.array([w, x, y, z])
    q_conj = np.array([w, -x, -y, -z])  # Conjugate for inverse

    # Convert vector to quaternion (0 + xi + yj + zk)
    v_quat = np.array([0, *vector])

    # Quaternion rotation: q * v * q^-1
    rotated = quaternion_multiply(quaternion_multiply(q, v_quat), q_conj)
    return rotated[1:]  # Return the vector part

def quaternion_multiply(q1, q2):
    """Multiply two quaternions."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def generate_branches(initial_vector, num_levels, azimuth_angle, bifurcation_angle):
    """Generate branching vectors with rotations."""
    lines = []  # To store line segments (start, end pairs)
    start_points = [np.array([0, 0, 0])]  # Start at origin
    directions = [initial_vector]

    for level in range(num_levels):
        new_start_points = []
        new_directions = []
        for start, direction in zip(start_points, directions):
            # Normalize direction to maintain consistent vector magnitudes
            direction = direction / np.linalg.norm(direction)

            # Rotate for the first branch
            branch1_dir = quaternion_rotate(
                direction, axis=[0, 0, 1], angle=azimuth_angle
            )
            branch1_dir = quaternion_rotate(
                branch1_dir, axis=np.cross(direction, [0, 0, 1]), angle=bifurcation_angle
            )
            branch1_end = start + branch1_dir
            lines.append([start, branch1_end])

            # Rotate for the second branch
            branch2_dir = quaternion_rotate(
                direction, axis=[0, 0, 1], angle=-azimuth_angle
            )
            branch2_dir = quaternion_rotate(
                branch2_dir, axis=np.cross(direction, [0, 0, 1]), angle=bifurcation_angle
            )
            branch2_end = start + branch2_dir
            lines.append([start, branch2_end])

            # Update new starts and directions
            new_start_points.extend([branch1_end, branch2_end])
            new_directions.extend([branch1_dir, branch2_dir])

        # Update start points and directions for the next level
        start_points = new_start_points
        directions = new_directions

    return np.array(lines)

def plot_branches(lines):
    """Plot the branches as 3D line segments."""
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each line segment
    for line in lines:
        start, end = line
        x_coords, y_coords, z_coords = zip(start, end)
        ax.plot(x_coords, y_coords, z_coords, color='green')

    # Adjust the axes for equal scaling
    all_points = np.concatenate(lines)
    x_limits = [np.min(all_points[:, 0]), np.max(all_points[:, 0])]
    y_limits = [np.min(all_points[:, 1]), np.max(all_points[:, 1])]
    z_limits = [np.min(all_points[:, 2]), np.max(all_points[:, 2])]

    ax.set_xlim(x_limits)
    ax.set_ylim(y_limits)
    ax.set_zlim(z_limits)

    # Force equal aspect ratio
    max_range = np.array([np.ptp(x_limits), np.ptp(y_limits), np.ptp(z_limits)]).max() / 2.0
    mid_x = np.mean(x_limits)
    mid_y = np.mean(y_limits)
    mid_z = np.mean(z_limits)
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Set axis labels and show plot
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title("3D Plant Branching")
    plt.show()

# Parameters
initial_vector = np.array([0, 0, 1])  # Initial upward vector
num_levels = 3                        # Number of branching levels
azimuth_angle = np.pi / 6             # Rotation around vertical (30 degrees)
bifurcation_angle = np.pi / 4         # Tilt from vertical (45 degrees)

# Generate and plot branches
lines = generate_branches(initial_vector, num_levels, azimuth_angle, bifurcation_angle)
plot_branches(lines)
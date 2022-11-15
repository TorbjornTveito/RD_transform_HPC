import matplotlib.pyplot as plt



def load_points():
    with open('./data/point_data/a0/3_points.dat', 'r') as fp:
        data = [line for line in fp]
    p1 = []
    p2 = []
    for line in data:
        if line[0]== '#':
            continue

        points = line.split(' ')
        p1x, p1y = points[0].split(',')
        p2x, p2y = points[1].split(',')
        p1.append((float(p1x), float(p1y)))
        p2.append((float(p2x), float(p2y)))

    return p1, p2

def plot_points():
    points1, points2 = load_points()

    points1_lat = [point[0] for point in points1]
    points1_lon = [point[1] for point in points1]
    points2_lat = [point[0] for point in points2]
    points2_lon = [point[1] for point in points2]

    plt.plot(points1_lon, points1_lat, '.',  c = "red")
    plt.plot(points2_lon, points2_lat, '.',  c = "blue")
    plt.show()
plot_points()
import sys
lines = sys.stdin.readlines()

points = [0 for i in range(len(lines))]
for i in range(len(lines)):
    points[i] = lines[i].split(',')
    points[i] = [float(coord) for coord in points[i]]

k = int(sys.argv[1])
max_iter = 200
if len(sys.argv) > 2:
    max_iter = int(sys.argv[2])
centroids = [0.0] * k
point_to_cluster = dict()
dim = len(points[0])

def initialize():
    for i in range(k):
        centroids[i] = points[i]
    for i in range(len(points)):
        point_to_cluster[tuple(points[i])] = -1


def distance(p,cluster):
    d = 0
    for i in range(dim):
        d += (p[i] - centroids[cluster][i]) * (p[i] - centroids[cluster][i])
    d = d ** 0.5
    return d


def set_cluster(p):
    min_index = 0
    distances = [0]*k
    for i in range(k):
        distances[i] = distance(p,i)
        if distances[i] < distances[min_index]:
            min_index = i
    point_to_cluster[tuple(p)] = min_index

def update_centroids():
    cluster_to_points = dict()
    for p in points:
        cluster = point_to_cluster[tuple(p)]
        if cluster in cluster_to_points:
            cluster_to_points[cluster].append(p)
        else:
            cluster_to_points[cluster] = [p]

    changed = False
    for i in range(k):
        new_centroid = cluster_mean(i,cluster_to_points)
        if (centroids[i] != new_centroid):
            centroids[i] = new_centroid
            changed = True
    return changed



def cluster_mean(cluster,c2p):
    center = [0.0]*dim
    val = 0
    for i in range(dim):
        for point in c2p[cluster]:
            val += point[i]
        center[i] = val/len(c2p[cluster])
        val = 0
    return center


def k_means_clustering():
    initialize()
    iters = 0
    while True:
        for p in points:
            set_cluster(p)
        changed = update_centroids()
        iters += 1
        if not changed:
            break
        if iters == max_iter:
            break


k_means_clustering()

for i in range(k):
    for j in range(dim):
        centroids[i][j] = '{0:.4f}'.format(centroids[i][j])


for cent in centroids:
    to_print = ','.join(cent)
    print(to_print)
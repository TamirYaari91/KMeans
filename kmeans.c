#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


double distance(double[], double *, int, int);

void set_cluster(int, int, double *, double *, int);

double *cluster_mean(int, int *, double *, int, int);

int update_centroids(int, int, double *, double *, int);

int equal(double *, double *, int);

int main(int argc, char *argv[]) {
    double n1, k_double;
    char c;
    int changed;
    double *centroids;
    double *points_to_cluster;
    int i, j, k, iters, max_iter;
    int dim = 0;
    int num_of_lines = 0;
    char* end_k;
    char* end_max_iter;


    k = strtol(argv[1],&end_k,10);
    if (*end_k != '\0' || k <= 0) {
        printf("Error in k argument!");
        return 0;
    }
    if (argc > 2) {
        max_iter = strtol(argv[2],&end_max_iter,10);
        if (*end_max_iter != '\0' || max_iter <= 0) {
            printf("Error in max_iter argument!");
            return 0;
        }
    } else {
        max_iter = 200;
    }
    while (scanf("%lf%c", &n1, &c) == 2) {
        if (num_of_lines == 0) {
            dim++;
        }
        if (c == '\n') {
            num_of_lines++;
        }
    }
    if (k >= num_of_lines) {
        printf("Error in k argument!");
        return 0;
    }
    centroids = (double *) malloc(k * (dim + 1) * sizeof(double));
    assert(centroids);
    points_to_cluster = (double *) malloc(num_of_lines * (dim + 1) * sizeof(double));
    assert(points_to_cluster);
    for (i = 0; i < num_of_lines * (dim + 1); ++i) {
        points_to_cluster[i] = 0.0;
    }
    i = 0;
    rewind(stdin);
    while (scanf("%lf%c", &n1, &c) == 2) {
        points_to_cluster[i] = n1;
        i++;
        if (c == '\n') {
            points_to_cluster[i] = -1;
            i++;
        }
    }
    j = 0;
    for (i = 1; i < k * (dim + 1) + 1; ++i) {
        if (i % (dim + 1) == 0) {
            continue;
        }
        centroids[j] = points_to_cluster[i - 1];
        j++;
    }
    iters = 0;
    while (1) {
        for (i = 0; i < num_of_lines; ++i) {
            set_cluster(i, k, points_to_cluster, centroids, dim);
        }
        changed = update_centroids(k, num_of_lines, points_to_cluster, centroids, dim);
        iters++;
        if (changed == 0 || iters == max_iter) {
            break;
        }
    }

    for (i = 0; i < (k * dim); ++i) {
        printf("%.4f", centroids[i]);
        if ((i + 1) % dim != 0) {
            printf(",");
        } else {
            printf("\n");
        }
    }
    free(centroids);
    free(points_to_cluster);
    return 0;
}

double distance(double *p, double *centroids, int cluster, int dim) {
    double d = 0;
    int i;
    for (i = 0; i < dim; ++i) {
        double multi;
        multi = (p[i] - *((centroids + cluster * dim) + i));
        d += multi * multi;
    }
    return d;
}

void set_cluster(int p_index, int k, double *point_to_cluster, double *centroids, int dim) {
    int i;
    int min_index = 0;
    double *distances;
    distances = (double *) malloc(k * sizeof(double));
    assert(distances);

    for (i = 0; i < k; ++i) {
        distances[i] = distance((point_to_cluster + p_index * (dim + 1)), centroids, i, dim);
        if (distances[i] < distances[min_index]) {
            min_index = i;
        }
    }
    point_to_cluster[(p_index * (dim + 1) + dim)] = min_index;
}

double *cluster_mean(int cluster, int *c2p, double *p2c, int dim, int num_of_points) {
    int size;
    double val;
    double p_val;
    int i, j;
    static double *center;
    center = (double *) malloc(dim * sizeof(double));
    assert(center);
    size = 0;
    val = 0.0;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < num_of_points; ++j) {
            int p_index = *(c2p + cluster * num_of_points + j);
            if (p_index < 0) {
                break;
            }
            size++;
            p_val = *(p2c + p_index * (dim + 1) + i);
            val += p_val;
        }
        center[i] = val / size;
        size = 0;
        val = 0.0;
    }
    return center;
}

int update_centroids(int k, int num_of_points, double *p2c, double *centroids, int dim) {
    int *c2p;
    int i, j, changed;
    c2p = (int *) malloc(k * num_of_points * sizeof(int));
    assert(c2p);

    for (i = 0; i < k; ++i) {
        for (j = 0; j < num_of_points; ++j) {
            c2p[i * num_of_points + j] = -1;
        }
    }

    for (i = 0; i < num_of_points; ++i) {
        int cluster = p2c[i * (dim + 1) + dim];
        for (j = 0; j < num_of_points; ++j) {
            if (c2p[(cluster * num_of_points) + j] == -1) {
                c2p[(cluster * num_of_points) + j] = i;

                break;
            }
        }
    }
    changed = 0;
    for (i = 0; i < k; ++i) {
        double *new_centroid = cluster_mean(i, c2p, p2c, dim, num_of_points);
        if (equal((centroids + dim * i), new_centroid, dim) == 0) {
            for (j = 0; j < dim; ++j) {
                *(centroids + dim * i + j) = new_centroid[j];
            }
            changed = 1;
        }
        free(new_centroid);
    }
    free(c2p);
    return changed;
}

int equal(double *arr1, double *arr2, int dim) {
    int i;
    for (i = 0; i < dim; ++i) {
        if (arr1[i] != arr2[i]) {
            return 0;
        }
    }
    return 1;
}

# Import modules
import os
import subprocess
import numpy as np
import scipy.sparse as sp
import matplotlib.ticker as tkr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pymeshlab as ml
import igl
import fcmatch
from multiprocessing import Process
from shutil import rmtree
from textwrap import wrap
from time import time as ptime
from psutil import cpu_count
from copy import deepcopy
import robust_laplacian

NEIGHBOURS = 2
DEFAULTKALOITER = 10
BORDERWIDTH = 6
BOXPLOTFIXLIM = True
ROOT = '/home/federico/Work/estimators/generated/'
FIGPATH = '/home/federico/Work/estimators/tests/figures/'
TIMEPATH = '/home/federico/Work/estimators/tests/times/'
OUTLPATH = '/home/federico/Work/estimators/tests/outliers/'
PRECPATH = '/home/federico/Work/estimators/tests/precomputed/'
FIGSIZE = (10, 10)
DRAWOUTLIERS = True
os.chdir(ROOT)


def inspect_values(*arrays, n=5):
    ''' Function for sanity check: print random pairs of values from several arrays
    Input:
        arrays: any number of arrays
        n: number of pairs to pick for each array
    Output: None, prints to stdout
    '''
    print()
    indices = np.random.randint(0, len(arrays[0]), size=(n))
    indices = sorted(indices)
    for a in arrays:
        assert len(arrays[0]) == len(a)
    for i in indices:
        print(i, end=' ')
        for a in arrays:
            print(a[i], end='\t')
        print()


def scalar_err(gt, est):
    ''' computes error on scalar data - currently using SMAPE
    Input:
        gt: ground truth array
        est: estimated array
    Output: array of errors
    '''
    den = np.abs(gt) + np.abs(est)
    den[den == 0] = 1
    return np.abs(gt-est) / den


def vector_norm(t):
    return np.sqrt(np.sum(np.square(t), axis=1))


def vector_err(gt, est):
    ''' computes error on vector data - currently using 3D angles
    Input:
        gt: ground truth array
        est: estimated array
    Output: array of errors
    '''
    ng = np.linalg.norm(gt, axis=1)
    ne = np.linalg.norm(est, axis=1)
    den = ng*ne
    den[den == 0] = np.finfo(float).eps
    cos = np.sum(gt*est, axis=1) / den
    theta = np.arccos(np.clip(cos, -1, 1))
    return theta / np.pi


def get_outliers(gt, est, vector=False, fname='outliers'):
    ''' writes some info about the outliers to a file for debugging
    Input:
        gt: ground truth array
        est: estimated array
        vector: whether this is vector or scalar data
        fname: file name
    Output: None, writes to a file
    '''
    if vector:
        err = vector_norm(gt-est)
    else:
        err = np.abs(gt - est)
    aserr = np.argsort(-err)
    serr = err[aserr]
    q1, q3 = np.percentile(serr, [25, 75])
    iqr = q3 - q1
    threshold = q3 + iqr*1.5
    with open(fname, 'w') as file:
        for i in aserr:
            if (err[i] > threshold):
                file.write(f'{i}\t{gt[i]}\t{est[i]}\t{err[i]}\n')


def log_time(filename, entryname, start, end):
    path = TIMEPATH + filename
    with open(path, 'a+') as file:
        file.write(str(entryname) + ' ')
        file.write(str(end-start))
        file.write('\n')


def print_results(results, parameters):
    for k in results.keys():
        er = '%s' % float('%.1g' % results[k])
        print(f'Method {k} (par {parameters[k]}): {er}')

# Scalar error measures
# ERR_EPSILON = .01


# def rmse(ground, est):
#     return np.average(np.square(ground - est)) ** .5


# def mape(ground, est):
#     return np.average(np.abs((ground-est)/(ground + ERR_EPSILON)))


# def smape(ground, est):
#     return 2 * np.average(np.abs(ground-est) /
#                           (np.abs(ground) + np.abs(est) + ERR_EPSILON))


# # Vector error measures
# def rmse_v(ground, est):
#     return np.average(np.sum(np.square(ground - est), axis=1)) ** .5


# def mape_v(ground, est):
#     return np.average(vector_norm(ground-est) /
#                       (vector_norm(ground) + ERR_EPSILON))


# def smape_v(ground, est):
#     return 2 * np.average(vector_norm(ground-est) /
#                           (vector_norm(ground) + vector_norm(est)
#                            + ERR_EPSILON))


def write_time(filename, value):
    with open(filename, 'w')as file:
        file.write(value)


# Format number of faces/vertices
def format_vf_number(n):
    s = ''
    # millions
    if n > 999999:
        n /= 1000000
        s = 'M'
    # thousands
    elif n > 999:
        n /= 1000
        s = 'K'
    n = round(n)
    return str(n) + s


# BoxPlot
def boxplot(datadict, foldname, plotname, outliers):
    data = {}
    for k, v in datadict.items():
        data[k] = v[~v.mask]

    s = ''
    if outliers:
        s = '.'
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_yscale('linear')
    # ax.set_yscale('function', functions=(lambda x: np.log(x+1), lambda x: np.exp(x)-1))

    if BOXPLOTFIXLIM:
        toplim = max([np.percentile(d, 95) for d in data.values()])
        margin = toplim / 50
        ax.set_ylim(bottom=-margin, top=toplim)

    plt.title(plotname)
    ax.boxplot(data.values(), sym=s, whis=[5, 95])

    if BOXPLOTFIXLIM:
        ax.yaxis.set_minor_locator(tkr.AutoMinorLocator(10))
        tix = ax.get_yticks()
        newtix = [0] + [x for x in tix if x > 0]
        ax.set_yticks(newtix)
        ax.tick_params(axis='x', bottom=False)
        ax.set_xticklabels(data.keys())

    ax.grid(axis='y')
    fig.savefig(FIGPATH + foldname + '/' + plotname,
                format='png', bbox_inches='tight', pad_inches=.1)
    plt.close()

# Bar chart
def barchart(data, foldname, plotname):
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_yscale('linear')
    # ax.set_yscale('function', functions=(lambda x: np.log(x+1), lambda x: np.exp(x)-1))

    plt.title(plotname)
    ax.bar(data.keys(), data.values())
    ax.tick_params(axis='x', bottom=False)
    ax.set_xticklabels(data.keys())

    ax.grid(axis='y')
    fig.savefig(FIGPATH + foldname + '/' + plotname,
                format='png', bbox_inches='tight', pad_inches=.1)
    plt.close()

# Generate border mask (irrelevant to point clouds)
def border_mask(verts, faces, width=1):
    # if there are no faces, there is no border
    if len(faces == 0): return np.full((len(verts),), False, dtype='bool')
    # original border
    mask = np.array(igl.is_border_vertex(verts, faces), dtype='bool')
    for i in range(1, width):
        # expand
        new_mask = np.copy(mask)
        for f in faces:
            border = False
            for v in f:
                border = (border or mask[v])
            for v in f:
                new_mask[v] = (border or new_mask[v])
        mask = new_mask
    return mask


def multiLinePlot(data, foldname, plotname = "", parameters = None):
    # data is a dict where keys are x levels, and values are dictionaries
    # where keys are method names and values are error lists
    x = list(data.keys())
    xs = [str(q) for q in x]
    hasLabels = parameters != None
    labelColours = []

    fig, ax = plt.subplots(figsize=FIGSIZE)
    plt.title(plotname)
    ax.grid(axis='y')
    methods = list(data[x[0]].keys())
    ymid = [[] for m in methods]
    ytop = [[] for m in methods]
    ybot = [[] for m in methods]
    for k in data.keys():
        for i,m in enumerate(methods):
            y = data[k][m]
            avg = np.average(y)
            var = np.var(y)
            ymid[i].append(avg)
            ytop[i].append(avg + var)
            ybot[i].append(avg - var)
    mtot = len(methods)
    colors = cm.get_cmap('tab10')
    for i,m in enumerate(methods):
        cc = colors(i/mtot,1)
        labelColours.append(cc)
        ax.plot(xs, ymid[i], label=m, color=cc)
        ax.fill_between(xs, ybot[i], ytop[i], alpha=.1, color=cc)
    ax.set_ylim(bottom=0, top=max([max(yt) for yt in ytop]))
    # ax.legend()
    ax.grid(which='major', axis='x', linestyle='--')

    if hasLabels:
        pdata = []
        for i in methods:
            aux = []
            for j in x:
                par = parameters[j][i]
                val = str(par) if (par != 0) else '-'
                aux.append(val)
            pdata.append(aux)
        plt.subplots_adjust(bottom=0.3)
        tab = plt.table(cellText = pdata, rowLabels=methods,
            rowColours=labelColours, colLabels=xs, loc='bottom')

        plt.tick_params(axis='x', which='both', bottom=False,
            top=False, labelbottom=False)

        plt.xlim(0-.5,len(x)-.5)


    fig.savefig(FIGPATH + foldname + '/' + plotname,
                format='png', bbox_inches='tight', pad_inches=.1)
    plt.close()


def mergeGraphs(title, fnames, remove=True):
    MERGECMD = 'montage'
    subprocess.call([MERGECMD] + fnames + [
        '-geometry', '+5+5',
        title])
    if remove:
        for f in fnames:
            os.remove(f)


# Mesh class

class My_Mesh:
    name = ""
    path = ""

    vNum = 0
    fNum = 0

    positions = []
    normals = []
    gauss_curv = []
    mean_curv = []
    faces = []

    b_mask = 0   # mask for border vertices

    scalarPath = ""
    scalarField = []

    ### Read mesh from PLY file (path relative to ROOT) ##
    def __init__(self, relpath):
        self.positions = []
        self.normals = []
        self.gauss_curv = []
        self.mean_curv = []
        self.faces = []
        self.relpath = relpath
        self.name = os.path.basename(self.relpath)
        self.abspath = os.path.join(ROOT, relpath + ".ply")

        with open(self.abspath, 'r') as f:
            line = f.readline().strip()
            while line[:] != "end_header":
                line = f.readline().strip()
                s = line.split()
                if s and s[0] == "element":
                    if s[1] == "vertex":
                        self.vNum = int(s[2])
                    elif s[1] == "face":
                        self.fNum = int(s[2])
            for v in range(self.vNum):
                line = f.readline()
                tokens = [float(t) for t in line.split()]
                self.positions.append(tokens[0:3])
                self.normals.append(tokens[3:6])
                # self.parametric.append(tokens[6:8])
                self.gauss_curv.append(tokens[8])
                self.mean_curv.append(tokens[9])
            for e in range(self.fNum):
                line = f.readline()
                tokens = [float(t) for t in line.split()]
                self.faces.append(tokens[1:4])

        self.positions = np.array(self.positions, dtype='float64')
        self.faces = np.array(self.faces, dtype='int32')
        self.b_mask = border_mask(self.positions, self.faces, BORDERWIDTH)

        self.normals = np.ma.array(self.normals,
                                   mask=np.repeat(
                                       self.b_mask[:, np.newaxis], 3, axis=1),
                                   dtype='float64')
        self.gauss_curv = np.ma.array(self.gauss_curv, mask=self.b_mask,
                                      dtype='float64')
        self.mean_curv = np.ma.array(self.mean_curv, mask=self.b_mask,
                                     dtype='float64')

    def is_point_cloud(self):
        return (len(self.faces) == 0)

    def read_scalar_field(self, suffix):
        self.scalarField = []
        self.scalarPath = os.path.join(ROOT, self.relpath + suffix)

        with open(self.scalarPath, 'r') as f:
            # ignore header if present
            starting_pos = f.tell()
            line = f.readline().strip()
            if len(line.split()) == 1:
                # if there's no header, put the line back to the start
                f.seek(starting_pos)

            for v in range(self.vNum):
                line = f.readline().strip()
                value = float(line)
                self.scalarField.append(value)

            self.scalarField = np.ma.array(self.scalarField, mask=self.b_mask,
                                           dtype='float64')

    def average_edge_length(self, K=6):
        if self.is_point_cloud():
            assert(K>0)
            # get knn
            v = self.positions
            pair_dist = igl.all_pairs_distances(v, v, False)
            knn = []
            N = len(v)
            for pt in range(N):
                knn.append([])
                for nb in range(N):
                    d = pair_dist[pt,nb]
                    for i in range(K):
                        if len(knn[-1]) <= i:
                            knn[-1].append((nb, pair_dist[pt,nb]))
                            break
                        else:
                            if 0 < d < knn[-1][i][1]:
                                knn[-1].insert(i, (nb, pair_dist[pt,nb]))
                                break
            assert(len(knn[pt]) >= K)
            return sum([p[1] for p in knn[pt][:K] for pt in range(N)])/(K*N)
        else:
            tot = 0
            cnt = 0
            dist = lambda x,y: np.linalg.norm(x-y)
            v = self.positions
            for f in self.faces:
                tot += dist(v[f[0]], v[f[1]])
                tot += dist(v[f[1]], v[f[2]])
                tot += dist(v[f[2]], v[f[0]])
                cnt += 3
            return tot/cnt

    def write(self):
        with open(self.abspath, 'w') as f:
            f.write('ply\n')
            f.write('format ascii 1.0\n')
            f.write(f'comment {self.name}\n')
            f.write(f'element vertex {len(self.positions)}\n')
            f.write('property double x\n')
            f.write('property double y\n')
            f.write('property double z\n')
            f.write('property double nx\n')
            f.write('property double ny\n')
            f.write('property double nz\n')
            f.write(f'element face {len(self.faces)}\n')
            f.write('property list uchar int vertex_indices\n')
            f.write('end_header\n')
            for v in range(len(self.positions)):
                n0 = self.normals[v][0] if not self.normals.mask[v][0] else 0
                n1 = self.normals[v][1] if not self.normals.mask[v][1] else 0
                n2 = self.normals[v][2] if not self.normals.mask[v][2] else 0
                f.write(str(self.positions[v][0]) + ' ' + str(self.positions[v][1]) + ' '
                    + str(self.positions[v][2]) + ' ' + str(n0) + ' '
                    + str(n1) + ' ' + str(n2) + '\n')
            for e in range(len(self.faces)):
                f.write('3 ' + str(self.faces[e][0]) + ' ' + str(self.faces[e][1])
                    + ' ' + str(self.faces[e][2]) + '\n')
    
    def writeOBJ(self):
        with open(self.abspath[:-3] + 'obj', 'w') as f:
            for v in range(len(self.positions)):
                f.write('v ' + str(self.positions[v][0]) + ' ' + str(self.positions[v][1]) + ' '
                    + str(self.positions[v][2]) + '\n')
            for e in range(len(self.faces)):
                f.write('f ' + str(self.faces[e][0]+1) + ' ' + str(self.faces[e][1]+1)
                    + ' ' + str(self.faces[e][2]+1) + '\n')
    
    def writeOFF(self):
        with open(self.abspath[:-3] + 'off', 'w') as f:
            f.write(f'OFF {self.vNum} {self.fNum} 0\n')
            for v in range(len(self.positions)):
                f.write(str(self.positions[v][0]) + ' ' + str(self.positions[v][1]) + ' '
                    + str(self.positions[v][2]) + '\n')
            for e in range(len(self.faces)):
                f.write('3 ' + str(self.faces[e][0]) + ' ' + str(self.faces[e][1])
                    + ' ' + str(self.faces[e][2]) + '\n')

    def add_noise(self, percentage):
        new_mesh = deepcopy(self)
        ael = new_mesh.average_edge_length()
        epsilon = np.random.normal(0, percentage*ael/100,
            size=new_mesh.positions.shape[0])
        new_mesh.positions += epsilon[:,None] * new_mesh.normals
        new_mesh.name = new_mesh.name + '_n' + str(percentage)
        new_mesh.relpath = new_mesh.relpath + '_n' + str(percentage)
        new_mesh.abspath = os.path.join(ROOT, new_mesh.relpath + ".ply")
        new_mesh.write()
        new_mesh.writeOBJ()
        new_mesh.writeOFF()
        return new_mesh


# def import_generic_mesh(path, fine_faces_min, coarse_faces):
#     fcmatch(path, fine_faces_min, coarse_faces, quiet = False)
#     bn, ex = os.path.splitext(p)
#     coarse = My_Mesh(bn + '_coarse' + ex)
#     fine = My_Mesh(bn + '_fine' + ex)
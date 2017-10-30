import numpy as np
from numpy import sin, cos, pi
from skimage import measure
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ctypes import sizeof, c_float

import struct
import gc

def save_ply(path, pts, colors=np.array([]), faces=np.array([])):
    """
    Saves a 3D mesh model to a PLY file.

    :param path: A path to the resulting PLY file.
    :param pts: nx3 ndarray
    :param colors: nx3 ndarray
    :param faces: mx3 ndarray
    """
    colors = np.array(colors)
    if colors.size != 0:
        assert(len(pts) == len(colors))

    valid_pts_count = 0
    for pt_id, pt in enumerate(pts):
        if not np.isnan(np.sum(pt)):
            valid_pts_count += 1

    f = open(path, 'w')
    f.write(
        'ply\n'
        'format ascii 1.0\n'
        #'format binary_little_endian 1.0\n'
        'element vertex ' + str(valid_pts_count) + '\n'
        'property float x\n'
        'property float y\n'
        'property float z\n'
    )
    if colors.size != 0:
        f.write(
            'property uchar red\n'
            'property uchar green\n'
            'property uchar blue\n'
            'property uchar alpha\n'
        )
    if faces.size != 0:
        f.write(
            'element face ' + str(len(faces)) + '\n'
            'property list uchar int vertex_indices\n'
        )
    f.write('end_header\n')

    for pt_id, pt in enumerate(pts):
        if not np.isnan(np.sum(pt)):
            f.write(' '.join(map(str, pt.squeeze().tolist())) + ' ')
            if colors.size != 0:
                clr = list(colors[pt_id])
                clr.append(1)
                f.write(' '.join(map(str, map(int, clr))))
            f.write('\n')
    for face in faces:
        f.write(' '.join(map(str, map(int, [len(face)] + list(face.squeeze())))) + ' ')
        f.write('\n')
    f.close()


def main():
    ## Load TSDF voxel grid from binary file
    fid = open('tsdf.bin', 'rb')
    buf = fid.read()
    fid.close()
    info = np.fromstring(buf[0:4*8], np.float32)
    voxelGridDimX = int(info[0])
    voxelGridDimY = int(info[1])
    voxelGridDimZ = int(info[2])
    voxel_grid_origin = info[3:6]
    voxel_size = info[6]

    tsdf = np.fromstring(buf[4*8:], np.float32)
    tsdf = tsdf.reshape(voxelGridDimX, voxelGridDimY, voxelGridDimZ)
    print "executing marching cube"
    verts, faces = measure.marching_cubes(tsdf, 0)
    # np.save('verts.npy', verts)
    # np.save('faces.npy', faces)
    # verts = np.load('verts.npy')
    # faces = np.load('faces.npy')
    fid_color = open('tsdf_color.bin', 'rb')
    buf_color = fid_color.read()
    fid_color.close()

    tsdf_color = np.fromstring(buf_color, np.uint8)
    tsdf_color = tsdf_color.reshape(voxelGridDimX, voxelGridDimY, voxelGridDimZ, 3)
    verts_idx = np.round(verts).astype(np.int32)
    colors = tsdf_color[verts_idx[:, 0], verts_idx[:, 1], verts_idx[:, 2], :]
    del tsdf_color
    gc.collect()
    verts[:, 2] = - verts[:, 2]
    verts = (voxel_grid_origin[np.newaxis, :] + verts * voxel_size)
    print "vertices color num : {}".format(len(colors))
    print "vertices num : {}".format(len(verts))
    print "faces num : {}".format(len(faces))
    print "saving ply mesh"
    save_ply('mesh.ply', verts, colors=colors, faces=faces)

if __name__ =='__main__':
    main()

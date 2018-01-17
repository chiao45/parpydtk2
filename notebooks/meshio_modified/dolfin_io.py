# -*- coding: utf-8 -*-
#
'''
I/O for DOLFIN's XML format, cf.
<https://people.sc.fsu.edu/~jburkardt/data/dolfin_xml/dolfin_xml.html>.

.. moduleauthor:: Nico Schlömer <nico.schloemer@gmail.com>
'''
import logging
import os
import re
import xml.etree.cElementTree as ET

import numpy


def _read_mesh(filename):
    dolfin_to_meshio_type = {
        'triangle': ('triangle', 3),
        'tetrahedron': ('tetra', 4),
        }

    # Use iterparse() to avoid loading the entire file via parse(). iterparse()
    # allows to discard elements (via clear()) after they have been processed.
    # See <https://stackoverflow.com/a/326541/353337>.
    for event, elem in ET.iterparse(filename, events=('start', 'end')):
        if event == 'end':
            continue

        if elem.tag == 'dolfin':
            # Don't be too strict with the assertion. Some meshe files don't
            # have the proper tags.
            # assert elem.attrib['nsmap'] \
            #     == '{\'dolfin\': \'https://fenicsproject.org/\'}'
            pass
        elif elem.tag == 'mesh':
            dim = int(elem.attrib['dim'])
            cell_type, npc = dolfin_to_meshio_type[elem.attrib['celltype']]
        elif elem.tag == 'vertices':
            points = numpy.empty((int(elem.attrib['size']), 3))
            keys = ['x', 'y']
            if dim == 2:
                points[:, 2] = 0.0
            else:
                assert dim == 3
                keys += ['z']
        elif elem.tag == 'vertex':
            k = int(elem.attrib['index'])
            points[k][:dim] = [float(elem.attrib[key]) for key in keys]
        elif elem.tag == 'cells':
            cells = {
                cell_type:
                numpy.empty((int(elem.attrib['size']), npc), dtype=int)
                }
        elif elem.tag in ['triangle', 'tetrahedron']:
            k = int(elem.attrib['index'])
            cells[cell_type][k] = [
                int(elem.attrib['v{}'.format(i)]) for i in range(npc)
                ]
        else:
            logging.warning('Unknown entry %s. Ignoring.', elem.tag)

        elem.clear()

    return points, cells, cell_type


def _read_cell_data(filename, cell_type):
    dolfin_type_to_numpy_type = {
        'int': numpy.dtype('int'),
        'float': numpy.dtype('float'),
        'uint': numpy.dtype('uint'),
        }

    cell_data = {cell_type: {}}
    dir_name = os.path.dirname(filename)
    if not os.path.dirname(filename):
        dir_name = os.getcwd()

    # Loop over all files in the same directory as `filename`.
    basename = os.path.splitext(filename)[0]
    for f in os.listdir(dir_name):
        # Check if there are files by the name "<filename>_*.xml"; if yes,
        # extract the * pattern and make it the name of the data set.
        out = re.match('{}_([^\\.]+)\\.xml'.format(basename), f)
        if not out:
            continue
        name = out.group(1)
        tree = ET.parse(f)
        root = tree.getroot()

        mesh_functions = list(root)
        assert len(mesh_functions) == 1
        mesh_function = mesh_functions[0]

        assert mesh_function.tag == 'mesh_function'
        size = int(mesh_function.attrib['size'])
        dtype = dolfin_type_to_numpy_type[mesh_function.attrib['type']]
        data = numpy.empty(size, dtype=dtype)
        for child in mesh_function:
            assert child.tag == 'entity'
            idx = int(child.attrib['index'])
            data[idx] = child.attrib['value']

        cell_data[cell_type][name] = data
    return cell_data


def read(filename):
    points, cells, cell_type = _read_mesh(filename)
    point_data = {}
    cell_data = _read_cell_data(filename, cell_type)
    field_data = {}
    return points, cells, point_data, cell_data, field_data


def _write_mesh(
        filename,
        points,
        cell_type,
        cells
        ):
    stripped_cells = {cell_type: cells[cell_type]}

    dolfin = ET.Element(
        'dolfin',
        nsmap={'dolfin': 'https://fenicsproject.org/'}
        )

    meshio_to_dolfin_type = {
        'triangle': 'triangle',
        'tetra': 'tetrahedron',
        }

    if len(cells) > 1:
        discarded_cells = list(cells.keys())
        discarded_cells.remove(cell_type)
        logging.warning(
          'DOLFIN XML can only handle one cell type at a time. '
          'Using %s, discarding %s.',
          cell_type, ', '.join(discarded_cells)
          )

    dim = 2 if all(points[:, 2] == 0) else 3

    mesh = ET.SubElement(
            dolfin,
            'mesh',
            celltype=meshio_to_dolfin_type[cell_type],
            dim=str(dim)
            )
    vertices = ET.SubElement(mesh, 'vertices', size=str(len(points)))

    for k, point in enumerate(points):
        ET.SubElement(
            vertices,
            'vertex',
            index=str(k),
            x=repr(point[0]),
            y=repr(point[1]),
            z=repr(point[2])
            )

    num_cells = 0
    for cls in stripped_cells.values():
        num_cells += len(cls)

    xcells = ET.SubElement(mesh, 'cells', size=str(num_cells))
    idx = 0
    for ct, cls in stripped_cells.items():
        for cell in cls:
            cell_entry = ET.SubElement(
                xcells,
                meshio_to_dolfin_type[ct],
                index=str(idx)
                )
            for k, c in enumerate(cell):
                cell_entry.attrib['v{}'.format(k)] = str(c)
            idx += 1

    tree = ET.ElementTree(dolfin)
    tree.write(filename)
    return


def _numpy_type_to_dolfin_type(dtype):
    types = ['int', 'uint', 'float']
    for t in types:
        # issubtype handles all of int8, int16, float64 etc.
        if numpy.issubdtype(dtype, numpy.dtype(t)):
            return t
    return None


def _write_cell_data(
        filename,
        dim,
        cell_data
        ):
    dolfin = ET.Element(
        'dolfin',
        nsmap={'dolfin': 'https://fenicsproject.org/'}
        )

    mesh_function = ET.SubElement(
            dolfin,
            'mesh_function',
            type=_numpy_type_to_dolfin_type(cell_data.dtype),
            dim=str(dim),
            size=str(len(cell_data))
            )

    for k, value in enumerate(cell_data):
        ET.SubElement(
            mesh_function,
            'entity',
            index=str(k),
            value=repr(value),
            )

    tree = ET.ElementTree(dolfin)
    tree.write(filename)
    return


def write(
        filename,
        points,
        cells,
        point_data=None,
        cell_data=None,
        field_data=None
        ):
    logging.warning(
        'Dolfin\'s XML is a legacy format. Consider using XDMF instead.'
        )

    point_data = {} if point_data is None else point_data
    cell_data = {} if cell_data is None else cell_data
    field_data = {} if field_data is None else field_data

    if 'tetra' in cells:
        cell_type = 'tetra'
    else:
        assert 'triangle' in cells
        cell_type = 'triangle'

    _write_mesh(filename, points, cell_type, cells)

    if cell_type in cell_data:
        for key, data in cell_data[cell_type].items():
            cell_data_filename = \
                '{}_{}.xml'.format(os.path.splitext(filename)[0], key)
            dim = 2 if all(points[:, 2] == 0) else 3
            _write_cell_data(cell_data_filename, dim, numpy.array(data))
    return

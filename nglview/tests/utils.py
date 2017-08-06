import os


def get_fn(fn):
    this_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(this_path, 'data', fn)


repr_dict = {
    '0': {
        '0': {
            'params': {
                'aspectRatio': 5,
                'assembly': 'default',
                'capped': True,
                'clipCenter': {
                    'x': 0,
                    'y': 0,
                    'z': 0
                },
                'clipNear': 0,
                'clipRadius': 0,
                'colorMode': 'hcl',
                'colorReverse': False,
                'colorScale': 'RdYlBu',
                'colorScheme': 'chainname',
                'colorValue': 9474192,
                'defaultAssembly': '',
                'depthWrite': True,
                'diffuse': 16777215,
                'disablePicking': False,
                'flatShaded': False,
                'lazy': False,
                'linewidth': 2,
                'matrix': {
                    'elements':
                    [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
                },
                'metalness': 0,
                'opacity': 1,
                'quality': 'high',
                'radialSegments': 20,
                'radius': 'sstruc',
                'roughness': 0.4,
                'scale': 0.7,
                'sele': '',
                'side': 'double',
                'smoothSheet': False,
                'subdiv': 6,
                'tension': None,
                'visible': True,
                'wireframe': False
            },
            'type': 'cartoon'
        },
        '1': {
            'params': {
                'aspectRatio': 1,
                'assembly': 'default',
                'bondScale': 0.4,
                'clipCenter': {
                    'x': 0,
                    'y': 0,
                    'z': 0
                },
                'clipNear': 0,
                'clipRadius': 0,
                'colorMode': 'hcl',
                'colorReverse': False,
                'colorScale': 'RdYlBu',
                'colorScheme': 'chainname',
                'colorValue': 9474192,
                'cylinderOnly': False,
                'defaultAssembly': '',
                'depthWrite': True,
                'diffuse': 16777215,
                'disableImpostor': False,
                'disablePicking': False,
                'flatShaded': False,
                'lazy': False,
                'lineOnly': False,
                'linewidth': 2,
                'matrix': {
                    'elements':
                    [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
                },
                'metalness': 0,
                'opacity': 1,
                'openEnded': True,
                'quality': 'high',
                'radialSegments': 20,
                'radius': 0.3,
                'roughness': 0.4,
                'scale': 1,
                'sele': '',
                'side': 'double',
                'sphereDetail': 2,
                'visible': True,
                'wireframe': False
            },
            'type': 'base'
        },
        '2': {
            'params': {
                'aspectRatio': 1.5,
                'assembly': 'default',
                'bondScale': 0.3,
                'bondSpacing': 0.75,
                'clipCenter': {
                    'x': 0,
                    'y': 0,
                    'z': 0
                },
                'clipNear': 0,
                'clipRadius': 0,
                'colorMode': 'hcl',
                'colorReverse': False,
                'colorScale': '',
                'colorScheme': 'element',
                'colorValue': 9474192,
                'cylinderOnly': False,
                'defaultAssembly': '',
                'depthWrite': True,
                'diffuse': 16777215,
                'disableImpostor': False,
                'disablePicking': False,
                'flatShaded': False,
                'lazy': False,
                'lineOnly': False,
                'linewidth': 2,
                'matrix': {
                    'elements':
                    [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
                },
                'metalness': 0,
                'multipleBond': 'off',
                'opacity': 1,
                'openEnded': True,
                'quality': 'high',
                'radialSegments': 20,
                'radius': 0.15,
                'roughness': 0.4,
                'scale': 2,
                'sele': 'ligand',
                'side': 'double',
                'sphereDetail': 2,
                'visible': True,
                'wireframe': False
            },
            'type': 'ball+stick'
        }
    }
}

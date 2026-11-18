"""Tests for pf.view(), the linked sequence + 3D viewer (the Python side; the
interactions were checked in a headless browser)."""

import json
import re
from pathlib import Path

import pytest

from pdbframe import PDBFrame
from pdbframe.structure.viewer.widget import _residues_of

DATA = Path(__file__).parent / 'data'


@pytest.fixture(scope='module')
def ref():
    return PDBFrame.read(DATA / '3tfk_rcsb.pdb')


@pytest.fixture(scope='module')
def model():
    return PDBFrame.read(DATA / '3tfk_af2_renamed.pdb')


def _data(view):
    return json.loads(re.search(r'var D = (\{.*?\});\n', view.html, re.S).group(1))


def _shown(d, m=0):
    return [c['label'] for c in d['models'][m]['chains'] if c['shown']]


def test_region_specs(ref):
    assert len(_residues_of(ref, 'C')) == 9
    assert len(_residues_of(ref, ('A', 50, 59))) == 10
    pocket = ref.near(ref.select(chains='C'), cutoff=4)
    assert len(_residues_of(ref, pocket)) == len(pocket.residues)
    # a selection's observed_index restarts at 1; regions map back by numbering
    assert ('A', 1) not in _residues_of(ref, ref.select(chains='A', residues=(50, 59)))
    assert len(_residues_of(ref, ['C', ('A', 50, 59)])) == 19


def test_widget_data(ref):
    d = _data(ref.view(regions={'peptide': 'C', 'helix': ('A', 50, 85)},
                       region_colors={'helix': 'orange'}))
    assert len(d['models']) == 1
    assert [c['label'] for c in d['models'][0]['chains']] == ref.chains
    assert sum(len(c['residues']) for c in d['models'][0]['chains']) == 623
    assert [r['name'] for r in d['regions']] == ['peptide', 'helix']
    assert d['regions'][1]['color'] == 'orange'
    assert d['models'][0]['pdb'].startswith('ATOM')


def test_sequence_option(ref):
    assert _shown(_data(ref.view())) == []                       # none by default
    assert _shown(_data(ref.view(sequence=True))) == ref.chains
    assert _shown(_data(ref.view(sequence=['A', 'C']))) == ['A', 'C']
    with pytest.raises(KeyError, match='not in'):
        ref.view(sequence=['Q'])


def test_sequences_only(ref):
    """sequence_only opens the sequences (all chains, or the ones named) and no 3D view."""
    d = _data(ref.view(sequence_only=True))
    assert d['sequence_only'] is True and _shown(d) == ref.chains
    assert _shown(_data(ref.view(sequence_only=True, sequence=['C']))) == ['C']
    assert _data(ref.view())['sequence_only'] is False
    assert 'Show 3D' in ref.view(sequence_only=True).html


def test_long_chain_ids(ref):
    d = _data(ref.rename_chains({'A': 'MHC'}).view(regions={'r': ('MHC', 50, 60)}))
    ch = d['models'][0]['chains'][0]
    assert ch['label'] == 'MHC' and len(ch['id']) == 1


def test_region_styles(ref):
    d = _data(ref.view(regions={'pep': 'C', 'helix': ('A', 50, 60)},
                       region_styles={'pep': 'stick'}))
    assert [r['extra'] for r in d['regions']] == ['stick', None]
    with pytest.raises(KeyError, match='not given'):
        ref.view(regions={'pep': 'C'}, region_styles={'other': 'stick'})


def test_compare_links_residues(ref, model):
    d = _data(ref.view(compare=model))
    assert len(d['models']) == 2
    assert [c['label'] for c in d['models'][1]['chains']] == ['H', 'P', 'X', 'Y']
    # each reference residue links to its partner in the model, and back
    assert len(d['links']) == 2 * 623
    # sequence= names the frame's chains; their partners in the model show too
    d2 = _data(ref.view(compare=model, sequence=['A', 'C']))
    assert _shown(d2, 1) == ['H', 'P']
    first = d['links']['0:A:1']
    assert first == ['1:A:1'] and d['links']['1:A:1'] == ['0:A:1']


def test_values(ref, model):
    lddt = model.metrics.lddt(ref, per_residue=True)
    d = _data(ref.view(values=lddt))
    assert d['values']['name'] == 'lddt' and (d['values']['lo'], d['values']['hi']) == (0, 1)
    res = d['models'][0]['chains'][0]['residues'][0]
    assert res[4].startswith('#') and 0 <= res[5] <= 1
    d = _data(model.view(values='plddt'))
    assert (d['values']['lo'], d['values']['hi']) == (50, 90)
    with pytest.raises(KeyError, match='needs columns'):
        ref.view(values={'x': [1]})


def test_interfaces(ref):
    d = _data(ref.view(interfaces=[('A', 'C')]))
    assert [r['name'] for r in d['regions']] == ['A–C interface']
    chains = {k.split(':')[1] for k in d['regions'][0]['keys']}
    assert chains == {'A', 'B'}       # internal IDs of chains A and C
    assert len(_data(ref.view(interfaces=True))['regions']) == 6


def test_options(ref):
    d = _data(ref.view(surface=True, axes=False, ligands=False))
    assert (d['surface'], d['axes'], d['ligands']) == (True, False, False)


def test_errors(ref):
    with pytest.raises(ValueError, match='no residues'):
        ref.view(regions={'empty': ('A', 900, 950)})
    with pytest.raises(ValueError, match='style'):
        ref.view(style='ribbon')


def test_page_includes_script_and_style(ref, tmp_path):
    v = ref.view()
    assert 'function makeAxes' in v.html and 'function copyAsCode' in v.html
    assert '#ROOT' not in v.html
    assert v.save(tmp_path / 'v.html').read_text().startswith('<!DOCTYPE html>')


# ------------------------------------------------- several structures (roadmap)

import numpy as np

import pdbframe
from pdbframe.structure.viewer import COMPARE_COLOR, STRUCTURE_COLORS


def _model_color(d, m):
    return d['models'][m]['chains'][0]['color']


def test_compare_list_dict_and_names(ref, model):
    d = _data(ref.view(compare=[model, model]))
    assert len(d['models']) == 3
    assert [m['name'] for m in d['models']] == ['3tfk_rcsb', '3tfk_af2_renamed (2)', '3tfk_af2_renamed (3)']
    d = _data(ref.view(compare={'AF2': model, 'AF2b': model}))
    assert [m['name'] for m in d['models']] == ['3tfk_rcsb', 'AF2', 'AF2b']
    d = _data(pdbframe.view([ref, model], names=['crystal', 'af2']))
    assert [m['name'] for m in d['models']] == ['crystal', 'af2']
    with pytest.raises(ValueError, match='names'):
        ref.view(compare=[model], names=['only one'])


def test_compare_colors(ref, model):
    assert _model_color(_data(ref.view(compare=model)), 1) == COMPARE_COLOR      # one: grey, as before
    d = _data(ref.view(compare=[model, model]))
    assert [_model_color(d, m) for m in (1, 2)] == STRUCTURE_COLORS[:2]          # several: one each
    d = _data(ref.view(compare=[model, model], compare_color='black'))
    assert {_model_color(d, m) for m in (1, 2)} == {'black'}
    d = _data(ref.view(compare=[model, model], colors=[{'A': 'red'}, 'blue', None]))
    assert d['models'][0]['chains'][0]['color'] == 'red' and _model_color(d, 1) == 'blue'


def test_per_structure_lists(ref, model):
    lddt = model.metrics.lddt(ref, per_residue=True)
    d = _data(ref.view(compare=[model], style=['cartoon', 'line'],
                       values=[lddt, 'plddt'],
                       regions=[{'pocket': ('A', 50, 60)}, {'pocket': ('H', 50, 60), 'tail': 'P'}]))
    assert [m['style'] for m in d['models']] == ['cartoon', 'line']
    assert [v['name'] for v in d['value_scales']] == ['lddt', 'plddt']
    pocket = [r for r in d['regions'] if r['name'] == 'pocket'][0]
    assert {k.split(':')[0] for k in pocket['keys']} == {'0', '1'}   # one region over both
    with pytest.raises(ValueError, match='2 structures'):
        ref.view(compare=[model], style=['cartoon'])


def test_links_through_the_lead(ref, model):
    d = _data(ref.view(compare=[model, model]))
    assert sorted(d['links']['0:A:1']) == ['1:A:1', '2:A:1']
    assert d['links']['2:A:1'] == ['0:A:1']


def test_aligned_sequence_rows(ref, model):
    # remove residues on each side, so both get gaps in the aligned rows
    lead = ref.select(chains='A', residues=(50, 55), invert=True)
    other = model.select(chains='H', residues=(80, 84), invert=True)
    d = _data(lead.view(compare=other))
    group = d['layout'][0]
    lead_row, other_row = group['rows']
    assert lead_row['chain'] == '0:A' and other_row['chain'].startswith('1:')
    assert len(lead_row['cells']) == len(other_row['cells'])
    assert None in lead_row['cells'] and None in other_row['cells']
    # paired residues share a column
    for a, b in zip(lead_row['cells'], other_row['cells']):
        if a is not None and b is not None:
            assert d['links'][f'0:A:{a}'] == [f"{other_row['chain']}:{b}"]


def test_region_toggles_and_groups(ref):
    d = _data(ref.view(regions={'h1': ('A', 50, 60), 'h2': ('A', 140, 150), 'pep': 'C'},
                       region_groups={'helices': ['h1', 'h2']}, hidden_regions=['helices']))
    assert d['region_groups'] == {'helices': ['h1', 'h2']}
    assert d['hidden_regions'] == ['h1', 'h2']
    with pytest.raises(KeyError, match='region_groups'):
        ref.view(regions={'h1': 'A'}, region_groups={'g': ['nope']})
    with pytest.raises(KeyError, match='hidden_regions'):
        ref.view(regions={'h1': 'A'}, hidden_regions=['nope'])


def test_orientations(ref):
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)
    d = _data(ref.view(orientations={'up': R, 'fit': (R, [1, 2, 3])}, orientation_center='frame'))
    assert d['orientations']['up'] == [R.tolist(), [0.0, 0.0, 0.0]]
    assert d['orientations']['fit'][1] == [1.0, 2.0, 3.0] and d['orientation_center'] == 'frame'
    with pytest.raises(ValueError, match='proper rotation'):
        ref.view(orientations={'bad': np.diag([1, 1, -1])})


def test_start_orientation(ref):
    R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    assert _data(ref.view())['orientation'] is None
    assert _data(ref.view(orientation=R))['orientation'] == [R.tolist(), [0.0, 0.0, 0.0]]
    d = _data(ref.view(orientations={'up': (R, [1, 2, 3])}, orientation='up'))
    assert d['orientation'] == [R.tolist(), [1.0, 2.0, 3.0]]
    with pytest.raises(ValueError, match='not one of the orientations'):
        ref.view(orientations={'up': R}, orientation='down')
    # the coordinates are not changed
    assert _data(ref.view(orientation=R))['models'] == _data(ref.view())['models']
    with pytest.raises(ValueError, match='orientation_center'):
        ref.view(orientation_center='middle')


def test_compact(ref):
    d = _data(ref.view(compact=True))
    atoms = {line[12:16].strip() for line in d['models'][0]['pdb'].splitlines() if line.startswith('ATOM')}
    assert atoms == {'N', 'CA', 'C', 'O'} and d['ligands'] is False


def test_values_must_match_the_structure(ref, model):
    lddt = model.metrics.lddt(ref, per_residue=True)      # keyed by the reference's chains
    with pytest.raises(ValueError, match='matches no residues'):
        ref.view(compare=[model], values=[None, lddt])


def _first_atom(d, m):
    line = next(l for l in d['models'][m]['pdb'].splitlines() if l.startswith('ATOM'))
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def test_superpose_option(ref, model):
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)
    moved = ref.transform(R, [20, 0, 0])                       # a rigidly moved copy
    as_is = _data(ref.view(compare=moved))
    assert np.linalg.norm(_first_atom(as_is, 1) - _first_atom(as_is, 0)) > 5   # default: untouched
    for option in (True, 'all'):
        fitted = _data(ref.view(compare=moved, superpose=option))
        np.testing.assert_allclose(_first_atom(fitted, 1), _first_atom(fitted, 0), atol=1e-2)
    on_a1 = _data(ref.view(compare=moved, superpose='A'))           # one chain ID
    np.testing.assert_allclose(_first_atom(on_a1, 1), _first_atom(on_a1, 0), atol=1e-2)
    assert ref.df.loc[0, 'x'] != moved.df.loc[0, 'x']           # the frames passed in are unchanged
    # fit on the reference's chain A only; the model's chain there is H
    on_a = _data(ref.view(compare=model, superpose=['A']))
    assert on_a['models'][1]['chains'][0]['label'] == 'H'
    with pytest.raises(KeyError, match='superpose'):
        ref.view(compare=model, superpose=['Q'])


def test_arrows(ref):
    ca = ref.coords[(ref.df['atom_name'] == 'CA').to_numpy()]
    d = _data(ref.view(arrows=[(ca[0], ca[1]),
                               {'start': ca[2], 'end': ca[3], 'color': 'red', 'label': 'axis',
                                'group': 'Axes', 'radius': 0.2}]))
    first, second = d['arrows']
    assert first['start'] == pytest.approx(ca[0].tolist()) and first['group'] == 'Arrows'
    assert first['label'] is None and first['color'] == '#000000'
    assert second == {'start': pytest.approx(ca[2].tolist()), 'end': pytest.approx(ca[3].tolist()),
                      'color': 'red', 'label': 'axis', 'radius': 0.2, 'group': 'Axes',
                      'head': True, 'info': None}
    assert _data(ref.view())['arrows'] == []
    assert len(_data(ref.view(arrows={'start': ca[0], 'end': ca[1]}))['arrows']) == 1
    with pytest.raises(ValueError, match='different points'):
        ref.view(arrows=[(ca[0], ca[0])])
    with pytest.raises(ValueError, match='3 numbers'):
        ref.view(arrows=[([1, 2], ca[0])])
    with pytest.raises(KeyError, match='unknown keys'):
        ref.view(arrows=[{'start': ca[0], 'end': ca[1], 'colour': 'red'}])


def test_arrows_without_a_head(ref):
    """A contact or a bond has no direction, so it is drawn as a plain cylinder."""
    ca = ref.coords[(ref.df['atom_name'] == 'CA').to_numpy()]
    d = _data(ref.view(arrows=[{'start': ca[0], 'end': ca[1], 'head': False, 'info': '3.4 A'}]))
    assert d['arrows'][0]['head'] is False and d['arrows'][0]['info'] == '3.4 A'


def test_spheres(ref):
    ca = ref.coords[(ref.df['atom_name'] == 'CA').to_numpy()]
    d = _data(ref.view(spheres=[(ca[0],),
                                {'center': ca[1], 'radius': 0.9, 'color': 'red',
                                 'label': 'P1', 'group': 'Nodes', 'info': 'chain C 1'}]))
    first, second = d['spheres']
    assert first['center'] == pytest.approx(ca[0].tolist())
    assert (first['radius'], first['group'], first['label']) == (0.6, 'Spheres', None)
    assert second == {'center': pytest.approx(ca[1].tolist()), 'radius': 0.9, 'color': 'red',
                      'label': 'P1', 'group': 'Nodes', 'info': 'chain C 1'}
    assert _data(ref.view())['spheres'] == []
    assert len(_data(ref.view(spheres={'center': ca[0]}))['spheres']) == 1
    with pytest.raises(ValueError, match='3 numbers'):
        ref.view(spheres=[([1, 2],)])
    with pytest.raises(KeyError, match='unknown keys'):
        ref.view(spheres=[{'center': ca[0], 'colour': 'red'}])


def test_opacity(ref):
    assert _data(ref.view())['opacity'] == 1.0
    assert _data(ref.view(opacity=0.4))['opacity'] == 0.4
    with pytest.raises(ValueError, match='between 0 and 1'):
        ref.view(opacity=1.5)

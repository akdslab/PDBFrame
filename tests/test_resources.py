"""The documentation and skills that come with the package: docs, install_skills, CLI."""

import tomllib
from pathlib import Path

import pytest

import pdbframe
from pdbframe.cli import main
from pdbframe.resources import SKILLS, skills

ROOT = Path(__file__).parents[1]


def test_docs():
    pages = pdbframe.docs()
    assert {'pdbframe', 'seqframe'} <= set(pages)
    assert pdbframe.docs('seqframe.md') == pdbframe.docs('seqframe') == pages['seqframe']
    with pytest.raises(FileNotFoundError, match='the pages are .*pdbframe'):
        pdbframe.docs('nowhere')


def test_skills_match_the_wheel():
    found = skills()
    assert list(found) == SKILLS
    for name, folder in found.items():
        assert f'name: {name}\n' in (folder / 'SKILL.md').read_text()
    config = tomllib.loads((ROOT / 'pyproject.toml').read_text())
    packaged = config['tool']['hatch']['build']['targets']['wheel']['force-include']
    assert packaged.pop('docs') == 'pdbframe/docs'
    assert packaged == {f'.claude/skills/{s}': f'pdbframe/skills/{s}' for s in SKILLS}


def test_install_skills(tmp_path):
    assert set(pdbframe.install_skills(tmp_path).values()) == {'installed'}
    assert set(pdbframe.install_skills(tmp_path).values()) == {'up to date'}
    edited = tmp_path / SKILLS[0] / 'SKILL.md'
    edited.write_text('my own version')
    assert pdbframe.install_skills(tmp_path)[SKILLS[0]].startswith('kept')
    assert pdbframe.install_skills(tmp_path, force=True)[SKILLS[0]] == 'updated'
    assert edited.read_text() == (skills()[SKILLS[0]] / 'SKILL.md').read_text()


def test_cli(tmp_path, capsys):
    assert main(['docs']) == 0
    assert 'seqframe' in capsys.readouterr().out
    assert main(['docs', 'seqframe', '--path']) == 0
    assert capsys.readouterr().out.strip() == str(pdbframe.docs('seqframe'))
    with pytest.raises(SystemExit, match='no page'):
        main(['docs', 'nowhere'])
    assert main(['skills', 'install', '--target', str(tmp_path)]) == 0
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(SKILLS)

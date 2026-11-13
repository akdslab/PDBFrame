"""Command line: the documentation and the Claude Code skills that come with pdbframe.

    pdbframe docs                     the pages, with their files
    pdbframe docs seqframe            print one (--path: its file)
    pdbframe skills install           copy the skills into ~/.claude/skills
    pdbframe skills install --project ... or into ./.claude/skills
"""

from __future__ import annotations

import argparse
import sys

__all__ = ['main']


def _docs(args):
    from .resources import docs
    try:
        page = None if args.name is None else docs(args.name)
    except FileNotFoundError as e:
        raise SystemExit(str(e))
    if page is None:
        for name, path in docs().items():
            print(f'{name:<12}{path}')
    elif args.path:
        print(page)
    else:
        print(page.read_text())


def _skills(args):
    from .resources import install_skills, skills
    if args.action == 'list':
        for name, path in skills().items():
            print(f'{name:<22}{path}')
        return
    for name, status in install_skills(args.target, project=args.project,
                                       force=args.force).items():
        print(f'{name:<22}{status}')


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog='pdbframe', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    q = sub.add_parser('docs', help='the documentation pages',
                       description='Without a name, lists the pages and their files; with '
                                   'one, prints it.')
    q.add_argument('name', nargs='?', help='a page: pdbframe or seqframe')
    q.add_argument('--path', action='store_true', help="print the page's file, not its text")
    q.set_defaults(func=_docs)
    q = sub.add_parser('skills', help='the Claude Code skills for this package',
                       description='install copies them into ~/.claude/skills (every '
                                   'project) or, with --project, ./.claude/skills.')
    q.add_argument('action', choices=['install', 'list'])
    q.add_argument('--project', action='store_true', help='into ./.claude/skills')
    q.add_argument('--target', help='another skills folder to copy into')
    q.add_argument('--force', action='store_true', help='replace skills that differ')
    q.set_defaults(func=_skills)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == '__main__':
    sys.exit(main())

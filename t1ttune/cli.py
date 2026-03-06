#! /usr/bin/env python3 

import argparse
import sys

from .scripts.__init__ import registry


def main():
    # Top-level parser — mostra i sottocomandi disponibili nell'help
    parser = argparse.ArgumentParser(
        prog="t1ttune",
        description="Main program. Select a subcommand to start.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_build_commands_epilog(),
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<command>",
    )

    # Registra dinamicamente ogni sottocomando trovato nel registry
    for name, cmd in registry.items():
        sub = subparsers.add_parser(
            name,
            help=cmd.SHORT_HELP,           # mostrato in: M --help
            description=cmd.DESCRIPTION,   # mostrato in: M T1 --help
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        cmd.add_arguments(sub)             # ogni comando aggiunge i propri argomenti

    # Se invocato senza sottocomando → mostra help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    registry[args.command].run(args)


def _build_commands_epilog() -> str:
    """Costruisce la sezione descrittiva dei comandi per l'help globale."""
    if not registry:
        return ""
    lines = ["Available commands:"]
    for name, cmd in registry.items():
        lines.append(f"  {name:<12} {cmd.SHORT_HELP}")
    lines.append("")
    lines.append("Use 't1ttune <command> --help' for details on each command.")
    return "\n".join(lines)

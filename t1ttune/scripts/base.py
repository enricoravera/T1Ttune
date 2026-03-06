#! /usr/bin/env python3

import argparse
from abc import ABC, abstractmethod


class BaseCommand(ABC):
    # Testo breve mostrato nell'elenco di `M --help`
    SHORT_HELP: str = "(nessuna descrizione)"

    # Testo esteso mostrato in `M <comando> --help`
    DESCRIPTION: str = ""

    @staticmethod
    @abstractmethod
    def add_arguments(parser: argparse.ArgumentParser) -> None:
        """Aggiunge gli argomenti specifici di questo comando al parser."""
        ...

    @staticmethod
    @abstractmethod
    def run(args: argparse.Namespace) -> None:
        """Esegue la logica del comando."""
        ...


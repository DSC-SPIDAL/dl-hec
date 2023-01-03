from dataclasses import dataclass

@dataclass(init=True)
class Data:
    """Class for keeping track of an item in an inventory."""
    a: str = None
    b: str = "b"
    c: str = "c"

    def print_values(self) -> None:
        print(self.a)
        print(self.b)
        print(self.c)
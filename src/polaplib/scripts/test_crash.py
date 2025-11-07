#!/usr/bin/env python3
# Version: v0.1.0
def f():
    g()


def g():
    raise RuntimeError("demo crash: Python")


if __name__ == "__main__":
    f()

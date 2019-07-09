"""A Mersenne twister implementation.

https://en.wikipedia.org/wiki/Mersenne_Twister
"""

from operator import xor
from typing import List


class SeedException(Exception):
    """Raised on a problem with seeding."""

    pass


class TwisterRand:
    """Set algo constants.

    (w, n, m, r) = (64, 312, 156, 31)
    a = B5026F5AA96619E916
    (u, d) = (29, 555555555555555516)
    (s, b) = (17, 71D67FFFEDA6000016)
    (t, c) = (37, FFF7EEE00000000016)
    l = 43
    """

    n: int = 312
    r: int = 31
    w: int = 64
    m: int = 156
    f: int = 1812433253
    u: int = 29
    d: int = 5555555555555555
    s: int = 17
    b: int = 0x71D67FFFEDA60000
    t: int = 37
    c: int = 0xFFF7EEE000000000
    l: int = 43
    a: int = 0xB5026F5AA96619E9

    def __init__(self):
        """Create a length n array to store the state of the generator. Set Masks."""
        self.MT: List = []
        for i in range(self.n):
            self.MT.append(0)

        self.index: int = self.n + 1
        self.LOWER_MASK: int = (1 << self.r) - 1  # That is, the binary number of r 1's
        self.UPPER_MASK: int = (1 << self.r)

    def seed(self, seed: int):
        """Initialize the generator from a seed."""
        self.index = self.n
        self.MT[0] = seed
        for idx in range(1, self.n - 1):
            self.MT[idx] = (
                self.f * xor(self.MT[idx - 1], (self.MT[idx - 1] >> (self.w - 2))) + idx
            )

    def twist(self):
        """Generate the next n values from the series x_i."""
        for i in range(0, self.n - 1):
            x: int = (self.MT[i] and self.UPPER_MASK) + (
                self.MT[(i + 1) % self.n] and self.LOWER_MASK
            )
            xA: int = x >> 1
            if (x % 2) != 0:
                xA = xor(xA, self.a)
            self.MT[i] = xor(self.MT[(i + self.m) % self.n], xA)
        self.index = 0

    def extract_number(self) -> int:
        """Extract a tempered value based on MT[index] calling twist() every n numbers."""
        if self.index >= self.n:
            if self.index > self.n:
                raise SeedException("Generator was never seeded")
            self.twist()

        y: int = self.MT[self.index]
        y = xor(y, ((y >> self.u) and self.d))
        y = xor(y, ((y << self.s) and self.b))
        y = xor(y, ((y << self.t) and self.c))
        y = xor(y, (y >> self.l))

        self.index += 1
        return y

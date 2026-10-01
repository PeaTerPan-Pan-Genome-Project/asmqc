"""Length statistics shared by modules."""


def nx(lengths: list[int], fraction: float) -> tuple[int, int]:
    """(Nx, Lx): the length and rank at which the cumulative sorted length
    first reaches `fraction` of the total. (0, 0) for no sequence."""
    total = sum(lengths)
    if total == 0:
        return 0, 0
    acc = 0
    for i, n in enumerate(sorted(lengths, reverse=True), 1):
        acc += n
        if acc * 1.0 >= fraction * total:
            return n, i
    raise AssertionError("unreachable")


def pieces_between(length: int, gaps: list[tuple[int, int]]) -> list[int]:
    """Lengths of the pieces of a sequence between gaps (1-based inclusive)."""
    out, pos = [], 0
    for start, end in sorted(gaps):
        if start - 1 > pos:
            out.append(start - 1 - pos)
        pos = max(pos, end)
    if length > pos:
        out.append(length - pos)
    return out

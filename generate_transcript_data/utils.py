from typing import Iterable, Iterator, TypeVar

T = TypeVar("T")


def batch_iterator(iterable: Iterable[T], batch_size: int) -> Iterator[list[T]]:
    """ Yields lists of up to batch_size items, eg for remote APIs that take many IDs per request """
    batch: list[T] = []
    for record in iterable:
        batch.append(record)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

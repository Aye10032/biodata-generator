import math
import random


def sample_positive_normal(
    rng: random.Random,
    mean: float,
    sd: float,
    minimum: int,
    maximum: int,
) -> int:
    if maximum < minimum:
        raise ValueError('maximum is smaller than minimum')
    if sd <= 0:
        return min(max(round(mean), minimum), maximum)
    for _ in range(100):
        value = round(rng.gauss(mean, sd))
        if minimum <= value <= maximum:
            return value
    return min(max(round(mean), minimum), maximum)


def weighted_index(rng: random.Random, weights: list[float]) -> int:
    total = sum(weights)
    if total <= 0:
        raise ValueError('weights must sum to a positive number')
    target = rng.random() * total
    cumulative = 0.0
    for index, weight in enumerate(weights):
        cumulative += weight
        if cumulative >= target:
            return index
    return len(weights) - 1


def zipf_weights(size: int, exponent: float = 1.1) -> list[float]:
    return [1.0 / math.pow(rank, exponent) for rank in range(1, size + 1)]


def random_dna(rng: random.Random, length: int) -> str:
    return ''.join(rng.choice('ACGT') for _ in range(length))

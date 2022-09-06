"""Functions to read and understand YAML configuration files"""
from strictyaml import ScalarValidator, Map, Seq, Str, Optional, Int, Float
from strictyaml.validators import Validator, SeqValidator
from dataclasses import dataclass


def parse_value(x: str, expected_unit: str = None, expected_length: int = -1):
    power = 0
    if expected_unit is None:
        value = x
    else:
        x_split = x.split()
        if len(x_split) > 2:
            # multiple comma-separated values? (hopefully as a literal tuple)
            x_split = ' '.join(x_split[:-1]), x_split[-1]
        if len(x_split) != 2:
            raise RuntimeError(f"{x} does not match '[value] [unit]' expected format")
        value, unit = x_split[0], x_split[1]
        if expected_unit is None:
            expected_unit = unit
        if unit != expected_unit:
            # attempt to support *very limited* power conversions
            pack = {'m': {'cm':  2, 'mm':  3, 'um': 6, 'µm': 6, 'nm': 9, 'å': 10, 'angstrom': 10},
                    'mm': {'m': -3, 'cm': -1, 'um': 3, 'µm': 3, 'nm': 6, 'å': 7, 'angstrom': 7},
                    's': {'ms':  3, 'µs':  6, 'us': 6},
                    'ms': {'s': -3, 'µs':  3, 'us': 3},
                    'µs': {'s': -6, 'ms': -3}
                    }
            if expected_unit not in pack:
                raise RuntimeError("No conversion for {unit} to {expected_unit")
            known = pack[expected_unit]
            if not unit.lower() in known:
                raise RuntimeError("Unknown {unit} for conversion to {expected_unit}")
            power = known[unit.lower()]
    # This is a not-safe thing to do if this code should ever run as part of a service, e.g. on a server
    # A malicious user could create an entry like '__import__("pathlib").Path().absolute()' to do nefarious things
    #
    # Running on your own machine as your own user, you can't do anything special; so just ensure
    # the input you provide doesn't attempt anything bad
    value = eval(value, {})
    if expected_length > -1:
        if isinstance(value, tuple) and len(value) != expected_length:
            raise RuntimeError(f"Expected length {expected_length} item but got {len(value) = }")
        elif not isinstance(value, tuple) and expected_length > 0:
            raise RuntimeError(f"Expected length{expected_length} item but got a scalar")
    power = 10 ** power
    if isinstance(value, tuple):
        value = tuple([v / power for v in value])
    else:
        value /= power
    return value


class List(ScalarValidator):
    def validate_scalar(self, chunk):
        return parse_value(chunk.contents)

    def to_yaml(self, data):
        return f"{data}"


class Val(ScalarValidator):
    def __init__(self, unit: str):
        self.unit = unit

    def validate_scalar(self, chunk):
        return parse_value(chunk.contents, self.unit)

    def to_yaml(self, data):
        return f"{data} {self.unit}"

class PairStrInt(ScalarValidator):
    def __repr__(self):
        return f"Pair(Str, Int)"

    def validate_scalar(self, chunk):
        pair = chunk.contents.split(',')  # do something fancy with parenthetical groups?
        assert len(pair) == 2, "the Pair must be provided two comma separated values"
        return pair[0], int(pair[1])

    def to_yaml(self, data):
        return f"{data[0]}, {data[1]}"


SOURCE_SCHEMA = Map({
    'name': Str(),
    'frequency': Val('Hz'),
    'duration': Val('s'),
    'velocities': Val('m/s'),
    'emission_delay': Val('s'),
})

FREQUENCY_SCHEMA = Map({
    'name': Str(),
    Optional('harmonics'): List(),
    Optional('highest_harmonic'): Int(),
    Optional('ratio'): PairStrInt(),
    Optional('value'): Int(),
})

SEGMENT_SCHEMA = Map({
    'name': Str(), 'length': Val('m'),
    Optional('guide'): Map({'velocities': Val('m/s'), 'short': Val('m'), 'long': Val('m')}),
})

CHOPPER_SCHEMA = Map({
    'name': Str(), 'position': Val('m'), 'opening': Val('degrees'), 'radius': Val('mm'),
    Optional('discs'): Int(), Optional('slots'): Int(),
    Optional('frequency'): Map({'name': Str(), Optional('multiplier'): Int()}),
    'aperture': Map({'width': Val('mm'), 'height': Val('mm'), Optional('offset'): Val('mm')})
})

SAMPLE_SCHEMA = Map({'position': Val('m')})

PRIMARY_SCHEMA = Map({
    'frequencies': Seq(FREQUENCY_SCHEMA),
    'path_segments': Seq(SEGMENT_SCHEMA),
    'choppers': Seq(CHOPPER_SCHEMA),
    'sample': SAMPLE_SCHEMA,
})

SCHEMA = Map({'name': Str(), 'source': SOURCE_SCHEMA, 'primary_spectrometer': PRIMARY_SCHEMA})


def load(filename):
    from pathlib import Path
    from strictyaml import load as syload
    text = Path(filename).read_text()
    yaml = syload(text, SCHEMA)
    return yaml.data


def load_flight_path(path, velocities):
    from .flightpaths import FlightPath, Guide
    length = path['length']
    if 'guide' in path:
        g = path['guide']
        return Guide(name=path['name'], velocity=g['velocities'], shortest=g['short'], longest=g['long'], nominal=(length, length))
    elif 'bragg' in path:
        raise NotImplementedError("Not implemented yet ...")
    else:
        return FlightPath(name=path['name'], velocity=velocities, nominal=(length, length))


def load_chopper(vals, harmonics):
    from numpy import pi, arange
    from .chopper import RectangularAperture, DiscChopper
    phase = 0.
    h = vals['aperture']['height']
    offset = vals['aperture'].get('offset', vals['radius'] - h)
    aperture = RectangularAperture(vals['aperture']['width'], h, offset)
    theta = vals['opening'] / 180 * pi
    slots = vals.get('slots', 1)
    windows = tuple([(z - theta / 2, z + theta / 2) for z in arange(slots) * 2 * pi / slots])

    freq_dict = vals.get('frequency', {})
    phase_to = freq_dict.get('name', 'Source'), freq_dict.get('multiplier', 1)
    frequency = harmonics[phase_to[0]] * phase_to[1]  # don't worry too much about his value yet

    return DiscChopper(name=vals['name'], radius=vals['radius'], frequency=frequency, phase_to=phase_to,
                       phase=phase, aperture=aperture, windows=windows, discs=vals.get('discs', 1))


def prime_factors(n):
    """Returns all the prime factors of a positive integer"""
    factors = []
    d = 2
    while n > 1:
        while n % d == 0:
            factors.append(d)
            n /= d
        d = d + 1
        if d*d > n:
            if n > 1: factors.append(n)
            break
    return factors


def divisors(n):
    return [d for d in range(n+1) if n % d == 0]


def load_frequencies(vals, base_frequency):
    from numpy import arange, ndarray
    from .frequencies import IndependentHarmonics, DependentHarmonics, Frequencies
    # build up the *Harmonics objects before creating the composite Frequency object

    values = {}
    objects = {}
    names = []
    for val in vals:
        name = val.get('name', 'UNDEFINED')
        names.append(name)
        value = val.get('value', 1)
        if 'harmonics' in val:
            harmonics = val['harmonics']
        else:
            highest = val.get('highest_harmonic', 1)
            harmonics = list(arange(highest) + 1)
        if value not in harmonics:
            value = min(harmonics)

        if 'ratio' in val:
            to, ratio = val['ratio']
            assert to in objects, f"the frequency {to} must be defined before {name}"
            allowed = {m: [d for d in divisors(m) if d in harmonics] for m in objects[to].allowed * ratio}
            obj = DependentHarmonics(name, objects[to], allowed)
        else:
            obj = IndependentHarmonics(name, harmonics)

        values[name] = value
        objects[name] = obj

    if 'Source' not in names:
        values['Source'] = 1
        objects['Source'] = IndependentHarmonics('Source', [1])
        names.append('Source')

    harmonics = [values[x] for x in names]
    composite = Frequencies(base_frequency, [objects[x] for x in names], harmonics)

    # use __setitem__ to verify that all harmonics are allowed
    for name, harmonic in zip(names, harmonics):
        composite[name] = harmonic

    return composite


def load_primary_spectrometer(filename):
    from numpy import std, array
    from .primary import PulsedSource, PrimarySpectrometer
    #
    data = load(filename)
    s = data['source']
    delay = s['emission_delay']
    duration = s['duration']
    velocities = s['velocities']
    if not hasattr(duration, '__len__'):
        duration = array((duration,))
    else:
        duration = array(duration)
    if not hasattr(delay, '__len__'):
        delay = array((delay,))
    else:
        delay = array(delay)
    if not hasattr(velocities, '__len__'):
        velocities = array((velocities,))
    else:
        velocities = array(velocities)
    ps = PulsedSource(frequency=s['frequency'], delay=delay, duration=duration, velocities=velocities)
    #
    primary = data['primary_spectrometer']
    #
    frequencies = load_frequencies(primary['frequencies'], ps.frequency)
    #
    paths = primary['path_segments']
    choppers = primary['choppers']
    assert len(paths) == len(choppers) + 1  # [Source] === Chopper === Chopper === Chopper --- [Sample]
    pairs = []
    for path, chopper in zip(paths[:-1], choppers):
        pairs.append((load_flight_path(path, (ps.slowest, ps.fastest)), load_chopper(chopper, frequencies)))
    #
    # There is a key data['primary_spectrometer']['sample'], which has information about the sample
    # but we only need/want the flight path information here
    sample = load_flight_path(paths[-1], (ps.slowest, ps.fastest))
    #
    return PrimarySpectrometer(ps, pairs, sample)
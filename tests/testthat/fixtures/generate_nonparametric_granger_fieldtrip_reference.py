"""Generate and collect the pinned FieldTrip nonparametric-GC reference."""

import json
import platform
import sys

import numpy as np
import scipy
from scipy.io import loadmat, savemat
from scipy.signal.windows import dpss


FIELDTRIP_REVISION = "4c553bda8fb238b0afc549c90388d9c86124d126"


def make_signals():
    rng = np.random.default_rng(20260728)
    n_keep = 4096
    burn = 512
    innovations = rng.standard_normal((n_keep + burn, 2))
    innovations[:, 1] = (
        0.2 * innovations[:, 0]
        + np.sqrt(1 - 0.2**2) * innovations[:, 1]
    )
    x = np.zeros(n_keep + burn)
    y = np.zeros(n_keep + burn)
    for index in range(2, n_keep + burn):
        x[index] = 0.55 * x[index - 1] + innovations[index, 0]
        y[index] = (
            0.35 * y[index - 1]
            + 0.45 * x[index - 2]
            + innovations[index, 1]
        )
    return x[burn:], y[burn:]


def multitaper_csd(x, y, sr, segment_length, n_fft, overlap, nw, n_tapers):
    x = x - np.mean(x)
    y = y - np.mean(y)
    step = max(1, int(np.floor(segment_length * (1 - overlap))))
    starts = np.arange(0, len(x) - segment_length + 1, step)
    tapers = dpss(segment_length, nw, Kmax=n_tapers, sym=True)
    spectrum = np.zeros((2, 2, n_fft), dtype=np.complex128)
    realizations = 0
    for start in starts:
        xs = x[start : start + segment_length]
        ys = y[start : start + segment_length]
        xs = xs - np.mean(xs)
        ys = ys - np.mean(ys)
        for taper in tapers:
            fx = np.fft.fft(xs * taper, n=n_fft)
            fy = np.fft.fft(ys * taper, n=n_fft)
            fourier = np.column_stack((fx, fy))
            normalization = sr * np.sum(taper**2)
            spectrum += np.einsum(
                "fi,fj->ijf", fourier, np.conjugate(fourier)
            ) / normalization
            realizations += 1
    spectrum /= realizations
    one_sided = spectrum[:, :, : n_fft // 2 + 1].copy()
    one_sided[:, :, 1:-1] *= 2
    one_sided = (
        one_sided + np.conjugate(np.swapaxes(one_sided, 0, 1))
    ) / 2
    return one_sided, realizations


def write_input(path):
    settings = {
        "sr": 256.0,
        "segment_length": 256,
        "n_fft": 256,
        "overlap": 0.5,
        "time_bandwidth": 3.0,
        "n_tapers": 5,
    }
    x, y = make_signals()
    spectrum, realizations = multitaper_csd(
        x,
        y,
        settings["sr"],
        settings["segment_length"],
        settings["n_fft"],
        settings["overlap"],
        settings["time_bandwidth"],
        settings["n_tapers"],
    )
    frequencies = np.arange(settings["n_fft"] // 2 + 1) * (
        settings["sr"] / settings["n_fft"]
    )
    savemat(
        path,
        {
            "x": x,
            "y": y,
            "S": spectrum,
            "freq": frequencies,
            "sr": settings["sr"],
            "segment_length": settings["segment_length"],
            "n_fft": settings["n_fft"],
            "overlap": settings["overlap"],
            "time_bandwidth": settings["time_bandwidth"],
            "n_tapers": settings["n_tapers"],
            "n_realizations": realizations,
        },
        do_compression=True,
    )


def nested_real(array):
    return np.real(np.asarray(array)).tolist()


def nested_imag(array):
    return np.imag(np.asarray(array)).tolist()


def collect(input_path, output_path, json_path):
    inputs = loadmat(input_path, squeeze_me=True)
    outputs = loadmat(output_path, squeeze_me=True)
    payload = {
        "x": np.asarray(inputs["x"]).reshape(-1).tolist(),
        "y": np.asarray(inputs["y"]).reshape(-1).tolist(),
        "sr": float(inputs["sr"]),
        "segment_length": int(inputs["segment_length"]),
        "n_fft": int(inputs["n_fft"]),
        "overlap": float(inputs["overlap"]),
        "time_bandwidth": float(inputs["time_bandwidth"]),
        "n_tapers": int(inputs["n_tapers"]),
        "n_realizations": int(inputs["n_realizations"]),
        "frequencies": np.asarray(inputs["freq"]).reshape(-1).tolist(),
        "csd_real": nested_real(inputs["S"]),
        "csd_imag": nested_imag(inputs["S"]),
        "gc_xy": np.asarray(outputs["gc_xy"]).reshape(-1).tolist(),
        "gc_yx": np.asarray(outputs["gc_yx"]).reshape(-1).tolist(),
        "transfer_real": nested_real(outputs["H"]),
        "transfer_imag": nested_imag(outputs["H"]),
        "sigma": nested_real(outputs["Z"]),
        "reconstructed_real": nested_real(outputs["Sf"]),
        "reconstructed_imag": nested_imag(outputs["Sf"]),
        "provenance": {
            "reference": (
                "FieldTrip sfactorization_wilson and "
                "ft_connectivity_granger"
            ),
            "fieldtrip_revision": FIELDTRIP_REVISION,
            "fieldtrip_commit_url": (
                "https://github.com/fieldtrip/fieldtrip/commit/"
                + FIELDTRIP_REVISION
            ),
            "sfactor_source_url": (
                "https://github.com/fieldtrip/fieldtrip/blob/"
                + FIELDTRIP_REVISION
                + "/connectivity/private/sfactorization_wilson.m"
            ),
            "granger_source_url": (
                "https://github.com/fieldtrip/fieldtrip/blob/"
                + FIELDTRIP_REVISION
                + "/connectivity/ft_connectivity_granger.m"
            ),
            "octave_version": str(outputs["octave_version"]),
            "python_version": platform.python_version(),
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "generation_date": "2026-07-28",
            "commands": [
                (
                    "python3 "
                    "generate_nonparametric_granger_fieldtrip_reference.py "
                    "input /tmp/wscb11-input.mat"
                ),
                (
                    "docker run --rm -v /tmp:/tmp -v $PWD:/fixture "
                    "-v /tmp/wscb11-fieldtrip:/fieldtrip "
                    "gnuoctave/octave:9.4.0 octave --quiet "
                    "/fixture/generate_nonparametric_granger_fieldtrip_"
                    "reference.m /tmp/wscb11-input.mat "
                    "/tmp/wscb11-output.mat /fieldtrip"
                ),
                (
                    "python3 "
                    "generate_nonparametric_granger_fieldtrip_reference.py "
                    "collect /tmp/wscb11-input.mat /tmp/wscb11-output.mat "
                    "/tmp/wscb11-reference.json"
                ),
            ],
        },
    }
    with open(json_path, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)


def main(argv):
    if len(argv) == 3 and argv[1] == "input":
        write_input(argv[2])
        return
    if len(argv) == 5 and argv[1] == "collect":
        collect(argv[2], argv[3], argv[4])
        return
    raise SystemExit(
        f"usage: {argv[0]} input INPUT.mat | "
        "collect INPUT.mat OUTPUT.mat OUTPUT.json"
    )


if __name__ == "__main__":
    main(sys.argv)

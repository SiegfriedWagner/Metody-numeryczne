import numpy as np
import argparse
from os.path import exists
from scipy.io.wavfile import read as readwav
from scipy.io.wavfile import write as writewav
from scipy.sparse import save_npz, load_npz, coo_matrix
from pathlib import Path
from haar import haar, recreate_signal
import sys


def cutout_haar_threshold(h, band, percentage=False):
    """
    Removes data that are not with out (-band, band) or (-max(data) * band, max(data) * band) if percentage is True
    for each data in haar result.
    """
    orig_band = band
    for index in range(len(h)):
        arr = h[index]
        if percentage:
            band = np.max(arr) * orig_band  # TODO: Better implementation, maybe with gaussian threshold
        logical = np.logical_not(np.logical_or((arr > band), (arr < -band)))
        print(f"Zeroed: {logical.size - np.sum(logical)} values")
        arr[logical] = 0
        h[index] = arr


def save_compressed_haar(dirname, haar_result, orginal_signal_size: int, original_sampling_rate: int):
    """
    Custom serializer that stores data as sparse matrices without zeros
    """
    directory = Path(dirname)
    if not directory.exists():
        directory.mkdir()
    with (directory / "info.txt").open('w') as f:
        f.write(str(orginal_signal_size))
        f.write("\n")
        f.write(str(original_sampling_rate))
    sparse = coo_matrix(haar_result[0])
    sparse.eliminate_zeros()
    save_npz(directory / "a.npz", sparse, 'coo')
    for index, detail in enumerate(haar_result[1:]):
        sparse = coo_matrix(detail)
        sparse.eliminate_zeros()
        save_npz(directory / f"d{index}.npz", sparse, 'coo')


def load_compressed_haar(dirname):
    """
    Custom deserializer that loads haar result stored as sparse matrices
    """
    directory = Path(dirname)
    with (directory / "info.txt").open('r') as f:
        original_size = int(f.readline())
        original_sampling_rate = int(f.readline())
    a = load_npz(directory / "a.npz").toarray()
    a = a.reshape(a.size)

    def load_detail(path):
        d = load_npz(path).toarray()
        return d.reshape(d.size)

    h = (a, *(load_detail(path) for path in directory.glob("d*.npz")))
    return h, original_size, original_sampling_rate


if __name__ == "__main__":
    parser = argparse.ArgumentParser(sys.argv[0])
    parser.description = "Compressed wav file"
    parser.add_argument("uncompressed_file", help="Path to uncompressed file [wav music]")
    parser.add_argument("compressed_directory", help="Path to compress directory")
    parser.add_argument("-x", "--extract", action="store_true", help="Extract compressed file")
    parser.add_argument("-o", "--overwrite", action="store_true", help="Allows overwriting file")
    parser.add_argument("-c", "--cutout", type=float, default=0.05, help="Cutout filter [0.0 - 1.0] (default: 0.05)")
    parser.add_argument("-n", "--haar", default=5, type=int, help="Number of haar wavelets to generate (default: 5)")
    args = parser.parse_args()
    if not args.extract:
        # compress
        if not args.overwrite and exists(args.compressed_directory):
            print("Compressed archive exists")
            exit(1)
        sampling_rate, data = readwav("guitar.wav")
        print(f"Read file sampling rate: {sampling_rate}, samples: {data.size}")
        l = int(np.log2(data.size))
        print(f"Cutting signal from {data.size} to {pow(2, l)}")
        data = data[:pow(2, l)]
        print("Calculating haar wave compression")
        h = haar(data, args.haar)
        cutout_haar_threshold(h, args.cutout, True)
        save_compressed_haar(args.compressed_directory, h, data.size, sampling_rate)
        print(f"Saved compressed signal to: {args.compressed_directory}")
    else:
        if not args.overwrite and exists(args.uncompressed_file):
            print("Uncompressed file exists")
            exit(1)
        print("Reading compressed archive")
        h, original_signal_length, sampling_rate = load_compressed_haar(args.compressed_directory)
        print(f"Signal samples: {original_signal_length}")
        print(f"Signal sampling rate: {sampling_rate}")
        print("Recreating signal from haar waves")
        signal = recreate_signal(h[0], h[1:], original_signal_length)
        writewav(args.uncompressed_file, sampling_rate, signal.astype(np.int16))
        print(f"Saved recreated signal to: {args.uncompressed_file}")

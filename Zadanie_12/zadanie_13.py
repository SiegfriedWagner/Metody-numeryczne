import argparse
import sys
import numpy as np
from PIL import Image
from haar import haar2d, recreate_signal2d
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = argparse.ArgumentParser(sys.argv[0])
    parser.description = "Compressed image"
    parser.add_argument("input_image", help="Path to input image")
    args = parser.parse_args()
    image = Image.open(args.input_image)
    data = np.asarray(image).astype(np.float)
    f1 = plt.figure(1)
    plt.imshow(data, cmap="gray", vmin=0, vmax=np.max(data))
    plt.title("Oryginalny obrazek")
    haar1 = haar2d(data, 1)
    f2 = plt.figure(2)
    plt.imshow(haar1, cmap="gray", vmin=0, vmax=np.max(haar1))
    plt.title("Haar 1")
    haar2 = haar2d(data, 2)
    f3 = plt.figure(3)
    plt.imshow(haar2, cmap="gray", vmin=0, vmax=np.max(haar2))
    plt.title("Haar 2")
    haar3 = haar2d(data, 3)
    f4 = plt.figure(4)
    plt.imshow(haar3, cmap="gray", vmin=0, vmax=np.max(haar3))
    plt.title("Haar 3")
    recreated = recreate_signal2d(haar1, 1)
    f5 = plt.figure(5)
    plt.imshow(recreated, cmap="gray", vmin=0, vmax=np.max(recreated))
    plt.title("Odtworzony obraz haar 1")
    recreated = recreate_signal2d(haar2, 2)
    f6 = plt.figure(6)
    plt.imshow(recreated, cmap="gray", vmin=0, vmax=np.max(recreated))
    plt.title("Odtworzony obraz haar 2")
    recreated = recreate_signal2d(haar3, 3)
    f7 = plt.figure(7)
    plt.imshow(recreated, cmap="gray", vmin=0, vmax=np.max(recreated))
    plt.title("Odtworzony obraz haar 3")
    plt.show()

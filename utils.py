import matplotlib.pyplot as plt

def savefig_bioinformatics(path_save):
    """Save image for publication in bioinformatics

    Args:
        path_save (os.Path): Filename without any extension.
    """

    # low resolution jpg
    path_save_jpeg = path_save + ".jpg"
    plt.savefig(path_save_jpeg, format="jpg")

    # high resolution eps, to be converted to tif
    path_save_eps = path_save + ".eps"
    plt.savefig(path_save_eps, format="eps", dpi=1200)

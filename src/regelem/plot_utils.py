import subprocess
import sys

def view_plot(filepath : str):
    """Attempts to opens an image file (.jpg, .png) of a plot in the system's default image viewer.

    Args:
        filepath (str): Filepath to the image file
    """
    print("Attempting to open in file viewer...");
    imageViewerFromCommandLine = {'linux':'xdg-open',
                                'win32':'explorer',
                                'darwin':'open'}[sys.platform]
    subprocess.run([imageViewerFromCommandLine, filepath])


def show_plot(filepath: str):
    """Alias for view_plot()
    """

    view_plot(filepath)
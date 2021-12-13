# InfoVis-FinalProject
Code repo for information vidualization final project

The main denoising function used in this project is mpdenoise.py
I include a demo notebook showing an example of how this denoising can be applied to diffusion data, this notebook requires the python library dipy in order to download example diffusion images. The file demo_tracts.py performs both denoising as well as some light visualization of fiber tracking, however tracking results cannot be visualized in a notebook context.
Code used for data processing (registration, segmentation, artifact correction, figure generation) can be found in the "processing" folder.

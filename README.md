# visiona
visiona is an optical coded target tracking library. It allows to precisely and accurately determine the relative position and orientation of a target with respect to a calibrated camera. 

![concentric circle coded target used by visiona](https://github.com/DavideACucci/visiona/blob/master/target/target_raster.png)

Some of the principles are described in details in the following publication, along with test results:

    @article{cucci2016accurate,
      AUTHOR = {Cucci, D. A.},
      TITLE = {ACCURATE OPTICAL TARGET POSE DETERMINATION FOR APPLICATIONS IN AERIAL PHOTOGRAMMETRY},
      JOURNAL = {ISPRS Annals of Photogrammetry, Remote Sensing and Spatial Information Sciences},
      VOLUME = {III-3},
      YEAR = {2016},
      PAGES = {257--262},
    }

## Current Limitations

The current implementation searches for a specific target whose code is known a-priori. The library is ready to detect and measure multiple different targets in the same image, but minor modifications to the code are required.

## Licence

The source code is released under a GPLv3 licence.

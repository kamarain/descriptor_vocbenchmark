  VOC Benchmarking Toolbox of Local Region Detectors and Descriptors
  ------------------------------------------------------------------

This toolbox contains Matlab functions to evaluate repeatability and
matching of local region detectors and descriptors over multiple
images of several visual object classes (faces, cars, motorbikes etc.)

The first results were published in [1] and the project continues by
the researchers in Lappeenranta University of Technology, Tampere
University of Technology and University of Southern Denmark.

For more information, read the Wiki page at:
https://bitbucket.org/kamarain/descriptor_vocbenchmark/wiki/Home

References:

[1] Lankinen, J.; Kangas, V.; Kamarainen, J.-K.. A comparison of local
feature detectors and descriptors for visual object categorization by
intra-class repeatability and matching, In 21th Int. Conf. on Pattern
Recognition (ICPR2012), Tsukuba Science City, Japan, 2012.

*Bugs*

None reported.


*Todo*

1. Unifying detectors_benchmark and descriptors_benchmark and, in
particular, Mikolajczyk's repeatability.m (our mikolajczyk_repeatability.m
and mikolajczyk_repeatability2.m) since the both pairs basically do the
same thing - having two different functions is ridiculous (and still
there because of historic reasons).


# Matlab Toolbox for Benchmarking Local Region Detectors and Descriptors

The main idea for this Toolbox derives from the two papers by Krystian Mikolajczyk:

* [A Comparison of Affine Region Detectors](http://www.robots.ox.ac.uk/~vgg/research/affine/det_eval_files/vibes_ijcv2004.pdf), International Journal of Computer Vision, 2005
* [A Performance Evaluation of Local Descriptors](http://research.microsoft.com/en-us/um/people/manik/projects/trade-off/papers/MikolajczykPAMI05.pdf), IEEE Transactions on Pattern Analysis and Machine Intelligence, 2005.

The main difference in our work is that while Mikolajczyk used images of same scenes under various geometric transformations and imaging distortions we use multiple instances of visual object classes (motorbikes, cars, human faces, etc.) Our work is motivated by the fact that local region detectors (such as [SIFT](http://en.wikipedia.org/wiki/Scale-invariant_feature_transform)) are popular tools in learning and classifying objects and object classes in images. In the context of visual object classification, our benchmark can be used to measure whether a method X finds regions of similar shape and location and a method Y can match them well over different instances of classes.

## Installation

First you need *descriptor_vocbenchmark* functionality from Bitbucket:
```
$ hg clone https://bitbucket.org/kamarain/descriptor_vocbenchmark
```

The code needs many helper functions available in the *mvprmatlab* repository
```
$ hg clone https://bitbucket.org/kamarain/mvprmatlab
```
which needs to be in Matlab path
```
#!matlab

>> addpath <MYDIRS>/mvprmatlab
```

To run the basic detectors (wrapped by mvpr_feature_extract_new.m) you need to download certain files. First you can download VLFeat toolbox from [www.vlfeat.org](www.vlfeat.org) and extract the tar ball somewhere and add to the Matlab path:
```
#!matlab

>> addpath <MYDIRS>/vlfeat-<VERSIO>/toolbox/
>> vl_setup
```
To use feature space detectors and descriptors by Mikolajczyk, download his [implementation](http://kahlan.eps.surrey.ac.uk/featurespace/web/) and extract both files (64-bit and 32-bit version) into path *~/<MY_PATH>/mvprmatlab/binaries*. Remember to check that there are correct paths for binaries in the file mvpr_feature_extract_new.m.

To run OpenCV implementation, you should have the latest version of it (at least 2.4.3). If you already have installed OpenCV successfully you can skip the following step:

```
$ cd ~/<MY_WORKING DIR>
$ git clone https://github.com/Itseez/opencv.git
$ cd opencv
$ mkdir build
$ cd build
$ cmake ..
$ make
```

After we have installed necessary OpenCV libraries, we can start to compile our descriptor's source code. Instructions are given below:   
```
$ cd ~/<MY_PATH>/mvprmatlab/opencv_descriptors/
$ mkdir build; cd build
$ cmake -DOpenCV_DIR=<OPENCV_EXTRACTION_DIR>/build/ ..
$ make
$ cp bin/opencv_descriptors ~/<MY_PATH>/mvprmatlab/binaries/
```


There is one more annoying thing to install, that is Oliver Woodford's great [export_fig](http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig) at the Matlab Central file exchange. That is used to save graphs in nice formats. Save and extract the zip somewhere and add to the Matlab path:
```
#!matlab

>> addpath <MYDIRS>/export_fig/
```
 This Toolbox calls also one C function which has to be compiled before we can run the program:
```
#!matlab

>> cd ~/<MY_PATH>/matlab/base 
>> mex c_eoverlap.cxx
```

## Benchmarking Detectors

To run evaluation experiments, you need a) images of objects and b) ground truth landmarks. Next, we demonstrate the workflow by giving a step by step example how to run the "debug" configuration settings
based basic experiment.

First you need to download the [Caltech-101](http://www.vision.caltech.edu/Image_Datasets/Caltech101/) images, outlines of the objects in the these images (Annotations.tar) and optional m-file (show_annotation.m
) which is a matlab script to view the annotations:
```
$ cd /home/<ME>/<MY_DATA_DIR>
$ mkdir caltech101
$ cd caltech101
$ wget http://www.vision.caltech.edu/Image_Datasets/Caltech101/101_ObjectCategories.tar.gz
$ tar zxfv 101_ObjectCategories.tar.gz
$ wget http://www.vision.caltech.edu/Image_Datasets/Caltech101/Annotations.tar
$ tar -xvf Annotations.tar
$ wget http://www.vision.caltech.edu/Image_Datasets/Caltech101/show_annotation.m
```
We have also made available manually annotated landmarks which are used in method evaluation: [landmarks.tar.gz](https://bitbucket.org/kamarain/descriptor_vocbenchmark/downloads/landmarks.tar.gz) - download and extract to the *caltech101* directory. You can also download two more challenging data set, [Randomized Caltech-101 Data Set](https://bitbucket.org/kamarain/descriptor_vocbenchmark/downloads/r-caltech101.tar.gz) and [Imagenet Data Set](https://bitbucket.org/kamarain/descriptor_vocbenchmark/downloads/Imagenet.tar.gz). Extract the content into your working directory.

Since you have your stuff in various directories, you need to edit the configuration file:
```
$ cp detectors_benchmark1_conf.m.template detectors_benchmark1_conf.m
$ <MYEDITOR> detectors_benchmark1_conf.m 
```
Add also the line conf.detectorOverlap = [40] (or value 50) to section *report2* in the file detectors_benchmark1_conf.m. When the paths etc. are correct (see Installation instruction), then in Matlab prompt just type
```
#!matlab

>> cd ~/<MY_PATH>/matlab
>> addpath base
>> detectors_benchmark1
```
Some notes before we start running the toolbox: The Annotation.tar contains incorrectly named folder names. Atleast you should rename the folders Faces_3, Airplanes_Side_2 and Motorbikes_16 to Faces_easy, airplanes and Motorbikes respectively. 

If everything goes nicely (see the output to hunt possible bugs and errors), then you should end up with a result graph looking somewhat in figure below. The result is obtained using Feature Space Implementaitions (fs_hesaff+fs_sift).

![Detector evaluation result (config:debug)](http://bitbucket.org/kamarain/descriptor_vocbenchmark/downloads/detector_evaluation_debug_plot.png)


## Benchmarking Descriptors

To run descriptros benchmark, we have to do a few similar step as we did with the detector section.
First make your own copy of the configuration file and then add the correct paths:

```
$ descriptors_benchmark1_conf.m.template descriptors_benchmark1_conf.m
$ <MYEDITOR> descriptors_benchmark1_conf.m
```
If you were already able to run detector benchmark successfully, then everything should be made-up and you can just type in Matlab prompt
```
#!matlab
>> descriptors_benchmark1
```
If Matlab gives an error, make sure you have done all the necessaries described in the detector section.
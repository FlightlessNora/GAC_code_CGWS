Title:      Matlab code for "Efficient Non-uniform Deblurring Based on Generalized Additive Convolution Model"  
Author:     Hong Deng <denghong_hit@163.com>  
Version:    1.0
Date:       May 15, 2017  
Copyright:  2013, Hong Deng

Matlab code for "Efficient Non-uniform Deblurring Based on Generalized Additive Convolution Model"
==========================================================

This package contains code to deblur nonuniform blurry image, using the model described in
the paper

Hong Deng, Dongwei Ren, David Zhang, Wangmeng Zuo, Hongzhi Zhang and Kuanquan Wang
"Efficient Non-uniform Deblurring Based on Generalized Additive Convolution Model". In EURASIP Journal on Advances in Signal Processing, 2016.

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) 
in an academic publication.

For algorithmic details, please refer to our paper.

----------------
How to use
----------------
Most of the code is in Matlab, while the code in this package makes use of QPC package (by Adrian Wills), Multivariate Gaussian Distribution (by Timothy Felty) and
bilateral filter (by Douglas R. Lanman),
which are subject to their own license. For functions 
originally distributed by the toolbox, please refer to the original readme from 
the toolbox for details.
I have test the code in Matlab 2014a under Windows. It should 
work under Matlab.

1. include code/subdirectory in your Matlab path
2. run Deblur_demo_CGWS_func(maxit,innerit,im_id,win_num) or GAC_demo.m

For the test of any image, you should prepare the estimated blur kernel from several regions of the input blurry image (save in "k_set") and their position (save in "loc"). 
Then you should be able to remove non-uniform blur as in the demo.

----------------
User specified parameter
----------------
There are only a few parameters need to be specified by users.

'win_num' - construct the file name of the .mat file, which decides the size of image patch used to estimate the blur kernel.
'im_id' - number in the image name you want to deblur.

----------------
IMPORTANT NOTE 
----------------
Usage of the executable file is permitted only for academic purpose. 
If you find any bug, please send report to <denghong_hit@163.com> 


## 3D facial landmark localization

###Instructions

* ./files include *.pts files and *.abs files. Please refer to meshSpin.m for detailed info.

* ./spin_imgs include spin image feature of the samples.

* main.m contains steps of the algorithm.

* calcSpinImages.m calculates the spin image feature of each vertex on input mesh.

* extract_abs.m is used to read info from *.abs file, including flag and 3D-coordinate.

* read_shape.m is used to read 2D-coordinate from *.pts file.

* meshSpin.m obtains the spin image feature and 3D-coordinate of input mesh, the spin image and 3D-coordinate of keypoints. Note: the *.pts parameter should be NULL for test samples, as their locations are unknown.

### Examples
Here are images including 5 keypoints and 68 keypoints, respectively.

![image of 5pt](https://qizAidon.github.com/3d_facial/5pt.png)

![image of 5pt](https://qizAidon.github.com/3d_facial/68pt.png)